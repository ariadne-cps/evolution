/***************************************************************************
 *            hybrid_evolution_tree.cpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Evolution.
 *
 *  Ariadne Evolution is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Evolution is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Evolution. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file src/hybrid_evolution_tree.cpp
 *  @brief Parent-first trajectory scheduling and history-based Orbit projection.
 */

#include "evolution/hybrid_evolution_tree.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <utility>

#include "evolution/enclosure_recombiner.hpp"
#include "function/formula.hpp"
#include "function/formula.tpl.hpp"
#include "logging/logging.hpp"
#include "threading/thread_manager.hpp"
#include "threading/workload.hpp"

namespace Ariadne {

OnParentCompletionRecombinationTriggerPolicy* OnParentCompletionRecombinationTriggerPolicy::clone() const {
    return new OnParentCompletionRecombinationTriggerPolicy(*this);
}

bool OnParentCompletionRecombinationTriggerPolicy::should_evaluate(
    const HybridEvolutionTree&, NodeId, HybridRecombinationEvaluationPoint point, std::size_t) const {
    return point == HybridRecombinationEvaluationPoint::PARENT_COMPLETE;
}

AfterEachStepRecombinationTriggerPolicy* AfterEachStepRecombinationTriggerPolicy::clone() const {
    return new AfterEachStepRecombinationTriggerPolicy(*this);
}

bool AfterEachStepRecombinationTriggerPolicy::should_evaluate(
    const HybridEvolutionTree&, NodeId, HybridRecombinationEvaluationPoint point, std::size_t) const {
    return point == HybridRecombinationEvaluationPoint::AFTER_STEP;
}

OnSynchronizationCompletionRecombinationTriggerPolicy*
    OnSynchronizationCompletionRecombinationTriggerPolicy::clone() const {
    return new OnSynchronizationCompletionRecombinationTriggerPolicy(*this);
}

bool OnSynchronizationCompletionRecombinationTriggerPolicy::should_evaluate(
    const HybridEvolutionTree&, NodeId, HybridRecombinationEvaluationPoint point, std::size_t) const {
    return point == HybridRecombinationEvaluationPoint::SYNCHRONIZATION_COMPLETE;
}

EveryNStepsRecombinationTriggerPolicy::EveryNStepsRecombinationTriggerPolicy(std::size_t period)
    : _period(period) {
    if (_period == 0u)
        throw std::invalid_argument("Recombination evaluation period must be positive");
}

EveryNStepsRecombinationTriggerPolicy* EveryNStepsRecombinationTriggerPolicy::clone() const {
    return new EveryNStepsRecombinationTriggerPolicy(*this);
}

bool EveryNStepsRecombinationTriggerPolicy::should_evaluate(
    const HybridEvolutionTree&, NodeId, HybridRecombinationEvaluationPoint point,
    std::size_t completed_steps) const {
    if (point == HybridRecombinationEvaluationPoint::AFTER_STEP)
        return completed_steps % _period == 0u;
    return (point == HybridRecombinationEvaluationPoint::PARENT_COMPLETE ||
            point == HybridRecombinationEvaluationPoint::SYNCHRONIZATION_COMPLETE) &&
           completed_steps % _period != 0u;
}

NoRecombinationSchedulingPolicy* NoRecombinationSchedulingPolicy::clone() const {
    return new NoRecombinationSchedulingPolicy(*this);
}

std::vector<std::vector<HybridEvolutionTree::NodeId>>
    NoRecombinationSchedulingPolicy::select(const HybridEvolutionTree&, NodeId) const {
    return {};
}

SiblingWindowRecombinationSchedulingPolicy::SiblingWindowRecombinationSchedulingPolicy(
    std::size_t max_group_size, double max_time_spread)
    : _max_group_size(max_group_size), _max_time_spread(max_time_spread) {
    if (_max_group_size < 2u || !std::isfinite(_max_time_spread) || _max_time_spread < 0.0)
        throw std::invalid_argument("Sibling-window recombination requires group size >= 2 and a finite nonnegative time spread");
}

SiblingWindowRecombinationSchedulingPolicy* SiblingWindowRecombinationSchedulingPolicy::clone() const {
    return new SiblingWindowRecombinationSchedulingPolicy(*this);
}

namespace {

// Radius diagnostics are evaluated only inside enabled logging macros.
double largest_source_radius(const List<Enclosure>& sources) {
    double radius = 0.0;
    for (const auto& source: sources)
        radius = std::max(radius, source.radius().raw().get_d());
    return radius;
}

double largest_output_radius(const EnclosureRecombinationResult& result) {
    double radius = 0.0;
    for (const auto& output: result.outputs)
        radius = std::max(radius, output.enclosure.radius().raw().get_d());
    return radius;
}

// Require the same transition family and identical symbolic auxiliary maps.
// This check is shared by the built-in policy and the backend's policy validation.
bool compatible_recombination_children(const HybridEvolutionTree& tree,
                                       HybridEvolutionTree::NodeId parent,
                                       HybridEvolutionTree::NodeId first,
                                       HybridEvolutionTree::NodeId next) {
    const auto& first_metadata = tree.metadata(first);
    const auto& next_metadata = tree.metadata(next);
    if (!first_metadata.origin || !next_metadata.origin)
        return false;
    const auto& parent_entries = tree.object(parent).entries;
    const auto& first_successor = parent_entries.at(first_metadata.origin->entry).result.value().successors.at(first_metadata.origin->successor);
    const auto& next_successor = parent_entries.at(next_metadata.origin->entry).result.value().successors.at(next_metadata.origin->successor);
    const auto& a = tree.object(first).entries.front().work;
    const auto& b = tree.object(next).entries.front().work;
    return first_successor.kind == HybridEvolutionSuccessorKind::TRANSITION &&
           next_successor.kind == HybridEvolutionSuccessorKind::TRANSITION &&
           first_successor.event == next_successor.event &&
           first_successor.event_kind == next_successor.event_kind &&
           a.phase == b.phase &&
           a.enclosure.location() == b.enclosure.location() &&
           a.enclosure.state_space() == b.enclosure.state_space() &&
           a.enclosure.auxiliary_space() == b.enclosure.auxiliary_space() &&
           a.enclosure.previous_events() == b.enclosure.previous_events() &&
           identical(make_formula(a.enclosure.continuous_set().auxiliary_mapping()),
                     make_formula(b.enclosure.continuous_set().auxiliary_mapping()));
}

// Same numerical and discrete meaning is required for future work. A
// validated recombiner combines the clock as part of the continuous set;
// overlapping time bounds are an additional conservative scheduler filter.
bool compatible_frontier_checkpoints(const HybridEvolutionTree& tree,
                                     HybridEvolutionCheckpointId first,
                                     HybridEvolutionCheckpointId next,
                                     bool require_time_overlap = true) {
    const auto& a = tree.object(first.trajectory).entries.at(first.entry).work;
    const auto& b = tree.object(next.trajectory).entries.at(next.entry).work;
    const auto ta = a.enclosure.time_range();
    const auto tb = b.enclosure.time_range();
    if (require_time_overlap &&
        (ta.lower_bound().raw().get_d() > tb.upper_bound().raw().get_d() ||
         tb.lower_bound().raw().get_d() > ta.upper_bound().raw().get_d()))
        return false;
    return a.phase == b.phase &&
           a.enclosure.location() == b.enclosure.location() &&
           a.enclosure.state_space() == b.enclosure.state_space() &&
           a.enclosure.auxiliary_space() == b.enclosure.auxiliary_space() &&
           a.enclosure.previous_events() == b.enclosure.previous_events() &&
           identical(make_formula(a.enclosure.continuous_set().auxiliary_mapping()),
                     make_formula(b.enclosure.continuous_set().auxiliary_mapping()));
}

} // namespace

AllCompatibleCheckpointRecombinationPolicy*
    AllCompatibleCheckpointRecombinationPolicy::clone() const {
    return new AllCompatibleCheckpointRecombinationPolicy(*this);
}

std::vector<std::vector<HybridEvolutionCheckpointId>>
    AllCompatibleCheckpointRecombinationPolicy::select(
        const HybridEvolutionTree& tree, const std::vector<CheckpointId>& pending) const {
    std::vector<std::vector<CheckpointId>> groups;
    for (const auto checkpoint: pending) {
        const auto group = std::find_if(groups.begin(), groups.end(), [&](const auto& candidates) {
            return compatible_frontier_checkpoints(tree, candidates.front(), checkpoint, false);
        });
        if (group == groups.end())
            groups.push_back({checkpoint});
        else
            group->push_back(checkpoint);
    }
    groups.erase(std::remove_if(groups.begin(), groups.end(), [](const auto& group) {
                     return group.size() < 2u;
                 }),
                 groups.end());
    return groups;
}

ImmediateCheckpointAdmissionPolicy* ImmediateCheckpointAdmissionPolicy::clone() const {
    return new ImmediateCheckpointAdmissionPolicy(*this);
}

bool ImmediateCheckpointAdmissionPolicy::admit_successor(
    const HybridEvolutionWorkItem&, const HybridEvolutionSuccessor&) const { return true; }

std::vector<HybridEvolutionCheckpointId> ImmediateCheckpointAdmissionPolicy::select_next(
    const HybridEvolutionTree&, const std::vector<CheckpointId>& waiting) const { return waiting; }

LocationCheckpointAdmissionPolicy* LocationCheckpointAdmissionPolicy::clone() const {
    return new LocationCheckpointAdmissionPolicy(*this);
}

bool LocationCheckpointAdmissionPolicy::admit_successor(
    const HybridEvolutionWorkItem& producer, const HybridEvolutionSuccessor& successor) const {
    return successor.kind != HybridEvolutionSuccessorKind::TRANSITION &&
           successor.work.enclosure.location() == producer.enclosure.location();
}

std::vector<HybridEvolutionCheckpointId> LocationCheckpointAdmissionPolicy::select_next(
    const HybridEvolutionTree& tree, const std::vector<CheckpointId>& waiting) const {
    std::vector<CheckpointId> selected;
    if (waiting.empty())
        return selected;
    const auto location = tree.object(waiting.front().trajectory).entries.at(waiting.front().entry).work.enclosure.location();
    for (const auto checkpoint: waiting)
        if (tree.object(checkpoint.trajectory).entries.at(checkpoint.entry).work.enclosure.location() == location)
            selected.push_back(checkpoint);
    return selected;
}

NoFrontierRecombinationPolicy* NoFrontierRecombinationPolicy::clone() const {
    return new NoFrontierRecombinationPolicy(*this);
}

std::vector<std::vector<HybridEvolutionCheckpointId>>
    NoFrontierRecombinationPolicy::select(
        const HybridEvolutionTree&, const std::vector<CheckpointId>&) const {
    return {};
}

FrontierWindowRecombinationPolicy::FrontierWindowRecombinationPolicy(
    std::size_t max_group_size, double max_time_spread)
    : _max_group_size(max_group_size), _max_time_spread(max_time_spread) {
    if (_max_group_size < 2u || !std::isfinite(_max_time_spread) || _max_time_spread < 0.0)
        throw std::invalid_argument("Frontier-window recombination requires group size >= 2 and finite nonnegative time spread");
}

FrontierWindowRecombinationPolicy* FrontierWindowRecombinationPolicy::clone() const {
    return new FrontierWindowRecombinationPolicy(*this);
}

std::vector<std::vector<HybridEvolutionCheckpointId>>
    FrontierWindowRecombinationPolicy::select(
        const HybridEvolutionTree& tree, const std::vector<CheckpointId>& pending) const {
    std::vector<std::vector<CheckpointId>> groups;
    std::vector<bool> used(pending.size(), false);
    for (std::size_t first = 0; first != pending.size(); ++first) {
        if (used[first])
            continue;
        const auto& initial = tree.object(pending[first].trajectory).entries.at(pending[first].entry).work.enclosure;
        const auto interval = initial.time_range();
        double lower = interval.lower_bound().raw().get_d();
        double upper = interval.upper_bound().raw().get_d();
        if (!std::isfinite(lower) || !std::isfinite(upper) || upper - lower > _max_time_spread)
            continue;
        std::vector<CheckpointId> group{pending[first]};
        for (std::size_t next = first + 1; next != pending.size() && group.size() < _max_group_size; ++next) {
            if (used[next] || !compatible_frontier_checkpoints(tree, pending[first], pending[next]))
                continue;
            const auto time = tree.object(pending[next].trajectory).entries.at(pending[next].entry).work.enclosure.time_range();
            const double candidate_lower = time.lower_bound().raw().get_d();
            const double candidate_upper = time.upper_bound().raw().get_d();
            if (!std::isfinite(candidate_lower) || !std::isfinite(candidate_upper) ||
                std::max(upper, candidate_upper) - std::min(lower, candidate_lower) > _max_time_spread)
                continue;
            lower = std::min(lower, candidate_lower);
            upper = std::max(upper, candidate_upper);
            group.push_back(pending[next]);
            used[next] = true;
        }
        if (group.size() >= 2u) {
            used[first] = true;
            groups.push_back(std::move(group));
        }
    }
    return groups;
}

ImmediateFrontierDispatchPolicy* ImmediateFrontierDispatchPolicy::clone() const {
    return new ImmediateFrontierDispatchPolicy(*this);
}

bool ImmediateFrontierDispatchPolicy::should_defer(
    const HybridEvolutionTree&, CheckpointId, const std::vector<CheckpointId>&,
    std::size_t) const {
    return false;
}

BoundedTimeLagFrontierDispatchPolicy::BoundedTimeLagFrontierDispatchPolicy(
    std::size_t max_deferrals, double max_time_lag)
    : _max_deferrals(max_deferrals), _max_time_lag(max_time_lag) {
    if (!std::isfinite(_max_time_lag) || _max_time_lag < 0.0)
        throw std::invalid_argument("Frontier dispatch requires a finite nonnegative time lag");
}

BoundedTimeLagFrontierDispatchPolicy* BoundedTimeLagFrontierDispatchPolicy::clone() const {
    return new BoundedTimeLagFrontierDispatchPolicy(*this);
}

bool BoundedTimeLagFrontierDispatchPolicy::should_defer(
    const HybridEvolutionTree& tree, CheckpointId candidate,
    const std::vector<CheckpointId>& others, std::size_t previous_deferrals) const {
    if (previous_deferrals >= _max_deferrals)
        return false;
    const auto candidate_time = tree.object(candidate.trajectory).entries.at(candidate.entry).work.enclosure.time_range();
    const double candidate_lower = candidate_time.lower_bound().raw().get_d();
    if (!std::isfinite(candidate_lower))
        return false;
    for (const auto other: others) {
        if (other.trajectory == candidate.trajectory)
            continue;
        const auto other_time = tree.object(other.trajectory).entries.at(other.entry).work.enclosure.time_range();
        const double other_lower = other_time.lower_bound().raw().get_d();
        if (std::isfinite(other_lower) && other_lower < candidate_lower &&
            candidate_lower - other_lower <= _max_time_lag)
            return true;
    }
    return false;
}

AcceptAllFrontierMergesPolicy* AcceptAllFrontierMergesPolicy::clone() const {
    return new AcceptAllFrontierMergesPolicy(*this);
}

bool AcceptAllFrontierMergesPolicy::accept(
    const HybridEvolutionTree&, const std::vector<CheckpointId>&,
    const EnclosureRecombinationResult&) const {
    return true;
}

RadiusTrendFrontierMergesPolicy::RadiusTrendFrontierMergesPolicy(
    double max_output_amplification, double max_continuation_growth)
    : _max_output_amplification(max_output_amplification),
      _max_continuation_growth(max_continuation_growth) {
    if (!std::isfinite(_max_output_amplification) || _max_output_amplification <= 0.0 ||
        !std::isfinite(_max_continuation_growth) || _max_continuation_growth <= 0.0)
        throw std::invalid_argument("Radius-trend merge limits must be finite and strictly positive");
}

RadiusTrendFrontierMergesPolicy* RadiusTrendFrontierMergesPolicy::clone() const {
    return new RadiusTrendFrontierMergesPolicy(*this);
}

bool RadiusTrendFrontierMergesPolicy::accept(
    const HybridEvolutionTree& tree, const std::vector<CheckpointId>& sources,
    const EnclosureRecombinationResult& result) const {
    if (sources.empty() || result.outputs.empty())
        return false;
    double maximum_input_radius = 0.0;
    for (const auto source: sources) {
        const auto& entries = tree.object(source.trajectory).entries;
        const auto& checkpoint = entries.at(source.entry);
        if (!checkpoint.origin || checkpoint.work.phase != HybridEvolutionPhase::EVOLVE_CONTINUOUSLY)
            return false;
        const auto& predecessor = entries.at(checkpoint.origin->entry);
        if (!predecessor.result ||
            predecessor.result->successors.at(checkpoint.origin->successor).kind !=
                HybridEvolutionSuccessorKind::CONTINUATION)
            return false;
        const auto previous_time = predecessor.work.enclosure.time_range();
        const auto current_time = checkpoint.work.enclosure.time_range();
        if (!(current_time.lower_bound().raw().get_d() >
              previous_time.lower_bound().raw().get_d()))
            return false;
        const double previous_radius = predecessor.work.enclosure.continuous_set().radius().raw().get_d();
        const double current_radius = checkpoint.work.enclosure.continuous_set().radius().raw().get_d();
        if (!std::isfinite(previous_radius) || !std::isfinite(current_radius) ||
            current_radius > _max_continuation_growth * previous_radius)
            return false;
        maximum_input_radius = std::max(maximum_input_radius, current_radius);
    }
    for (const auto& output: result.outputs) {
        const double radius = output.enclosure.radius().raw().get_d();
        if (!std::isfinite(radius) ||
            radius > _max_output_amplification * maximum_input_radius)
            return false;
    }
    return true;
}

std::vector<std::vector<HybridEvolutionTree::NodeId>>
    SiblingWindowRecombinationSchedulingPolicy::select(const HybridEvolutionTree& tree, NodeId parent) const {
    std::vector<std::vector<NodeId>> groups;
    std::vector<NodeId> group;
    double group_start = 0.0;
    double group_end = 0.0;
    const auto flush = [&] {
        if (group.size() >= 2u)
            groups.push_back(group);
        group.clear();
    };
    for (const auto child: tree.children(parent)) {
        const auto& metadata = tree.metadata(child);
        if (metadata.retired_for_recombination || metadata.complete || !metadata.origin ||
            tree.object(child).entries.size() != 1u ||
            tree.object(child).entries.front().result) {
            flush();
            continue;
        }
        const auto& enclosure = tree.object(child).entries.front().work.enclosure;
        const auto time = enclosure.time_range();
        const double start = time.lower_bound().raw().get_d();
        const double end = time.upper_bound().raw().get_d();
        if (!std::isfinite(start) || !std::isfinite(end) || end - start > _max_time_spread) {
            flush();
            continue;
        }
        if (!group.empty() &&
            (group.size() == _max_group_size ||
             !compatible_recombination_children(tree, parent, group.front(), child) ||
             std::max(group_end, end) - std::min(group_start, start) > _max_time_spread))
            flush();
        if (group.empty()) {
            group_start = start;
            group_end = end;
        } else {
            group_start = std::min(group_start, start);
            group_end = std::max(group_end, end);
        }
        group.push_back(child);
    }
    flush();
    return groups;
}

AllCompatibleSiblingRecombinationSchedulingPolicy*
    AllCompatibleSiblingRecombinationSchedulingPolicy::clone() const {
    return new AllCompatibleSiblingRecombinationSchedulingPolicy(*this);
}

std::vector<std::vector<HybridEvolutionTree::NodeId>>
    AllCompatibleSiblingRecombinationSchedulingPolicy::select(
        const HybridEvolutionTree& tree, NodeId parent) const {
    std::vector<std::vector<NodeId>> groups;
    for (const auto child: tree.children(parent)) {
        const auto& metadata = tree.metadata(child);
        if (metadata.retired_for_recombination || metadata.complete || !metadata.origin ||
            tree.object(child).entries.size() != 1u ||
            tree.object(child).entries.front().result)
            continue;

        bool grouped = false;
        for (auto& group: groups) {
            if (compatible_recombination_children(tree, parent, group.front(), child)) {
                group.push_back(child);
                grouped = true;
                break;
            }
        }
        if (!grouped)
            groups.push_back({child});
    }
    groups.erase(std::remove_if(groups.begin(), groups.end(),
                                [](const auto& group) { return group.size() < 2u; }),
                 groups.end());
    return groups;
}

/** @brief Run state, kept out of the installed backend class's representation. */
class TreeHybridEvolutionBackend::Impl {
  public:
    HybridEvolutionTree history;
    HybridEvolutionStatistics statistics;
    std::unique_ptr<const HybridRecombinationSchedulingPolicy> policy;
    std::unique_ptr<const HybridFrontierRecombinationPolicy> frontier_policy;
    std::unique_ptr<const HybridFrontierDispatchPolicy> frontier_dispatch;
    std::unique_ptr<const HybridFrontierMergeAcceptancePolicy> frontier_acceptance;
    std::unique_ptr<const HybridCheckpointAdmissionPolicy> admission;
    std::unique_ptr<const HybridRecombinationTriggerPolicy> trigger;
    bool use_frontier = false;
    std::size_t global_steps = 0u;
    std::unique_ptr<const EnclosureRecombinerInterface> recombiner;
    EnclosureRecombinationContext context;
    Semantics semantics = Semantics::LOWER;
    Bool started = false;
    Bool complete = false;
    Bool released = false;

    /** @brief Create a trajectory seed without an implicit Orbit contribution. */
    explicit Impl(const HybridEnclosure& initial,
                  const HybridRecombinationSchedulingPolicy& grouping,
                  const HybridRecombinationTriggerPolicy& timing,
                  const EnclosureRecombinerInterface& method,
                  const EnclosureRecombinationContext& limits)
        : policy(grouping.clone()), trigger(timing.clone()),
          recombiner(method.clone()), context(limits) {
        HybridEvolutionTrajectory trajectory;
        trajectory.entries.push_back({
            {initial, HybridEvolutionPhase::PROCESS_STARTING_EVENTS},
            std::nullopt,
            std::nullopt
        });
        history.create_root(std::move(trajectory), HybridEvolutionTrajectoryMetadata{});
    }

    /** @brief Initialize a globally scheduled frontier with independent policy state. */
    explicit Impl(const HybridEnclosure& initial,
                  const HybridFrontierRecombinationPolicy& grouping,
                  const HybridRecombinationTriggerPolicy& timing,
                  const EnclosureRecombinerInterface& method,
                  const EnclosureRecombinationContext& limits,
                  const HybridFrontierDispatchPolicy& dispatch,
                  const HybridFrontierMergeAcceptancePolicy& acceptance)
        : frontier_policy(grouping.clone()), frontier_dispatch(dispatch.clone()),
          frontier_acceptance(acceptance.clone()),
          admission(std::make_unique<ImmediateCheckpointAdmissionPolicy>()), trigger(timing.clone()),
          use_frontier(true), recombiner(method.clone()), context(limits) {
        HybridEvolutionTrajectory trajectory;
        trajectory.entries.push_back({
            {initial, HybridEvolutionPhase::PROCESS_STARTING_EVENTS},
            std::nullopt,
            std::nullopt
        });
        HybridEvolutionTrajectoryMetadata metadata;
        metadata.pending_entries = 1u;
        history.create_root(std::move(trajectory), std::move(metadata));
    }

    /** @brief Classify successors using explicit provenance and whole-source evidence. */
    static Bool derives_trajectory(const HybridEvolutionSuccessor& successor, std::size_t successor_count, Bool whole_trajectory) {
        switch (successor.kind) {
            case HybridEvolutionSuccessorKind::TRANSITION:
                return !(whole_trajectory && successor_count == 1 &&
                         successor.resolution == HybridEvolutionTransitionResolution::FULL_SOURCE &&
                         (successor.event_kind == EventKind::URGENT || successor.event_kind == EventKind::IMPACT));
            case HybridEvolutionSuccessorKind::CONTINUATION:
            case HybridEvolutionSuccessorKind::SUBDIVISION:
            case HybridEvolutionSuccessorKind::RECONDITIONING:
                return false;
            default:
                throw std::invalid_argument("Unknown hybrid evolution successor kind");
        }
    }

    /** @brief Count valid contributions without maintaining a second authoritative Orbit. */
    Void count(const HybridEvolutionStepResult& result) {
        for (const auto& contribution: result.contributions) {
            switch (contribution.kind) {
                case HybridEvolutionSetKind::REACH:
                    ++statistics.reach;
                    break;
                case HybridEvolutionSetKind::INTERMEDIATE:
                    ++statistics.intermediate;
                    break;
                case HybridEvolutionSetKind::FINAL:
                    ++statistics.final;
                    break;
                default:
                    throw std::invalid_argument("Unknown hybrid evolution contribution kind");
            }
        }
    }

    /** @brief Replace accepted sibling groups at safe work-publication boundaries.
     *  No child has started evolving before its parent completes; even during
     *  per-step evaluation, only unevolved child checkpoints are selected.
     */
    Void schedule_recombination(HybridEvolutionTree::NodeId parent,
                                std::vector<HybridEvolutionTree::NodeId>& waiting_children) {
        LOGGING_SCOPE_CREATE;
        if (!policy->enables_recombination())
            return;
        const auto groups = policy->select(history, parent);
        // Validate all proposed groups before modifying the genealogy.
        std::vector<HybridEvolutionTree::NodeId> selected;
        for (const auto& group: groups) {
            if (group.size() < 2u)
                throw std::logic_error("Recombination policy selected a group smaller than two");
            for (const auto source: group) {
                if (!history.contains(source) || history.parent(source) != parent ||
                    history.metadata(source).retired_for_recombination ||
                    history.metadata(source).complete ||
                    history.object(source).entries.size() != 1u ||
                    history.object(source).entries.front().result ||
                    std::find(selected.begin(), selected.end(), source) != selected.end())
                    throw std::logic_error("Recombination policy selected invalid or overlapping children");
                if (source != group.front() &&
                    !compatible_recombination_children(history, parent, group.front(), source))
                    throw std::logic_error("Recombination policy selected incompatible siblings");
                selected.push_back(source);
            }
        }
        for (const auto& group: groups) {
            List<Enclosure> sources;
            for (const auto source: group)
                sources.append(history.object(source).entries.front().work.enclosure.continuous_set());
            const auto result = recombiner->recombine(sources, context);
            // An unsuccessful numerical proposal leaves all children runnable.
            if (result.outputs.empty() || result.outputs.size() >= group.size() ||
                !result.discarded_sources.empty() || !result.radius_limit_satisfied) {
                LOGGING_PRINTLN_AT(4, "RECOMBINE rejected parent-first"
                                          << " sources=" << group.size()
                                          << " outputs=" << result.outputs.size()
                                          << " discarded=" << result.discarded_sources.size()
                                          << " radius_ok=" << result.radius_limit_satisfied);
                continue;
            }
            std::vector<bool> represented(group.size(), false);
            bool complete_coverage = true;
            for (const auto& output: result.outputs) {
                if (output.source_indices.empty())
                    complete_coverage = false;
                for (const auto index: output.source_indices) {
                    if (index >= group.size()) {
                        complete_coverage = false;
                    } else {
                        represented[index] = true;
                    }
                }
            }
            if (!std::all_of(represented.begin(), represented.end(), [](bool covered) { return covered; }) ||
                !complete_coverage) {
                LOGGING_PRINTLN_AT(4, "RECOMBINE rejected parent-first: incomplete provenance"
                                          << " sources=" << group.size() << " outputs=" << result.outputs.size());
                continue;
            }
            LOGGING_PRINTLN("RECOMBINE parent-first"
                            << " location=" << history.object(group.front()).entries.front().work.enclosure.location()
                            << " t=" << history.object(group.front()).entries.front().work.enclosure.time_range()
                            << " sources=" << group.size() << " outputs=" << result.outputs.size()
                            << " max_source_radius=" << largest_source_radius(sources)
                            << " max_output_radius=" << largest_output_radius(result));
            // Commit all replacements before retiring the original sources.
            // Invalid or ineffective results leave every original child runnable.
            for (const auto& output: result.outputs) {
                const auto representative = group[output.source_indices.front()];
                auto replacement = history.object(representative).entries.front().work;
                replacement.enclosure.replace_continuous_set(output.enclosure);
                HybridEvolutionTrajectory trajectory;
                trajectory.entries.push_back({replacement, std::nullopt, std::nullopt});
                HybridEvolutionTrajectoryMetadata metadata;
                for (const auto index: output.source_indices)
                    metadata.recombined_sources.push_back(group[index]);
                const auto created = history.add_child(parent, std::move(trajectory), std::move(metadata));
                // Accepted replacements remain blocked like their sources.
                // Only completion of this scheduling group can admit them.
                waiting_children.push_back(created);
            }
            for (const auto source: group)
                history.metadata(source).retired_for_recombination = true;
        }
    }

    /** @brief Find the closest existing common ancestor for a multi-branch replacement. */
    HybridEvolutionTree::NodeId common_ancestor(
        const std::vector<HybridEvolutionCheckpointId>& sources) const {
        const auto first = sources.front().trajectory;
        std::vector<HybridEvolutionTree::NodeId> candidates{first};
        const auto ancestors = history.ancestors(first);
        candidates.insert(candidates.end(), ancestors.begin(), ancestors.end());
        for (const auto candidate: candidates) {
            bool shared = true;
            for (const auto source: sources) {
                if (source.trajectory == candidate)
                    continue;
                const auto parents = history.ancestors(source.trajectory);
                if (std::find(parents.begin(), parents.end(), candidate) == parents.end()) {
                    shared = false;
                    break;
                }
            }
            if (shared)
                return candidate;
        }
        throw std::logic_error("Frontier checkpoints have no common genealogy root");
    }

    /** @brief Recombine compatible pending checkpoints from any part of the tree.
     *  A processed entry is never replaced; all selected source IDs remain in
     *  append-only history with their original result-less checkpoint records.
     */
    Void schedule_frontier_recombination(std::deque<HybridEvolutionCheckpointId>& ready) {
        LOGGING_SCOPE_CREATE;
        if (!frontier_policy->enables_recombination())
            return;
        std::vector<HybridEvolutionCheckpointId> pending;
        for (const auto checkpoint: ready) {
            const auto& trajectory = history.object(checkpoint.trajectory);
            if (trajectory.entries.at(checkpoint.entry).result)
                continue;
            if (std::find(trajectory.retired_entries.begin(), trajectory.retired_entries.end(),
                          checkpoint.entry) == trajectory.retired_entries.end())
                pending.push_back(checkpoint);
        }
        if (pending.size() < 2u)
            return;
        const auto groups = frontier_policy->select(history, pending);
        std::vector<HybridEvolutionCheckpointId> selected;
        for (const auto& group: groups) {
            if (group.size() < 2u)
                throw std::logic_error("Frontier policy selected fewer than two checkpoints");
            for (const auto checkpoint: group) {
                if (!history.contains(checkpoint.trajectory) ||
                    std::find(pending.begin(), pending.end(), checkpoint) == pending.end() ||
                    std::find(selected.begin(), selected.end(), checkpoint) != selected.end())
                    throw std::logic_error("Frontier policy selected unknown or overlapping checkpoint");
                if (checkpoint != group.front() &&
                    !compatible_frontier_checkpoints(history, group.front(), checkpoint,
                                                     frontier_policy->requires_time_overlap()))
                    throw std::logic_error("Frontier policy selected incompatible checkpoints");
                selected.push_back(checkpoint);
            }
        }
        for (const auto& group: groups) {
            // A merged enclosure can be subdivided again immediately by the
            // numerical evolver. Do not let those pieces recombine at the same
            // physical time: every source descending from a previous merge
            // must have made strictly positive progress since that merge.
            // Held replacements cannot execute or resplit before admission;
            // they may absorb later arrivals without advancing physical time.
            bool progressed = true;
            double merge_time = -std::numeric_limits<double>::infinity();
            for (const auto checkpoint: group) {
                const auto& source_metadata = history.metadata(checkpoint.trajectory);
                const auto time = history.object(checkpoint.trajectory).entries.at(checkpoint.entry).work.enclosure.time_range();
                const double lower = time.lower_bound().raw().get_d();
                const double upper = time.upper_bound().raw().get_d();
                if (!std::isfinite(lower) || !std::isfinite(upper) ||
                    (admission->recombination_queue() == HybridRecombinationQueue::READY &&
                     source_metadata.recombination_time_barrier &&
                     !(lower > *source_metadata.recombination_time_barrier)))
                    progressed = false;
                merge_time = std::max(merge_time, upper);
            }
            if (!progressed)
                continue;
            List<Enclosure> sources;
            for (const auto checkpoint: group)
                sources.append(history.object(checkpoint.trajectory).entries.at(checkpoint.entry).work.enclosure.continuous_set());
            const auto result = recombiner->recombine(sources, context);
            if (result.outputs.empty() || result.outputs.size() >= group.size() ||
                !result.discarded_sources.empty() || !result.radius_limit_satisfied) {
                LOGGING_PRINTLN_AT(4, "RECOMBINE rejected frontier"
                                          << " sources=" << group.size()
                                          << " outputs=" << result.outputs.size()
                                          << " discarded=" << result.discarded_sources.size()
                                          << " radius_ok=" << result.radius_limit_satisfied);
                continue;
            }
            std::vector<bool> represented(group.size(), false);
            bool full_coverage = true;
            for (const auto& output: result.outputs) {
                if (output.source_indices.empty())
                    full_coverage = false;
                for (const auto index: output.source_indices) {
                    if (index >= group.size()) {
                        full_coverage = false;
                    } else {
                        represented[index] = true;
                    }
                }
            }
            if (!full_coverage ||
                !std::all_of(represented.begin(), represented.end(), [](bool covered) { return covered; })) {
                LOGGING_PRINTLN_AT(4, "RECOMBINE rejected frontier: incomplete provenance"
                                          << " sources=" << group.size() << " outputs=" << result.outputs.size());
                continue;
            }
            if (!frontier_acceptance->accept(history, group, result)) {
                LOGGING_PRINTLN_AT(4, "RECOMBINE rejected frontier: quality gate"
                                          << " sources=" << group.size() << " outputs=" << result.outputs.size());
                continue;
            }
            const auto& first = history.object(group.front().trajectory).entries.at(group.front().entry).work.enclosure;
            LOGGING_PRINTLN("RECOMBINE frontier"
                            << " location=" << first.location()
                            << " t=" << first.time_range()
                            << " sources=" << group.size() << " outputs=" << result.outputs.size()
                            << " max_source_radius=" << largest_source_radius(sources)
                            << " max_output_radius=" << largest_output_radius(result));
            const auto parent = common_ancestor(group);
            for (const auto& output: result.outputs) {
                const auto representative = group[output.source_indices.front()];
                auto work = history.object(representative.trajectory)
                                .entries.at(representative.entry)
                                .work;
                work.enclosure.replace_continuous_set(output.enclosure);
                HybridEvolutionTrajectory trajectory;
                trajectory.entries.push_back({work, std::nullopt, std::nullopt});
                HybridEvolutionTrajectoryMetadata metadata;
                metadata.pending_entries = 1u;
                metadata.recombination_time_barrier = merge_time;
                for (const auto index: output.source_indices)
                    metadata.recombined_checkpoints.push_back(group[index]);
                const auto created = history.add_child(parent, std::move(trajectory), std::move(metadata));
                ready.push_back({created, 0u});
            }
            for (const auto checkpoint: group) {
                history.object(checkpoint.trajectory).retired_entries.push_back(checkpoint.entry);
                auto& metadata = history.metadata(checkpoint.trajectory);
                if (metadata.pending_entries == 0u)
                    throw std::logic_error("Invalid checkpoint pending count");
                --metadata.pending_entries;
                if (metadata.pending_entries == 0u)
                    metadata.complete = true;
            }
        }
    }

    /** @brief Execute admissible checkpoints with serialized result publication.
     *  @details The same scheduler handles zero or several worker threads.
     *  Every ready checkpoint has one DynamicWorkload token; the token
     *  takes a checkpoint from the protected queue before beginning numerical
     *  work. In parent-first mode, discrete children are registered in an
     *  explicit waiting group, separate from the genealogy, until the last
     *  in-flight or queued continuation of their parent has published its
     *  result. In frontier mode, only ready (not in-flight)
     *  checkpoints can be grouped or retired for recombination.
     */
    Void run_scheduled(const StepFunction& advance) {
        struct Token {};
        struct ActiveGroup {
            HybridEvolutionTree::NodeId node;
            // Counts admissible + in-flight checkpoints, not waiting children.
            std::size_t outstanding = 1u;
            Bool whole_trajectory = true;
            std::size_t completed_steps = 0u;
        };
        struct WaitingGroup {
            HybridEvolutionTree::NodeId node;
            // Ordered discrete descendants, including accepted replacements.
            std::vector<HybridEvolutionTree::NodeId> children;
        };
        using Workload = DynamicWorkload<Token>;
        const auto root = history.root().value();
        std::mutex schedule_mutex;
        std::deque<HybridEvolutionCheckpointId> ready{
            {root, 0u}
        };
        std::deque<HybridEvolutionCheckpointId> held_frontier;
        std::size_t in_flight = 0u;
        std::vector<std::pair<HybridEvolutionCheckpointId, std::size_t>> deferrals;
        std::vector<ActiveGroup> active;
        std::vector<WaitingGroup> waiting;
        if (!use_frontier) {
            active.push_back({root, 1u, true, 0u});
            waiting.push_back({root, {}});
        }

        Workload workload(
            [](const Token&, SharedPointer<ProgressIndicator>) {},
            [&](Workload::Access& producer, const Token&) {
                HybridEvolutionCheckpointId checkpoint;
                std::optional<HybridEvolutionWorkItem> work;
                HybridEvolutionStatistics snapshot;
                {
                    std::lock_guard<std::mutex> lock(schedule_mutex);
                    // Tokens are allowed to outlive ready checkpoints retired
                    // by a concurrent recombination: such tokens do no work.
                    while (!ready.empty()) {
                        const auto candidate = ready.front();
                        ready.pop_front();
                        const auto& trajectory = history.object(candidate.trajectory);
                        if (std::find(trajectory.retired_entries.begin(),
                                      trajectory.retired_entries.end(), candidate.entry) !=
                            trajectory.retired_entries.end())
                            continue;
                        if (trajectory.entries.at(candidate.entry).result)
                            throw std::logic_error("Checkpoint cannot be executed twice");

                        if (use_frontier) {
                            const auto deferred = std::find_if(deferrals.begin(), deferrals.end(),
                                                               [&](const auto& entry) { return entry.first == candidate; });
                            const auto previous = deferred == deferrals.end() ? 0u : deferred->second;
                            std::vector<HybridEvolutionCheckpointId> alternatives;
                            for (const auto other: ready) {
                                const auto& other_trajectory = history.object(other.trajectory);
                                if (!other_trajectory.entries.at(other.entry).result &&
                                    std::find(other_trajectory.retired_entries.begin(),
                                              other_trajectory.retired_entries.end(), other.entry) ==
                                        other_trajectory.retired_entries.end())
                                    alternatives.push_back(other);
                            }
                            if (!alternatives.empty() &&
                                frontier_dispatch->should_defer(history, candidate, alternatives, previous)) {
                                if (deferred == deferrals.end())
                                    deferrals.push_back({candidate, 1u});
                                else
                                    ++deferred->second;
                                ready.push_back(candidate);
                                continue;
                            }
                            if (deferred != deferrals.end())
                                deferrals.erase(deferred);
                        }
                        checkpoint = candidate;
                        work.emplace(trajectory.entries.at(candidate.entry).work);
                        snapshot = statistics;
                        ++in_flight;
                        break;
                    }
                }
                if (!work)
                    return;

                // Numerical evolution never holds the scheduler lock.
                // DynamicWorkload runs this call directly at concurrency zero
                // or distributes it among available workers otherwise.
                HybridEvolutionStepResult result = advance(*work, snapshot);

                std::lock_guard<std::mutex> lock(schedule_mutex);
                const auto enqueue = [&](HybridEvolutionCheckpointId next) {
                    ready.push_back(next);
                    producer.append(Token{});
                };
                const auto is_retired = [&](HybridEvolutionCheckpointId id) {
                    const auto& trajectory = history.object(id.trajectory);
                    return std::find(trajectory.retired_entries.begin(), trajectory.retired_entries.end(),
                                     id.entry) != trajectory.retired_entries.end();
                };
                const auto evaluate_frontier = [&] {
                    const bool ready_queue = admission->recombination_queue() == HybridRecombinationQueue::READY;
                    auto& candidates = ready_queue ? ready : held_frontier;
                    const auto previous_size = candidates.size();
                    schedule_frontier_recombination(candidates);
                    // Recombined replacements append new ready checkpoints;
                    // existing retired tokens may still be in the workload.
                    if (ready_queue) {
                        for (std::size_t i = previous_size; i != candidates.size(); ++i)
                            producer.append(Token{});
                        ready.erase(std::remove_if(ready.begin(), ready.end(), is_retired), ready.end());
                    }
                };

                count(result);
                --in_flight;
                if (use_frontier) {
                    const auto publish_successor = [&](HybridEvolutionCheckpointId next,
                                                       const HybridEvolutionSuccessor& successor) {
                        if (admission->admit_successor(*work, successor))
                            enqueue(next);
                        else
                            held_frontier.push_back(next);
                    };
                    for (const auto& successor: result.successors) {
                        if (successor.kind == HybridEvolutionSuccessorKind::SUBDIVISION)
                            history.metadata(checkpoint.trajectory).whole_trajectory = false;
                    }
                    for (std::size_t ordinal = 0; ordinal != result.successors.size(); ++ordinal) {
                        const auto& successor = result.successors[ordinal];
                        const HybridEvolutionEntryOrigin origin{checkpoint.entry, ordinal};
                        const bool split = successor.kind == HybridEvolutionSuccessorKind::SUBDIVISION;
                        if (split || derives_trajectory(
                                         successor, result.successors.size(),
                                         history.metadata(checkpoint.trajectory).whole_trajectory)) {
                            HybridEvolutionTrajectory child;
                            child.entries.push_back({successor.work, std::nullopt, std::nullopt});
                            HybridEvolutionTrajectoryMetadata metadata;
                            metadata.origin = origin;
                            metadata.pending_entries = 1u;
                            metadata.recombination_time_barrier =
                                history.metadata(checkpoint.trajectory).recombination_time_barrier;
                            if (split)
                                metadata.whole_trajectory = false;
                            const auto created = history.add_child(checkpoint.trajectory,
                                                                   std::move(child), std::move(metadata));
                            publish_successor({created, 0u}, successor);
                            if (successor.kind == HybridEvolutionSuccessorKind::TRANSITION &&
                                successor.event_kind != EventKind::PERMISSIVE)
                                history.metadata(checkpoint.trajectory).whole_trajectory = false;
                        } else {
                            auto& entries = history.object(checkpoint.trajectory).entries;
                            const auto index = entries.size();
                            entries.push_back({successor.work, origin, std::nullopt});
                            ++history.metadata(checkpoint.trajectory).pending_entries;
                            publish_successor({checkpoint.trajectory, index}, successor);
                        }
                    }
                    history.object(checkpoint.trajectory).entries.at(checkpoint.entry).result.emplace(std::move(result));
                    auto& metadata = history.metadata(checkpoint.trajectory);
                    if (metadata.pending_entries == 0u)
                        throw std::logic_error("Invalid checkpoint completion count");
                    --metadata.pending_entries;
                    const bool finished = metadata.pending_entries == 0u;
                    if (finished)
                        metadata.complete = true;
                    ++global_steps;
                    if (frontier_policy->enables_recombination() &&
                        trigger->should_evaluate(history, checkpoint.trajectory,
                                                 HybridRecombinationEvaluationPoint::AFTER_STEP, global_steps))
                        evaluate_frontier();
                    if (finished && frontier_policy->enables_recombination() &&
                        trigger->should_evaluate(history, checkpoint.trajectory,
                                                 HybridRecombinationEvaluationPoint::PARENT_COMPLETE, global_steps))
                        evaluate_frontier();
                    // Retired tokens may still exist, but are not numerical
                    // producers. Close the batch only after every claimed work
                    // item has published all of its successors.
                    if (ready.empty() && in_flight == 0u) {
                        if (frontier_policy->enables_recombination() &&
                            trigger->should_evaluate(history, checkpoint.trajectory,
                                                     HybridRecombinationEvaluationPoint::SYNCHRONIZATION_COMPLETE, global_steps))
                            evaluate_frontier();
                        held_frontier.erase(std::remove_if(held_frontier.begin(), held_frontier.end(), is_retired),
                                            held_frontier.end());
                        const std::vector<HybridEvolutionCheckpointId> pending(held_frontier.begin(), held_frontier.end());
                        const auto admitted = admission->select_next(history, pending);
                        std::vector<HybridEvolutionCheckpointId> unique;
                        for (const auto next: admitted) {
                            if (std::find(pending.begin(), pending.end(), next) == pending.end() ||
                                std::find(unique.begin(), unique.end(), next) != unique.end())
                                throw std::logic_error("Admission policy selected unknown or duplicate waiting work");
                            unique.push_back(next);
                        }
                        if (!pending.empty() && admitted.empty())
                            throw std::logic_error("Admission policy left waiting work without a runnable batch");
                        for (const auto next: admitted) {
                            held_frontier.erase(std::find(held_frontier.begin(), held_frontier.end(), next));
                            enqueue(next);
                        }
                    }
                } else {
                    auto stage = std::find_if(active.begin(), active.end(),
                                              [&](const auto& item) { return item.node == checkpoint.trajectory; });
                    auto held = std::find_if(waiting.begin(), waiting.end(),
                                             [&](const auto& item) { return item.node == checkpoint.trajectory; });
                    if (stage == active.end() || held == waiting.end() || stage->outstanding == 0u)
                        throw std::logic_error("Checkpoint has no active scheduling group");
                    for (const auto& successor: result.successors)
                        if (successor.kind == HybridEvolutionSuccessorKind::SUBDIVISION)
                            stage->whole_trajectory = false;
                    for (std::size_t ordinal = 0; ordinal != result.successors.size(); ++ordinal) {
                        const auto& successor = result.successors[ordinal];
                        const HybridEvolutionEntryOrigin origin{checkpoint.entry, ordinal};
                        if (derives_trajectory(successor, result.successors.size(), stage->whole_trajectory)) {
                            HybridEvolutionTrajectory child;
                            child.entries.push_back({successor.work, std::nullopt, std::nullopt});
                            HybridEvolutionTrajectoryMetadata metadata;
                            metadata.origin = origin;
                            const auto created = history.add_child(
                                checkpoint.trajectory, std::move(child), std::move(metadata));
                            // The waiting queue, rather than a genealogy scan,
                            // owns the scheduler's list of blocked successors.
                            held->children.push_back(created);
                            if (successor.event_kind != EventKind::PERMISSIVE)
                                stage->whole_trajectory = false;
                        } else {
                            auto& entries = history.object(checkpoint.trajectory).entries;
                            const auto index = entries.size();
                            entries.push_back({successor.work, origin, std::nullopt});
                            ++stage->outstanding;
                            enqueue({checkpoint.trajectory, index});
                        }
                    }
                    history.object(checkpoint.trajectory).entries.at(checkpoint.entry).result.emplace(std::move(result));
                    ++stage->completed_steps;
                    if (policy->enables_recombination() &&
                        trigger->should_evaluate(history, checkpoint.trajectory,
                                                 HybridRecombinationEvaluationPoint::AFTER_STEP, stage->completed_steps))
                        schedule_recombination(checkpoint.trajectory, held->children);
                    if (--stage->outstanding == 0u) {
                        // All producers of this parent's discrete children
                        // have published, including those spawned by splits.
                        history.metadata(checkpoint.trajectory).complete = true;
                        if (policy->enables_recombination() &&
                            trigger->should_evaluate(history, checkpoint.trajectory,
                                                     HybridRecombinationEvaluationPoint::PARENT_COMPLETE,
                                                     history.object(checkpoint.trajectory).entries.size()))
                            schedule_recombination(checkpoint.trajectory, held->children);
                        // Last publisher closes the active group and drains its
                        // waiting queue atomically. Retired sources stay only
                        // in the genealogy, never in the ready queue.
                        const auto ready_children = std::move(held->children);
                        waiting.erase(held);
                        active.erase(stage);
                        for (const auto child: ready_children) {
                            if (history.metadata(child).retired_for_recombination)
                                continue;
                            active.push_back({child, 1u, true, 0u});
                            waiting.push_back({child, {}});
                            enqueue({child, 0u});
                        }
                    }
                }
            });
        workload.append(Token{}).process();
    }

    /** @brief Append one trajectory's explicit contributions in creation order. */
    Void project_trajectory(HybridEvolutionTree::NodeId node, Orbit<HybridEnclosure>& orbit) const {
        const auto& trajectory = history.object(node);
        for (std::size_t index = 0; index != trajectory.entries.size(); ++index) {
            const auto& entry = trajectory.entries[index];
            if (!entry.result) {
                if (std::find(trajectory.retired_entries.begin(), trajectory.retired_entries.end(),
                              index) != trajectory.retired_entries.end())
                    continue;
                throw std::logic_error("Cannot project an unprocessed checkpoint into Orbit");
            }
            for (const auto& contribution: entry.result.value().contributions) {
                switch (contribution.kind) {
                    case HybridEvolutionSetKind::REACH:
                        orbit.adjoin_reach(contribution.enclosure);
                        break;
                    case HybridEvolutionSetKind::INTERMEDIATE:
                        orbit.adjoin_intermediate(contribution.enclosure);
                        break;
                    case HybridEvolutionSetKind::FINAL:
                        orbit.adjoin_final(contribution.enclosure);
                        break;
                    default:
                        throw std::invalid_argument("Unknown hybrid evolution contribution kind");
                }
            }
        }
    }
};

TreeHybridEvolutionBackend::TreeHybridEvolutionBackend(
    const HybridEnclosure& initial, const HybridCheckpointAdmissionPolicy& admission,
    const HybridFrontierRecombinationPolicy& grouping,
    const HybridRecombinationTriggerPolicy& timing, const EnclosureRecombinerInterface& method,
    const EnclosureRecombinationContext& limits, const HybridFrontierDispatchPolicy& dispatch,
    const HybridFrontierMergeAcceptancePolicy& acceptance)
    : TreeHybridEvolutionBackend(initial, grouping, timing, method, limits, dispatch, acceptance) {
    _impl->admission.reset(admission.clone());
}

TreeHybridEvolutionBackend::TreeHybridEvolutionBackend(
    const HybridEnclosure& initial, const HybridRecombinationSchedulingPolicy& grouping,
    const HybridRecombinationTriggerPolicy& timing, const EnclosureRecombinerInterface& method,
    const EnclosureRecombinationContext& limits)
    : _impl(std::make_unique<Impl>(initial, grouping, timing, method, limits)) {}

TreeHybridEvolutionBackend::TreeHybridEvolutionBackend(
    const HybridEnclosure& initial, const HybridFrontierRecombinationPolicy& grouping,
    const HybridRecombinationTriggerPolicy& timing, const EnclosureRecombinerInterface& method,
    const EnclosureRecombinationContext& limits,
    const HybridFrontierDispatchPolicy& dispatch,
    const HybridFrontierMergeAcceptancePolicy& acceptance)
    : _impl(std::make_unique<Impl>(initial, grouping, timing, method, limits, dispatch, acceptance)) {}

TreeHybridEvolutionBackend::~TreeHybridEvolutionBackend() = default;

Void TreeHybridEvolutionBackend::set_semantics(Semantics semantics) {
    if (_impl->started)
        throw std::logic_error("Cannot change recombination semantics after the run starts");
    _impl->semantics = semantics;
}

Void TreeHybridEvolutionBackend::run(const StepFunction& advance) {
    if (_impl->started) {
        throw std::logic_error("A hybrid evolution backend can only run once");
    }
    if ((_impl->use_frontier ? _impl->frontier_policy->enables_recombination()
                             : _impl->policy->enables_recombination()) &&
        _impl->semantics != Semantics::UPPER)
        throw std::invalid_argument("Recombination is only sound for upper semantics");
    _impl->started = true;
    // The scheduler is identical for sequential and concurrent execution:
    // DynamicWorkload controls how many ready checkpoints run at once.
    _impl->run_scheduled(advance);
    _impl->complete = true;
}

Orbit<HybridEnclosure> TreeHybridEvolutionBackend::release_orbit() {
    if (!_impl->complete || _impl->released) {
        throw std::logic_error("No completed hybrid Orbit is available");
    }
    const auto root = _impl->history.root().value();
    Orbit<HybridEnclosure> result(_impl->history.object(root).entries.front().work.enclosure);
    _impl->project_trajectory(root, result);
    for (const auto node: _impl->history.descendants(root)) {
        if (!_impl->history.metadata(node).retired_for_recombination)
            _impl->project_trajectory(node, result);
    }
    _impl->released = true;
    return result;
}

const HybridEvolutionTree& TreeHybridEvolutionBackend::tree() const noexcept {
    return _impl->history;
}

TreeHybridEvolutionBackendFactory::TreeHybridEvolutionBackendFactory()
    : TreeHybridEvolutionBackendFactory(NoRecombinationSchedulingPolicy{}) {}

TreeHybridEvolutionBackendFactory::TreeHybridEvolutionBackendFactory(
    const HybridCheckpointAdmissionPolicy& admission, const HybridFrontierRecombinationPolicy& grouping,
    const HybridRecombinationTriggerPolicy& timing, const EnclosureRecombinerInterface& method,
    const EnclosureRecombinationContext& limits, const HybridFrontierDispatchPolicy& dispatch,
    const HybridFrontierMergeAcceptancePolicy& acceptance)
    : TreeHybridEvolutionBackendFactory(grouping, timing, method, limits, dispatch, acceptance) {
    _admission.reset(admission.clone());
}

TreeHybridEvolutionBackendFactory::TreeHybridEvolutionBackendFactory(
    const HybridRecombinationSchedulingPolicy& grouping,
    const HybridRecombinationTriggerPolicy& timing,
    const EnclosureRecombinerInterface& method,
    const EnclosureRecombinationContext& limits)
    : _grouping(grouping.clone()), _trigger(timing.clone()),
      _recombiner(method.clone()), _context(limits) {}

TreeHybridEvolutionBackendFactory::TreeHybridEvolutionBackendFactory(
    const HybridFrontierRecombinationPolicy& grouping,
    const HybridRecombinationTriggerPolicy& timing,
    const EnclosureRecombinerInterface& method,
    const EnclosureRecombinationContext& limits,
    const HybridFrontierDispatchPolicy& dispatch,
    const HybridFrontierMergeAcceptancePolicy& acceptance)
    : _frontier_grouping(grouping.clone()), _frontier_dispatch(dispatch.clone()),
      _frontier_acceptance(acceptance.clone()), _trigger(timing.clone()),
      _recombiner(method.clone()), _context(limits) {}

TreeHybridEvolutionBackendFactory* TreeHybridEvolutionBackendFactory::clone() const {
    return new TreeHybridEvolutionBackendFactory(*this);
}

std::unique_ptr<HybridEvolutionBackendInterface> TreeHybridEvolutionBackendFactory::create(const HybridEnclosure& initial, const Real&) const {
    if (_admission)
        return std::make_unique<TreeHybridEvolutionBackend>(
            initial, *_admission, *_frontier_grouping, *_trigger, *_recombiner, _context,
            *_frontier_dispatch, *_frontier_acceptance);
    if (_frontier_grouping)
        return std::make_unique<TreeHybridEvolutionBackend>(
            initial, *_frontier_grouping, *_trigger, *_recombiner, _context,
            *_frontier_dispatch, *_frontier_acceptance);
    return std::make_unique<TreeHybridEvolutionBackend>(
        initial, *_grouping, *_trigger, *_recombiner, _context);
}

} // namespace Ariadne
