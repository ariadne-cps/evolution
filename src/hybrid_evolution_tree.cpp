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
#include <stdexcept>
#include <utility>

#include "evolution/enclosure_recombiner.hpp"
#include "function/formula.hpp"
#include "function/formula.tpl.hpp"

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
    return point == HybridRecombinationEvaluationPoint::PARENT_COMPLETE &&
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
    const auto& first_successor = parent_entries.at(first_metadata.origin->entry).result.value()
                                      .successors.at(first_metadata.origin->successor);
    const auto& next_successor = parent_entries.at(next_metadata.origin->entry).result.value()
                                     .successors.at(next_metadata.origin->successor);
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
                                    HybridEvolutionCheckpointId next) {
    const auto& a = tree.object(first.trajectory).entries.at(first.entry).work;
    const auto& b = tree.object(next.trajectory).entries.at(next.entry).work;
    const auto ta = a.enclosure.time_range();
    const auto tb = b.enclosure.time_range();
    if (ta.lower_bound().raw().get_d() > tb.upper_bound().raw().get_d() ||
        tb.lower_bound().raw().get_d() > ta.upper_bound().raw().get_d())
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
            const auto time = tree.object(pending[next].trajectory).entries.at(pending[next].entry)
                                  .work.enclosure.time_range();
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
    const auto candidate_time = tree.object(candidate.trajectory).entries.at(candidate.entry)
                                    .work.enclosure.time_range();
    const double candidate_lower = candidate_time.lower_bound().raw().get_d();
    if (!std::isfinite(candidate_lower))
        return false;
    for (const auto other: others) {
        if (other.trajectory == candidate.trajectory)
            continue;
        const auto other_time = tree.object(other.trajectory).entries.at(other.entry)
                                    .work.enclosure.time_range();
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
          frontier_acceptance(acceptance.clone()), trigger(timing.clone()),
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

    /** @brief Complete all pieces of one trajectory; leave its children unprocessed. */
    Void process_trajectory(HybridEvolutionTree::NodeId node, const StepFunction& advance) {
        std::deque<std::size_t> pending{0};
        Bool whole_trajectory = true;
        std::size_t completed_steps = 0u;
        while (!pending.empty()) {
            const auto entry_index = pending.front();
            pending.pop_front();
            auto result = advance(history.object(node).entries.at(entry_index).work, statistics);
            count(result);
            for (const auto& successor: result.successors) {
                if (successor.kind == HybridEvolutionSuccessorKind::SUBDIVISION) {
                    whole_trajectory = false;
                }
            }
            for (std::size_t ordinal = 0; ordinal != result.successors.size(); ++ordinal) {
                const auto& successor = result.successors[ordinal];
                const HybridEvolutionEntryOrigin origin{entry_index, ordinal};
                if (derives_trajectory(successor, result.successors.size(), whole_trajectory)) {
                    HybridEvolutionTrajectory child;
                    child.entries.push_back({successor.work, std::nullopt, std::nullopt});
                    HybridEvolutionTrajectoryMetadata metadata;
                    metadata.origin = origin;
                    history.add_child(node, std::move(child), std::move(metadata));
                    // After an urgent branch, a later complete crossing concerns
                    // the remaining piece, not the original undivided trajectory.
                    if (successor.event_kind != EventKind::PERMISSIVE) {
                        whole_trajectory = false;
                    }
                } else {
                    auto& entries = history.object(node).entries;
                    const auto next_index = entries.size();
                    entries.push_back({successor.work, origin, std::nullopt});
                    pending.push_back(next_index);
                }
            }
            // Appending continuations can relocate entries: reacquire by index.
            history.object(node).entries.at(entry_index).result.emplace(std::move(result));
            ++completed_steps;
            if (policy->enables_recombination() &&
                trigger->should_evaluate(history, node, HybridRecombinationEvaluationPoint::AFTER_STEP,
                                         completed_steps))
                schedule_recombination(node);
        }
        history.metadata(node).complete = true;
    }

    /** @brief Replace accepted sibling groups at safe work-publication boundaries.
     *  No child has started evolving before its parent completes; even during
     *  per-step evaluation, only unevolved child checkpoints are selected.
     */
    Void schedule_recombination(HybridEvolutionTree::NodeId parent) {
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
            // Replacements must provide complete provenance and reduce the
            // number of unevolved trajectories, not just change representations.
            if (result.outputs.empty() || result.outputs.size() >= group.size() ||
                !result.discarded_sources.empty() || !result.radius_limit_satisfied)
                continue;
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
                !complete_coverage)
                continue;
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
                history.add_child(parent, std::move(trajectory), std::move(metadata));
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
                    !compatible_frontier_checkpoints(history, group.front(), checkpoint))
                    throw std::logic_error("Frontier policy selected incompatible checkpoints");
                selected.push_back(checkpoint);
            }
        }
        for (const auto& group: groups) {
            // A merged enclosure can be subdivided again immediately by the
            // numerical evolver. Do not let those pieces recombine at the same
            // physical time: every source descending from a previous merge
            // must have made strictly positive progress since that merge.
            bool progressed = true;
            double merge_time = -std::numeric_limits<double>::infinity();
            for (const auto checkpoint: group) {
                const auto& source_metadata = history.metadata(checkpoint.trajectory);
                const auto time = history.object(checkpoint.trajectory).entries.at(checkpoint.entry)
                                      .work.enclosure.time_range();
                const double lower = time.lower_bound().raw().get_d();
                const double upper = time.upper_bound().raw().get_d();
                if (!std::isfinite(lower) || !std::isfinite(upper) ||
                    (source_metadata.recombination_time_barrier &&
                     !(lower > *source_metadata.recombination_time_barrier)))
                    progressed = false;
                merge_time = std::max(merge_time, upper);
            }
            if (!progressed)
                continue;
            List<Enclosure> sources;
            for (const auto checkpoint: group)
                sources.append(history.object(checkpoint.trajectory).entries.at(checkpoint.entry)
                                   .work.enclosure.continuous_set());
            const auto result = recombiner->recombine(sources, context);
            if (result.outputs.empty() || result.outputs.size() >= group.size() ||
                !result.discarded_sources.empty() || !result.radius_limit_satisfied)
                continue;
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
                !std::all_of(represented.begin(), represented.end(), [](bool covered) { return covered; }))
                continue;
            if (!frontier_acceptance->accept(history, group, result))
                continue;
            const auto parent = common_ancestor(group);
            for (const auto& output: result.outputs) {
                const auto representative = group[output.source_indices.front()];
                auto work = history.object(representative.trajectory)
                                .entries.at(representative.entry).work;
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

    /** @brief Interleave work items from all branches, including split enclosures. */
    Void run_frontier(const StepFunction& advance) {
        const auto root = history.root().value();
        std::deque<HybridEvolutionCheckpointId> ready{{root, 0u}};
        std::vector<std::pair<HybridEvolutionCheckpointId, std::size_t>> deferrals;
        while (!ready.empty()) {
            const auto checkpoint = ready.front();
            ready.pop_front();
            const auto& retired = history.object(checkpoint.trajectory).retired_entries;
            if (std::find(retired.begin(), retired.end(), checkpoint.entry) != retired.end())
                continue;
            if (history.object(checkpoint.trajectory).entries.at(checkpoint.entry).result)
                throw std::logic_error("Checkpoint cannot be executed twice");
            // Only defer when a different, still-runnable work item is queued.
            // Keeping the per-item count across queue rotations guarantees
            // bounded waiting even when all candidates request deferral.
            const auto deferred = std::find_if(deferrals.begin(), deferrals.end(),
                [&](const auto& entry) { return entry.first == checkpoint; });
            const auto previous_deferrals = deferred == deferrals.end() ? 0u : deferred->second;
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
                frontier_dispatch->should_defer(history, checkpoint, alternatives, previous_deferrals)) {
                if (deferred == deferrals.end())
                    deferrals.push_back({checkpoint, 1u});
                else
                    ++deferred->second;
                ready.push_back(checkpoint);
                continue;
            }
            if (deferred != deferrals.end())
                deferrals.erase(deferred);
            auto result = advance(history.object(checkpoint.trajectory).entries.at(checkpoint.entry).work,
                                  statistics);
            count(result);
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
                    const auto created = history.add_child(checkpoint.trajectory, std::move(child),
                                                           std::move(metadata));
                    ready.push_back({created, 0u});
                    if (successor.kind == HybridEvolutionSuccessorKind::TRANSITION &&
                        successor.event_kind != EventKind::PERMISSIVE)
                        history.metadata(checkpoint.trajectory).whole_trajectory = false;
                } else {
                    auto& entries = history.object(checkpoint.trajectory).entries;
                    const auto index = entries.size();
                    entries.push_back({successor.work, origin, std::nullopt});
                    ++history.metadata(checkpoint.trajectory).pending_entries;
                    ready.push_back({checkpoint.trajectory, index});
                }
            }
            history.object(checkpoint.trajectory).entries.at(checkpoint.entry).result.emplace(std::move(result));
            auto& metadata = history.metadata(checkpoint.trajectory);
            if (metadata.pending_entries == 0u)
                throw std::logic_error("Invalid checkpoint completion count");
            --metadata.pending_entries;
            const bool trajectory_completed = metadata.pending_entries == 0u;
            if (trajectory_completed)
                metadata.complete = true;
            ++global_steps;
            if (frontier_policy->enables_recombination() &&
                trigger->should_evaluate(history, checkpoint.trajectory,
                    HybridRecombinationEvaluationPoint::AFTER_STEP, global_steps))
                schedule_frontier_recombination(ready);
            if (trajectory_completed && frontier_policy->enables_recombination() &&
                trigger->should_evaluate(history, checkpoint.trajectory,
                    HybridRecombinationEvaluationPoint::PARENT_COMPLETE, global_steps))
                schedule_frontier_recombination(ready);
        }
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
    if (_impl->use_frontier) {
        _impl->run_frontier(advance);
        _impl->complete = true;
        return;
    }
    std::deque<HybridEvolutionTree::NodeId> ready{_impl->history.root().value()};
    while (!ready.empty()) {
        const auto node = ready.front();
        ready.pop_front();
        _impl->process_trajectory(node, advance);
        // The trigger can request a final evaluation for any pending children.
        if (_impl->policy->enables_recombination() &&
            _impl->trigger->should_evaluate(_impl->history, node,
                HybridRecombinationEvaluationPoint::PARENT_COMPLETE,
                _impl->history.object(node).entries.size()))
            _impl->schedule_recombination(node);
        for (const auto child: _impl->history.children(node)) {
            if (!_impl->history.metadata(child).retired_for_recombination)
                ready.push_back(child);
        }
    }
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
    if (_frontier_grouping)
        return std::make_unique<TreeHybridEvolutionBackend>(
            initial, *_frontier_grouping, *_trigger, *_recombiner, _context,
            *_frontier_dispatch, *_frontier_acceptance);
    return std::make_unique<TreeHybridEvolutionBackend>(
        initial, *_grouping, *_trigger, *_recombiner, _context);
}

} // namespace Ariadne
