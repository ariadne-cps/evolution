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
 *  @brief Policy-driven checkpoint scheduling and history-based Orbit projection.
 */

#include "evolution/hybrid_evolution_tree.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <iterator>
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

BeforeAdmissionRecombinationTriggerPolicy* BeforeAdmissionRecombinationTriggerPolicy::clone() const {
    return new BeforeAdmissionRecombinationTriggerPolicy(*this);
}

bool BeforeAdmissionRecombinationTriggerPolicy::should_evaluate(
    const HybridEvolutionTree&, NodeId, HybridRecombinationEvaluationPoint point, std::size_t, HybridRecombinationQueue) const {
    return point == HybridRecombinationEvaluationPoint::BEFORE_ADMISSION;
}

AfterEachStepRecombinationTriggerPolicy* AfterEachStepRecombinationTriggerPolicy::clone() const {
    return new AfterEachStepRecombinationTriggerPolicy(*this);
}

bool AfterEachStepRecombinationTriggerPolicy::should_evaluate(
    const HybridEvolutionTree&, NodeId, HybridRecombinationEvaluationPoint point, std::size_t, HybridRecombinationQueue) const {
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
    std::size_t completed_steps, HybridRecombinationQueue queue) const {
    return (point == HybridRecombinationEvaluationPoint::BEFORE_ADMISSION && queue == HybridRecombinationQueue::WAITING) ||
           (point == HybridRecombinationEvaluationPoint::AFTER_STEP && completed_steps % _period == 0u);
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

// Follow replacement provenance to retain the original successor meaning.
const HybridEvolutionSuccessor* seed_successor(const HybridEvolutionTree& tree,
                                               HybridEvolutionCheckpointId checkpoint) {
    while (!tree.metadata(checkpoint.trajectory).recombined_checkpoints.empty())
        checkpoint = tree.metadata(checkpoint.trajectory).recombined_checkpoints.front();
    const auto& metadata = tree.metadata(checkpoint.trajectory);
    const auto parent = tree.parent(checkpoint.trajectory);
    if (checkpoint.entry != 0u || !metadata.origin || !parent)
        return nullptr;
    return &tree.object(*parent).entries.at(metadata.origin->entry).result.value().successors.at(metadata.origin->successor);
}

bool compatible_scope(const HybridEvolutionTree& tree, HybridEvolutionCheckpointId first,
                      HybridEvolutionCheckpointId next, HybridRecombinationScope scope) {
    if (scope == HybridRecombinationScope::ANY)
        return true;
    if (first.entry != 0u || next.entry != 0u ||
        tree.parent(first.trajectory) != tree.parent(next.trajectory))
        return false;
    const auto* a = seed_successor(tree, first);
    const auto* b = seed_successor(tree, next);
    return a && b && a->kind == b->kind && a->event == b->event && a->event_kind == b->event_kind;
}

// A continuous family contains a seed and all descendants born by subdivision,
// but excludes discrete descendants, which wait for that whole family.
HybridEvolutionTree::NodeId continuous_family(const HybridEvolutionTree& tree,
                                              HybridEvolutionTree::NodeId node) {
    while (const auto parent = tree.parent(node)) {
        const auto& origin = tree.metadata(node).origin;
        if (!origin || tree.object(*parent).entries.at(origin->entry).result.value().successors.at(origin->successor).kind !=
                           HybridEvolutionSuccessorKind::SUBDIVISION)
            break;
        node = *parent;
    }
    return node;
}

bool family_complete(const HybridEvolutionTree& tree, HybridEvolutionTree::NodeId root) {
    std::vector<HybridEvolutionTree::NodeId> pending{root};
    while (!pending.empty()) {
        const auto node = pending.back();
        pending.pop_back();
        if (!tree.metadata(node).complete)
            return false;
        for (const auto child: tree.children(node)) {
            const auto& origin = tree.metadata(child).origin;
            if (origin && tree.object(node).entries.at(origin->entry).result.value().successors.at(origin->successor).kind ==
                              HybridEvolutionSuccessorKind::SUBDIVISION)
                pending.push_back(child);
        }
    }
    return true;
}

// Same numerical and discrete meaning is required for future work. A
// validated recombiner combines the clock as part of the continuous set;
// overlapping time bounds are an additional conservative scheduler filter.
bool compatible_checkpoints(const HybridEvolutionTree& tree,
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

AllCompatibleRecombinationPolicy::AllCompatibleRecombinationPolicy(HybridRecombinationScope scope)
    : _scope(scope) {}

AllCompatibleRecombinationPolicy*
    AllCompatibleRecombinationPolicy::clone() const {
    return new AllCompatibleRecombinationPolicy(*this);
}

std::vector<std::vector<HybridEvolutionCheckpointId>>
    AllCompatibleRecombinationPolicy::select(
        const HybridEvolutionTree& tree, const std::vector<CheckpointId>& pending) const {
    std::vector<std::vector<CheckpointId>> groups;
    for (const auto checkpoint: pending) {
        const auto group = std::find_if(groups.begin(), groups.end(), [&](const auto& candidates) {
            return compatible_scope(tree, candidates.front(), checkpoint, _scope) &&
                   compatible_checkpoints(tree, candidates.front(), checkpoint, false);
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

ParentCheckpointAdmissionPolicy* ParentCheckpointAdmissionPolicy::clone() const {
    return new ParentCheckpointAdmissionPolicy(*this);
}

bool ParentCheckpointAdmissionPolicy::admit_successor(
    const HybridEvolutionWorkItem&, const HybridEvolutionSuccessor& successor, bool derived) const {
    return !derived || successor.kind != HybridEvolutionSuccessorKind::TRANSITION;
}

std::vector<HybridEvolutionCheckpointId> ParentCheckpointAdmissionPolicy::select_next(
    const HybridEvolutionTree& tree, const std::vector<CheckpointId>& waiting, bool) const {
    std::vector<CheckpointId> selected;
    std::vector<std::pair<HybridEvolutionTree::NodeId, bool>> families;
    for (const auto checkpoint: waiting) {
        bool eligible = true;
        std::vector<CheckpointId> sources{checkpoint};
        while (!sources.empty() && eligible) {
            const auto source = sources.back();
            sources.pop_back();
            const auto& provenance = tree.metadata(source.trajectory).recombined_checkpoints;
            if (!provenance.empty()) {
                sources.insert(sources.end(), provenance.begin(), provenance.end());
                continue;
            }
            const auto parent = tree.parent(source.trajectory);
            if (!parent)
                throw std::logic_error("Waiting descendant has no parent");
            const auto family = continuous_family(tree, *parent);
            auto found = std::find_if(families.begin(), families.end(), [&](const auto& item) { return item.first == family; });
            if (found == families.end()) {
                families.push_back({family, family_complete(tree, family)});
                found = std::prev(families.end());
            }
            eligible = found->second;
        }
        if (eligible)
            selected.push_back(checkpoint);
    }
    return selected;
}

ImmediateCheckpointAdmissionPolicy* ImmediateCheckpointAdmissionPolicy::clone() const {
    return new ImmediateCheckpointAdmissionPolicy(*this);
}

bool ImmediateCheckpointAdmissionPolicy::admit_successor(
    const HybridEvolutionWorkItem&, const HybridEvolutionSuccessor&, bool) const { return true; }

std::vector<HybridEvolutionCheckpointId> ImmediateCheckpointAdmissionPolicy::select_next(
    const HybridEvolutionTree&, const std::vector<CheckpointId>& waiting, bool) const { return waiting; }

LocationCheckpointAdmissionPolicy* LocationCheckpointAdmissionPolicy::clone() const {
    return new LocationCheckpointAdmissionPolicy(*this);
}

bool LocationCheckpointAdmissionPolicy::admit_successor(
    const HybridEvolutionWorkItem& producer, const HybridEvolutionSuccessor& successor, bool) const {
    return successor.kind != HybridEvolutionSuccessorKind::TRANSITION &&
           successor.work.enclosure.location() == producer.enclosure.location();
}

std::vector<HybridEvolutionCheckpointId> LocationCheckpointAdmissionPolicy::select_next(
    const HybridEvolutionTree& tree, const std::vector<CheckpointId>& waiting, bool batch_complete) const {
    std::vector<CheckpointId> selected;
    if (!batch_complete)
        return selected;
    if (waiting.empty())
        return selected;
    const auto location = tree.object(waiting.front().trajectory).entries.at(waiting.front().entry).work.enclosure.location();
    for (const auto checkpoint: waiting)
        if (tree.object(checkpoint.trajectory).entries.at(checkpoint.entry).work.enclosure.location() == location)
            selected.push_back(checkpoint);
    return selected;
}

NoRecombinationPolicy* NoRecombinationPolicy::clone() const {
    return new NoRecombinationPolicy(*this);
}

std::vector<std::vector<HybridEvolutionCheckpointId>>
    NoRecombinationPolicy::select(
        const HybridEvolutionTree&, const std::vector<CheckpointId>&) const {
    return {};
}

WindowRecombinationPolicy::WindowRecombinationPolicy(
    std::size_t max_group_size, double max_time_spread, HybridRecombinationScope scope)
    : _max_group_size(max_group_size), _max_time_spread(max_time_spread), _scope(scope) {
    if (_max_group_size < 2u || !std::isfinite(_max_time_spread) || _max_time_spread < 0.0)
        throw std::invalid_argument("Window recombination requires group size >= 2 and finite nonnegative time spread");
}

WindowRecombinationPolicy* WindowRecombinationPolicy::clone() const {
    return new WindowRecombinationPolicy(*this);
}

std::vector<std::vector<HybridEvolutionCheckpointId>>
    WindowRecombinationPolicy::select(
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
            if (used[next])
                continue;
            if (_scope == HybridRecombinationScope::SIBLINGS &&
                tree.parent(pending[first].trajectory) != tree.parent(pending[next].trajectory))
                continue;
            if (!compatible_scope(tree, pending[first], pending[next], _scope) ||
                !compatible_checkpoints(tree, pending[first], pending[next], requires_time_overlap())) {
                if (_scope == HybridRecombinationScope::SIBLINGS)
                    break;
                continue;
            }
            const auto time = tree.object(pending[next].trajectory).entries.at(pending[next].entry).work.enclosure.time_range();
            const double candidate_lower = time.lower_bound().raw().get_d();
            const double candidate_upper = time.upper_bound().raw().get_d();
            if (!std::isfinite(candidate_lower) || !std::isfinite(candidate_upper) ||
                std::max(upper, candidate_upper) - std::min(lower, candidate_lower) > _max_time_spread) {
                if (_scope == HybridRecombinationScope::SIBLINGS)
                    break;
                continue;
            }
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

ImmediateCheckpointDispatchPolicy* ImmediateCheckpointDispatchPolicy::clone() const {
    return new ImmediateCheckpointDispatchPolicy(*this);
}

bool ImmediateCheckpointDispatchPolicy::should_defer(
    const HybridEvolutionTree&, CheckpointId, const std::vector<CheckpointId>&,
    std::size_t) const {
    return false;
}

BoundedTimeLagCheckpointDispatchPolicy::BoundedTimeLagCheckpointDispatchPolicy(
    std::size_t max_deferrals, double max_time_lag)
    : _max_deferrals(max_deferrals), _max_time_lag(max_time_lag) {
    if (!std::isfinite(_max_time_lag) || _max_time_lag < 0.0)
        throw std::invalid_argument("Checkpoint dispatch requires a finite nonnegative time lag");
}

BoundedTimeLagCheckpointDispatchPolicy* BoundedTimeLagCheckpointDispatchPolicy::clone() const {
    return new BoundedTimeLagCheckpointDispatchPolicy(*this);
}

bool BoundedTimeLagCheckpointDispatchPolicy::should_defer(
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

AcceptAllRecombinationsPolicy* AcceptAllRecombinationsPolicy::clone() const {
    return new AcceptAllRecombinationsPolicy(*this);
}

bool AcceptAllRecombinationsPolicy::accept(
    const HybridEvolutionTree&, const std::vector<CheckpointId>&,
    const EnclosureRecombinationResult&) const {
    return true;
}

RadiusTrendRecombinationAcceptancePolicy::RadiusTrendRecombinationAcceptancePolicy(
    double max_output_amplification, double max_continuation_growth)
    : _max_output_amplification(max_output_amplification),
      _max_continuation_growth(max_continuation_growth) {
    if (!std::isfinite(_max_output_amplification) || _max_output_amplification <= 0.0 ||
        !std::isfinite(_max_continuation_growth) || _max_continuation_growth <= 0.0)
        throw std::invalid_argument("Radius-trend merge limits must be finite and strictly positive");
}

RadiusTrendRecombinationAcceptancePolicy* RadiusTrendRecombinationAcceptancePolicy::clone() const {
    return new RadiusTrendRecombinationAcceptancePolicy(*this);
}

bool RadiusTrendRecombinationAcceptancePolicy::accept(
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

/** @brief Run state, kept out of the installed backend class's representation. */
class TreeHybridEvolutionBackend::Impl {
  public:
    HybridEvolutionTree history;
    HybridEvolutionStatistics statistics;
    HybridEvolutionProgressTracker progress_tracker;
    mutable std::mutex schedule_mutex;
    std::unique_ptr<const HybridRecombinationGroupingPolicy> grouping_policy;
    std::unique_ptr<const HybridCheckpointDispatchPolicy> dispatch_policy;
    std::unique_ptr<const HybridRecombinationAcceptancePolicy> acceptance_policy;
    std::unique_ptr<const HybridCheckpointAdmissionPolicy> admission;
    std::unique_ptr<const HybridRecombinationTriggerPolicy> trigger;
    std::size_t global_steps = 0u;
    std::unique_ptr<const EnclosureRecombinerInterface> recombiner;
    EnclosureRecombinationContext context;
    Semantics semantics = Semantics::LOWER;
    Bool started = false;
    Bool complete = false;
    Bool released = false;

    /** @brief Initialize checkpoint history and independent scheduling policies. */
    explicit Impl(const HybridEnclosure& initial,
                  const HybridCheckpointAdmissionPolicy& admission_policy,
                  const HybridRecombinationGroupingPolicy& grouping,
                  const HybridRecombinationTriggerPolicy& timing,
                  const EnclosureRecombinerInterface& method,
                  const EnclosureRecombinationContext& limits,
                  const HybridCheckpointDispatchPolicy& dispatch,
                  const HybridRecombinationAcceptancePolicy& acceptance)
        : grouping_policy(grouping.clone()), dispatch_policy(dispatch.clone()),
          acceptance_policy(acceptance.clone()),
          admission(admission_policy.clone()), trigger(timing.clone()),
          recombiner(method.clone()), context(limits) {
        if (grouping.enables_recombination() && !acceptance.supports_queue(admission->recombination_queue()))
            throw std::invalid_argument("Merge acceptance policy does not support the admission policy pending queue");
        HybridEvolutionTrajectory trajectory;
        trajectory.entries.push_back({
            {initial, HybridEvolutionPhase::PROCESS_STARTING_EVENTS},
            std::nullopt,
            std::nullopt
        });
        HybridEvolutionTrajectoryMetadata metadata;
        metadata.pending_entries = 1u;
        history.create_root(std::move(trajectory), std::move(metadata));
        progress_tracker.add_pending(initial);
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
        throw std::logic_error("Checkpoints have no common genealogy root");
    }

    /** @brief Recombine compatible pending checkpoints from any part of the tree.
     *  A processed entry is never replaced; all selected source IDs remain in
     *  append-only history with their original result-less checkpoint records.
     */
    Void schedule_recombination(std::deque<HybridEvolutionCheckpointId>& ready) {
        LOGGING_SCOPE_CREATE;
        if (!grouping_policy->enables_recombination())
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
        LOGGING_PRINTLN("RECOMBINE inspect"
                        << " progress=" << 100.0 * progress_tracker.snapshot().completion_rate()
                        << " pending=" << pending.size()
                        << " queue=" << (admission->recombination_queue() == HybridRecombinationQueue::READY ? "ready" : "waiting"));
        const auto groups = grouping_policy->select(history, pending);
        if (groups.empty())
            LOGGING_PRINTLN("RECOMBINE skipped: no compatible groups");
        std::vector<HybridEvolutionCheckpointId> selected;
        for (const auto& group: groups) {
            if (group.size() < 2u)
                throw std::logic_error("Grouping policy selected fewer than two checkpoints");
            for (const auto checkpoint: group) {
                if (!history.contains(checkpoint.trajectory) ||
                    std::find(pending.begin(), pending.end(), checkpoint) == pending.end() ||
                    std::find(selected.begin(), selected.end(), checkpoint) != selected.end())
                    throw std::logic_error("Grouping policy selected unknown or overlapping checkpoint");
                if (checkpoint != group.front() &&
                    !compatible_checkpoints(history, group.front(), checkpoint,
                                            grouping_policy->requires_time_overlap()))
                    throw std::logic_error("Grouping policy selected incompatible checkpoints");
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
            double earliest_time = std::numeric_limits<double>::infinity();
            SizeType max_parameters = 0u;
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
                earliest_time = std::min(earliest_time, lower);
                max_parameters = std::max(max_parameters,
                                          history.object(checkpoint.trajectory).entries.at(checkpoint.entry).work.enclosure.number_of_parameters());
            }
            if (!progressed) {
                LOGGING_PRINTLN("RECOMBINE skipped: no time progress"
                                << " sources=" << group.size());
                continue;
            }
            List<Enclosure> sources;
            for (const auto checkpoint: group)
                sources.append(history.object(checkpoint.trajectory).entries.at(checkpoint.entry).work.enclosure.continuous_set());
            const auto location = history.object(group.front().trajectory).entries.at(group.front().entry).work.enclosure.location();
            const auto source_time = history.object(group.front().trajectory).entries.at(group.front().entry).work.enclosure.time_range();
            LOGGING_PRINTLN("RECOMBINE begin"
                            << " progress=" << 100.0 * progress_tracker.snapshot().completion_rate()
                            << " location=" << location
                            << " sources=" << group.size()
                            << " time=[" << earliest_time << "," << merge_time << "]"
                            << " max_parameters=" << max_parameters);
            const auto merge_start = std::chrono::steady_clock::now();
            EnclosureRecombinationResult result;
            try {
                result = recombiner->recombine(sources, context);
            } catch (...) {
                LOGGING_PRINTLN("RECOMBINE failed"
                                << " location=" << location
                                << " sources=" << group.size()
                                << " elapsed_seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - merge_start).count());
                throw;
            }
            const double elapsed_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - merge_start).count();
            if (result.outputs.empty() || result.outputs.size() >= group.size() ||
                !result.discarded_sources.empty() || !result.radius_limit_satisfied) {
                LOGGING_PRINTLN("RECOMBINE rejected"
                                << " sources=" << group.size()
                                << " outputs=" << result.outputs.size()
                                << " discarded=" << result.discarded_sources.size()
                                << " radius_ok=" << result.radius_limit_satisfied
                                << " elapsed_seconds=" << elapsed_seconds);
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
                LOGGING_PRINTLN("RECOMBINE rejected: incomplete provenance"
                                << " sources=" << group.size() << " outputs=" << result.outputs.size()
                                << " elapsed_seconds=" << elapsed_seconds);
                continue;
            }
            if (!acceptance_policy->accept(history, group, result)) {
                LOGGING_PRINTLN("RECOMBINE rejected: quality gate"
                                << " sources=" << group.size() << " outputs=" << result.outputs.size()
                                << " elapsed_seconds=" << elapsed_seconds);
                continue;
            }
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
                progress_tracker.add_pending(history.object(created).entries.front().work.enclosure);
            }
            for (const auto checkpoint: group) {
                progress_tracker.remove_pending(history.object(checkpoint.trajectory).entries.at(checkpoint.entry).work.enclosure);
                history.object(checkpoint.trajectory).retired_entries.push_back(checkpoint.entry);
                auto& metadata = history.metadata(checkpoint.trajectory);
                if (metadata.pending_entries == 0u)
                    throw std::logic_error("Invalid checkpoint pending count");
                --metadata.pending_entries;
                if (metadata.pending_entries == 0u)
                    metadata.complete = true;
            }
            LOGGING_PRINTLN("RECOMBINE accepted"
                            << " progress=" << 100.0 * progress_tracker.snapshot().completion_rate()
                            << " location=" << location
                            << " t=" << source_time
                            << " sources=" << group.size() << " outputs=" << result.outputs.size()
                            << " max_source_radius=" << largest_source_radius(sources)
                            << " max_output_radius=" << largest_output_radius(result)
                            << " output_time=" << result.outputs.front().enclosure.time_function().range()
                            << " elapsed_seconds=" << elapsed_seconds);
        }
    }

    /** @brief Execute admissible checkpoints with serialized result publication.
     *  @details The same scheduler handles zero or several worker threads.
     *  Every ready checkpoint has one DynamicWorkload token; the token
     *  takes a checkpoint from the protected queue before beginning numerical
     *  work. All admission policies share ready, waiting and in-flight state.
     *  Only unprocessed pending checkpoints can be grouped or retired; a
     *  numerical callback never holds the publication lock.
     */
    Void run_scheduled(const StepFunction& advance) {
        struct Token {};
        using Workload = DynamicWorkload<Token>;
        const auto root = history.root().value();
        std::deque<HybridEvolutionCheckpointId> ready{
            {root, 0u}
        };
        std::deque<HybridEvolutionCheckpointId> waiting;
        std::size_t in_flight = 0u;
        std::vector<std::pair<HybridEvolutionCheckpointId, std::size_t>> deferrals;
        Workload workload(
            [this](const Token&, SharedPointer<ProgressIndicator> indicator) {
                std::lock_guard<std::mutex> lock(schedule_mutex);
                const auto state = progress_tracker.snapshot();
                indicator->update_current(state.complete ? 1.0 : std::min(0.99, state.completion_rate()));
                indicator->update_final(1.0);
            },
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

                        if (dispatch_policy->enables_deferral()) {
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
                                dispatch_policy->should_defer(history, candidate, alternatives, previous)) {
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
                const auto evaluate_recombination = [&] {
                    const bool ready_queue = admission->recombination_queue() == HybridRecombinationQueue::READY;
                    auto& candidates = ready_queue ? ready : waiting;
                    const auto previous_size = candidates.size();
                    schedule_recombination(candidates);
                    // Recombined replacements append new ready checkpoints;
                    // existing retired tokens may still be in the workload.
                    if (ready_queue) {
                        for (std::size_t i = previous_size; i != candidates.size(); ++i)
                            producer.append(Token{});
                        ready.erase(std::remove_if(ready.begin(), ready.end(), is_retired), ready.end());
                    }
                };

                count(result);
                progress_tracker.publish(*work, result);
                --in_flight;
                {
                    const auto publish_successor = [&](HybridEvolutionCheckpointId next,
                                                       const HybridEvolutionSuccessor& successor) {
                        if (admission->admit_successor(*work, successor, next.trajectory != checkpoint.trajectory))
                            enqueue(next);
                        else
                            waiting.push_back(next);
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
                    if (metadata.pending_entries == 0u)
                        metadata.complete = true;
                    ++global_steps;
                    bool evaluated = false;
                    if (grouping_policy->enables_recombination() &&
                        trigger->should_evaluate(history, checkpoint.trajectory,
                                                 HybridRecombinationEvaluationPoint::AFTER_STEP, global_steps, admission->recombination_queue())) {
                        evaluate_recombination();
                        evaluated = true;
                    }
                    const bool batch_complete = ready.empty() && in_flight == 0u;
                    const auto select_admitted = [&] {
                        std::vector<HybridEvolutionCheckpointId> pending;
                        for (const auto id: waiting)
                            if (!is_retired(id))
                                pending.push_back(id);
                        return admission->select_next(history, pending, batch_complete);
                    };
                    // Parent admission may release one completed family while
                    // other families are still ready or running. Location
                    // admission returns no eligible work until the batch drains.
                    auto admitted = select_admitted();
                    if (grouping_policy->enables_recombination() && !evaluated &&
                        (admission->recombination_queue() == HybridRecombinationQueue::READY || !admitted.empty()) &&
                        trigger->should_evaluate(history, checkpoint.trajectory,
                                                 HybridRecombinationEvaluationPoint::BEFORE_ADMISSION, global_steps, admission->recombination_queue())) {
                        evaluate_recombination();
                        admitted = select_admitted();
                    }
                    {
                        waiting.erase(std::remove_if(waiting.begin(), waiting.end(), is_retired),
                                      waiting.end());
                        const std::vector<HybridEvolutionCheckpointId> pending(waiting.begin(), waiting.end());
                        std::vector<HybridEvolutionCheckpointId> unique;
                        for (const auto next: admitted) {
                            if (std::find(pending.begin(), pending.end(), next) == pending.end() ||
                                std::find(unique.begin(), unique.end(), next) != unique.end())
                                throw std::logic_error("Admission policy selected unknown or duplicate waiting work");
                            unique.push_back(next);
                        }
                        if (batch_complete && !pending.empty() && admitted.empty())
                            throw std::logic_error("Admission policy left waiting work without a runnable batch");
                        if (!admitted.empty()) {
                            LOGGING_RUN_AT(1, const auto state = progress_tracker.snapshot());
                            LOGGING_PRINTLN_AT(1, "ADMIT"
                                                      << " progress=" << 100.0 * state.completion_rate()
                                                      << " location=" << history.object(admitted.front().trajectory).entries.at(admitted.front().entry).work.enclosure.location()
                                                      << " time=" << history.object(admitted.front().trajectory).entries.at(admitted.front().entry).work.enclosure.time_range()
                                                      << " checkpoints=" << admitted.size()
                                                      << " ready=" << ready.size() + admitted.size()
                                                      << " waiting=" << pending.size() - admitted.size()
                                                      << " active=" << in_flight
                                                      << " credited_time=" << state.credited_time
                                                      << " time_estimate_available=" << state.remaining_time.has_value()
                                                      << " remaining_time=" << (state.remaining_time ? *state.remaining_time : 0.0));
                        }
                        for (const auto next: admitted) {
                            waiting.erase(std::find(waiting.begin(), waiting.end(), next));
                            enqueue(next);
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
    const HybridRecombinationGroupingPolicy& grouping,
    const HybridRecombinationTriggerPolicy& timing, const EnclosureRecombinerInterface& method,
    const EnclosureRecombinationContext& limits, const HybridCheckpointDispatchPolicy& dispatch,
    const HybridRecombinationAcceptancePolicy& acceptance)
    : _impl(std::make_unique<Impl>(initial, admission, grouping, timing, method, limits, dispatch, acceptance)) {}

TreeHybridEvolutionBackend::~TreeHybridEvolutionBackend() = default;

Void TreeHybridEvolutionBackend::set_semantics(Semantics semantics) {
    if (_impl->started)
        throw std::logic_error("Cannot change recombination semantics after the run starts");
    _impl->semantics = semantics;
}

Void TreeHybridEvolutionBackend::set_time_horizon(const Real& final_time) {
    if (_impl->started)
        throw std::logic_error("Cannot change the progress horizon after the run starts");
    _impl->progress_tracker = HybridEvolutionProgressTracker(final_time.get_d());
    _impl->progress_tracker.add_pending(_impl->history.object(_impl->history.root().value()).entries.front().work.enclosure);
}

HybridEvolutionProgress TreeHybridEvolutionBackend::progress() const {
    std::lock_guard<std::mutex> lock(_impl->schedule_mutex);
    return _impl->progress_tracker.snapshot();
}

Void TreeHybridEvolutionBackend::run(const StepFunction& advance) {
    if (_impl->started) {
        throw std::logic_error("A hybrid evolution backend can only run once");
    }
    if (_impl->grouping_policy->enables_recombination() &&
        _impl->semantics != Semantics::UPPER)
        throw std::invalid_argument("Recombination is only sound for upper semantics");
    _impl->started = true;
    // The scheduler is identical for sequential and concurrent execution:
    // DynamicWorkload controls how many ready checkpoints run at once.
    _impl->run_scheduled(advance);
    std::lock_guard<std::mutex> lock(_impl->schedule_mutex);
    _impl->progress_tracker.finish();
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
        _impl->project_trajectory(node, result);
    }
    _impl->released = true;
    return result;
}

const HybridEvolutionTree& TreeHybridEvolutionBackend::tree() const noexcept {
    return _impl->history;
}

TreeHybridEvolutionBackendFactory::TreeHybridEvolutionBackendFactory(
    const HybridCheckpointAdmissionPolicy& admission, const HybridRecombinationGroupingPolicy& grouping,
    const HybridRecombinationTriggerPolicy& timing, const EnclosureRecombinerInterface& method,
    const EnclosureRecombinationContext& limits, const HybridCheckpointDispatchPolicy& dispatch,
    const HybridRecombinationAcceptancePolicy& acceptance)
    : _admission(admission.clone()), _grouping(grouping.clone()), _dispatch(dispatch.clone()),
      _acceptance(acceptance.clone()), _trigger(timing.clone()), _recombiner(method.clone()), _context(limits) {
    if (grouping.enables_recombination() && !acceptance.supports_queue(admission.recombination_queue()))
        throw std::invalid_argument("Merge acceptance policy does not support the admission policy pending queue");
}

TreeHybridEvolutionBackendFactory* TreeHybridEvolutionBackendFactory::clone() const {
    return new TreeHybridEvolutionBackendFactory(*this);
}

std::unique_ptr<HybridEvolutionBackendInterface> TreeHybridEvolutionBackendFactory::create(const HybridEnclosure& initial, const Real& final_time) const {
    auto backend = std::make_unique<TreeHybridEvolutionBackend>(
        initial, *_admission, *_grouping, *_trigger, *_recombiner, _context, *_dispatch, *_acceptance);
    backend->set_time_horizon(final_time);
    return backend;
}

} // namespace Ariadne
