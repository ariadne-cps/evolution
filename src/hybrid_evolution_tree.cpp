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

} // namespace

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


/** @brief Run state, kept out of the installed backend class's representation. */
class TreeHybridEvolutionBackend::Impl {
  public:
    HybridEvolutionTree history;
    HybridEvolutionStatistics statistics;
    std::unique_ptr<const HybridRecombinationSchedulingPolicy> policy;
    Semantics semantics = Semantics::LOWER;
    Bool started = false;
    Bool complete = false;
    Bool released = false;

    /** @brief Create a trajectory seed without an implicit Orbit contribution. */
    explicit Impl(const HybridEnclosure& initial, const HybridRecombinationSchedulingPolicy& scheduling_policy)
        : policy(scheduling_policy.clone()) {
        HybridEvolutionTrajectory trajectory;
        trajectory.entries.push_back({
            {initial, HybridEvolutionPhase::PROCESS_STARTING_EVENTS},
            std::nullopt,
            std::nullopt
        });
        history.create_root(std::move(trajectory), HybridEvolutionTrajectoryMetadata{});
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
        }
        history.metadata(node).complete = true;
    }

    /** @brief Replace accepted sibling groups while keeping their unevolved sources inspectable.
     *  This runs strictly after parent completion and before any child dispatch.
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
            const auto result = ReferenceEnclosureRecombiner().recombine(sources);
            // Reject an ineffective or incomplete replacement without losing any work.
            if (result.outputs.size() != 1u || !result.discarded_sources.empty() ||
                !result.radius_limit_satisfied ||
                result.outputs.front().source_indices.size() != group.size())
                continue;
            Bool complete_coverage = true;
            for (std::size_t index = 0; index != group.size(); ++index)
                if (result.outputs.front().source_indices[index] != index)
                    complete_coverage = false;
            if (!complete_coverage)
                continue;
            auto replacement = history.object(group.front()).entries.front().work;
            replacement.enclosure.replace_continuous_set(result.outputs.front().enclosure);
            HybridEvolutionTrajectory trajectory;
            trajectory.entries.push_back({replacement, std::nullopt, std::nullopt});
            HybridEvolutionTrajectoryMetadata metadata;
            metadata.recombined_sources = group;
            history.add_child(parent, std::move(trajectory), std::move(metadata));
            for (const auto source: group)
                history.metadata(source).retired_for_recombination = true;
        }
    }

    /** @brief Append one trajectory's explicit contributions in creation order. */
    Void project_trajectory(HybridEvolutionTree::NodeId node, Orbit<HybridEnclosure>& orbit) const {
        for (const auto& entry: history.object(node).entries) {
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
    const HybridEnclosure& initial, const HybridRecombinationSchedulingPolicy& policy)
    : _impl(std::make_unique<Impl>(initial, policy)) {}

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
    if (_impl->policy->enables_recombination() && _impl->semantics != Semantics::UPPER)
        throw std::invalid_argument("Sibling recombination is only sound for upper semantics");
    _impl->started = true;
    std::deque<HybridEvolutionTree::NodeId> ready{_impl->history.root().value()};
    while (!ready.empty()) {
        const auto node = ready.front();
        ready.pop_front();
        _impl->process_trajectory(node, advance);
        // Only after publishing completion may policy inspect all born children.
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
    : _policy(std::make_shared<NoRecombinationSchedulingPolicy>()) {}

TreeHybridEvolutionBackendFactory::TreeHybridEvolutionBackendFactory(
    const HybridRecombinationSchedulingPolicy& policy)
    : _policy(policy.clone()) {}

TreeHybridEvolutionBackendFactory* TreeHybridEvolutionBackendFactory::clone() const {
    return new TreeHybridEvolutionBackendFactory(*this);
}

std::unique_ptr<HybridEvolutionBackendInterface> TreeHybridEvolutionBackendFactory::create(const HybridEnclosure& initial, const Real&) const {
    return std::make_unique<TreeHybridEvolutionBackend>(initial, *_policy);
}

} // namespace Ariadne
