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

#include <deque>
#include <stdexcept>
#include <utility>

namespace Ariadne {

/** @brief Run state, kept out of the installed backend class's representation. */
class TreeHybridEvolutionBackend::Impl {
  public:
    HybridEvolutionTree history;
    HybridEvolutionStatistics statistics;
    Bool started = false;
    Bool complete = false;
    Bool released = false;

    /** @brief Create a trajectory seed without an implicit Orbit contribution. */
    explicit Impl(const HybridEnclosure& initial) {
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
                    history.add_child(node, std::move(child), HybridEvolutionTrajectoryMetadata{false, origin});
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

TreeHybridEvolutionBackend::TreeHybridEvolutionBackend(const HybridEnclosure& initial)
    : _impl(std::make_unique<Impl>(initial)) {}

TreeHybridEvolutionBackend::~TreeHybridEvolutionBackend() = default;

Void TreeHybridEvolutionBackend::run(const StepFunction& advance) {
    if (_impl->started) {
        throw std::logic_error("A hybrid evolution backend can only run once");
    }
    _impl->started = true;
    std::deque<HybridEvolutionTree::NodeId> ready{_impl->history.root().value()};
    while (!ready.empty()) {
        const auto node = ready.front();
        ready.pop_front();
        _impl->process_trajectory(node, advance);
        // Parent completion is published after all of its children have been generated.
        for (const auto child: _impl->history.children(node)) {
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
        _impl->project_trajectory(node, result);
    }
    _impl->released = true;
    return result;
}

const HybridEvolutionTree& TreeHybridEvolutionBackend::tree() const noexcept {
    return _impl->history;
}

TreeHybridEvolutionBackendFactory* TreeHybridEvolutionBackendFactory::clone() const {
    return new TreeHybridEvolutionBackendFactory(*this);
}

std::unique_ptr<HybridEvolutionBackendInterface> TreeHybridEvolutionBackendFactory::create(const HybridEnclosure& initial, const Real&) const {
    return std::make_unique<TreeHybridEvolutionBackend>(initial);
}

} // namespace Ariadne
