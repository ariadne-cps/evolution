/***************************************************************************
 *            hybrid_evolution_backend.cpp
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

/** @file src/hybrid_evolution_backend.cpp
 *  @brief Flat hybrid execution backend and step-local result recording.
 */

#include "evolution/hybrid_evolution_backend.hpp"

#include <algorithm>
#include <mutex>
#include <stdexcept>
#include <utility>

#include "threading/workload.hpp"

namespace Ariadne {

Void HybridEvolutionStepResult::record_reach(const HybridEnclosure& set) {
    contributions.push_back({set, HybridEvolutionSetKind::REACH});
}

Void HybridEvolutionStepResult::record_intermediate(const HybridEnclosure& set) {
    contributions.push_back({set, HybridEvolutionSetKind::INTERMEDIATE});
}

Void HybridEvolutionStepResult::record_final(const HybridEnclosure& set) {
    contributions.push_back({set, HybridEvolutionSetKind::FINAL});
}

Void HybridEvolutionStepResult::continue_with(const HybridEnclosure& set, HybridEvolutionSuccessorKind kind) {
    successors.push_back({
        {set, HybridEvolutionPhase::EVOLVE_CONTINUOUSLY},
        kind,
        std::nullopt,
        std::nullopt,
        std::nullopt
    });
}

Void HybridEvolutionStepResult::transition_to(const HybridEnclosure& set, DiscreteEvent event, EventKind event_kind,
                                              std::optional<CrossingKind> crossing_kind,
                                              HybridEvolutionTransitionResolution resolution) {
    successors.push_back({
        {set, HybridEvolutionPhase::PROCESS_STARTING_EVENTS},
        HybridEvolutionSuccessorKind::TRANSITION,
        event,
        event_kind,
        crossing_kind,
        resolution
    });
}

namespace {

/** @brief Per-run flat state; numerical callbacks never see its mutable Orbit. */
class FlatHybridEvolutionBackend final : public HybridEvolutionBackendInterface {
    std::unique_ptr<Orbit<HybridEnclosure>> _orbit;
    double _progress_start;
    double _progress_duration;
    std::mutex _mutex;
    Bool _started = false;
    Bool _complete = false;

    /** @brief Read published counters under the same lock used for publication. */
    HybridEvolutionStatistics statistics() {
        std::lock_guard<std::mutex> lock(_mutex);
        return {_orbit->reach().size(), _orbit->intermediate().size(), _orbit->final().size()};
    }

    /** @brief Publish contributions before making any successor eligible to run. */
    Void publish(const HybridEvolutionStepResult& result) {
        std::lock_guard<std::mutex> lock(_mutex);
        for (const auto& contribution: result.contributions) {
            switch (contribution.kind) {
                case HybridEvolutionSetKind::REACH:
                    _orbit->adjoin_reach(contribution.enclosure);
                    break;
                case HybridEvolutionSetKind::INTERMEDIATE:
                    _orbit->adjoin_intermediate(contribution.enclosure);
                    break;
                case HybridEvolutionSetKind::FINAL:
                    _orbit->adjoin_final(contribution.enclosure);
                    break;
                default:
                    throw std::invalid_argument("Unknown hybrid evolution contribution kind");
            }
        }
    }

  public:
    /** @brief Preserve the input and the existing progress-indicator time scale. */
    FlatHybridEvolutionBackend(const HybridEnclosure& initial, const Real& final_time)
        : _orbit(std::make_unique<Orbit<HybridEnclosure>>(initial)),
          _progress_start(initial.time_range().lower_bound().get_d()),
          _progress_duration(final_time.get_d() - _progress_start) {}

    /** @brief Reuse DynamicWorkload execution while committing complete local results. */
    Void run(const StepFunction& advance) override {
        if (_started) {
            throw std::logic_error("A hybrid evolution backend can only run once");
        }
        _started = true;
        using Workload = DynamicWorkload<HybridEvolutionWorkItem>;
        Workload workload(
            [this](const HybridEvolutionWorkItem& item, SharedPointer<ProgressIndicator> indicator) {
                if (_progress_duration > 0.0) {
                    const double elapsed = item.enclosure.time_range().lower_bound().get_d() - _progress_start;
                    indicator->update_current(std::clamp(elapsed, 0.0, _progress_duration));
                    indicator->update_final(_progress_duration);
                } else {
                    indicator->update_current(1.0);
                    indicator->update_final(1.0);
                }
            },
            [this, &advance](Workload::Access& pending, const HybridEvolutionWorkItem& item) {
                const auto result = advance(item, statistics());
                publish(result);
                // No lock is held while submitting potentially concurrent work.
                for (const auto& successor: result.successors) {
                    pending.append(successor.work);
                }
            });
        workload.append({_orbit->initial(), HybridEvolutionPhase::PROCESS_STARTING_EVENTS}).process();
        _complete = true;
    }

    /** @brief Transfer only a successfully completed result, preventing repeated release. */
    Orbit<HybridEnclosure> release_orbit() override {
        if (!_complete || !_orbit) {
            throw std::logic_error("No completed hybrid Orbit is available");
        }
        auto result = std::move(_orbit);
        return std::move(*result);
    }
};

} // namespace

FlatHybridEvolutionBackendFactory* FlatHybridEvolutionBackendFactory::clone() const {
    return new FlatHybridEvolutionBackendFactory(*this);
}

std::unique_ptr<HybridEvolutionBackendInterface> FlatHybridEvolutionBackendFactory::create(const HybridEnclosure& initial, const Real& final_time) const {
    return std::make_unique<FlatHybridEvolutionBackend>(initial, final_time);
}

} // namespace Ariadne
