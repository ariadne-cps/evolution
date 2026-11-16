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
#include <cmath>
#include <limits>
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

double HybridEvolutionProgress::completion_rate() const {
    if (complete)
        return 1.0;
    if (!remaining_time || !std::isfinite(credited_time) || !std::isfinite(*remaining_time))
        return 0.0;
    const long double total = static_cast<long double>(credited_time) + *remaining_time;
    return total > 0.0L ? static_cast<double>(std::clamp(static_cast<long double>(credited_time) / total, 0.0L, 1.0L)) : 0.0;
}

HybridEvolutionProgressTracker::HybridEvolutionProgressTracker(std::optional<double> horizon)
    : _horizon(horizon && std::isfinite(*horizon) ? horizon : std::nullopt) {}

std::optional<long double> HybridEvolutionProgressTracker::remaining(const HybridEnclosure& enclosure) const {
    if (!_horizon)
        return std::nullopt;
    const double lower = enclosure.time_range().lower_bound().get_d();
    if (!std::isfinite(lower))
        return std::nullopt;
    return std::max(0.0L, static_cast<long double>(*_horizon) - lower);
}

std::optional<long double> HybridEvolutionProgressTracker::duration() const {
    if (!_horizon || !_initial_time)
        return std::nullopt;
    return std::max(0.0L, static_cast<long double>(*_horizon) - *_initial_time);
}

std::optional<long double> HybridEvolutionProgressTracker::elapsed(const HybridEnclosure& enclosure) const {
    const auto span = duration();
    if (!span)
        return std::nullopt;
    const double lower = enclosure.time_range().lower_bound().get_d();
    if (!std::isfinite(lower))
        return std::nullopt;
    return std::clamp(static_cast<long double>(lower) - *_initial_time, 0.0L, *span);
}

Void HybridEvolutionProgressTracker::add_pending(const HybridEnclosure& enclosure) {
    if (_complete)
        throw std::logic_error("Cannot register work after progress completion");
    if (!_origin_recorded) {
        const double initial = enclosure.time_range().lower_bound().get_d();
        if (std::isfinite(initial))
            _initial_time = initial;
        _origin_recorded = true;
    }
    const auto residual = remaining(enclosure);
    const auto prefix = elapsed(enclosure);
    if (prefix)
        _pending_credit += *prefix;
    if (residual)
        _remaining += *residual;
    else
        ++_unknown_pending;
    ++_pending;
}

Void HybridEvolutionProgressTracker::remove_pending(const HybridEnclosure& enclosure) {
    if (_pending == 0u)
        throw std::logic_error("Cannot retire a checkpoint from an empty progress tracker");
    const auto residual = remaining(enclosure);
    const auto prefix = elapsed(enclosure);
    if (prefix)
        _pending_credit = std::max(0.0L, _pending_credit - *prefix);
    if (residual)
        _remaining = std::max(0.0L, _remaining - *residual);
    else if (_unknown_pending != 0u)
        --_unknown_pending;
    else
        throw std::logic_error("Unregistered unknown-time checkpoint");
    if (--_pending == 0u) {
        _remaining = 0.0L;
        _pending_credit = 0.0L;
    }
}

Void HybridEvolutionProgressTracker::publish(const HybridEvolutionWorkItem& work, const HybridEvolutionStepResult& result) {
    remove_pending(work.enclosure);
    if (std::isfinite(result.continuous_evolution_time) && result.continuous_evolution_time >= 0.0)
        _processed += result.continuous_evolution_time;
    else
        _valid_processed = false;
    for (const auto& successor: result.successors)
        add_pending(successor.work.enclosure);
    if (result.successors.empty()) {
        if (const auto span = duration())
            _terminal_credit += *span;
    }
}

Void HybridEvolutionProgressTracker::finish() {
    if (_pending != 0u)
        throw std::logic_error("Cannot complete progress while checkpoints remain pending");
    _complete = true;
}

HybridEvolutionProgress HybridEvolutionProgressTracker::snapshot() const {
    std::optional<double> remaining_estimate;
    const long double credit = _pending_credit + _terminal_credit;
    if (duration() && _unknown_pending == 0u && _remaining <= std::numeric_limits<double>::max() &&
        credit <= std::numeric_limits<double>::max())
        remaining_estimate = static_cast<double>(_remaining);
    const double processed = _valid_processed ? static_cast<double>(_processed) : std::numeric_limits<double>::quiet_NaN();
    return {processed, remaining_estimate, _pending, _complete, static_cast<double>(credit)};
}

namespace {

/** @brief Per-run flat state; numerical callbacks never see its mutable Orbit. */
class FlatHybridEvolutionBackend final : public HybridEvolutionBackendInterface {
    std::unique_ptr<Orbit<HybridEnclosure>> _orbit;
    HybridEvolutionProgressTracker _progress;
    mutable std::mutex _mutex;
    Bool _started = false;
    Bool _complete = false;

    /** @brief Read published counters under the same lock used for publication. */
    HybridEvolutionStatistics statistics() {
        std::lock_guard<std::mutex> lock(_mutex);
        return {_orbit->reach().size(), _orbit->intermediate().size(), _orbit->final().size()};
    }

    /** @brief Publish contributions before making any successor eligible to run. */
    Void publish(const HybridEvolutionWorkItem& work, const HybridEvolutionStepResult& result) {
        std::lock_guard<std::mutex> lock(_mutex);
        _progress.publish(work, result);
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
    /** @brief Preserve the input and initialise accounting for its remaining evolution. */
    FlatHybridEvolutionBackend(const HybridEnclosure& initial, const Real& final_time)
        : _orbit(std::make_unique<Orbit<HybridEnclosure>>(initial)),
          _progress(final_time.get_d()) { _progress.add_pending(initial); }

    /** @brief Read both queued and executing work under the publication lock. */
    HybridEvolutionProgress progress() const override {
        std::lock_guard<std::mutex> lock(_mutex);
        return _progress.snapshot();
    }

    /** @brief Reuse DynamicWorkload execution while committing complete local results. */
    Void run(const StepFunction& advance) override {
        if (_started) {
            throw std::logic_error("A hybrid evolution backend can only run once");
        }
        _started = true;
        using Workload = DynamicWorkload<HybridEvolutionWorkItem>;
        Workload workload(
            [this](const HybridEvolutionWorkItem&, SharedPointer<ProgressIndicator> indicator) {
                const auto state = progress();
                indicator->update_current(state.complete ? 1.0 : std::min(0.99, state.completion_rate()));
                indicator->update_final(1.0);
            },
            [this, &advance](Workload::Access& pending, const HybridEvolutionWorkItem& item) {
                const auto result = advance(item, statistics());
                publish(item, result);
                // No lock is held while submitting potentially concurrent work.
                for (const auto& successor: result.successors) {
                    pending.append(successor.work);
                }
            });
        workload.append({_orbit->initial(), HybridEvolutionPhase::PROCESS_STARTING_EVENTS}).process();
        std::lock_guard<std::mutex> lock(_mutex);
        _progress.finish();
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
