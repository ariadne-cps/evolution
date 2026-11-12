/***************************************************************************
 *            inclusion_evolver_state.hpp
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

/** @file src/inclusion_evolver_state.hpp
 *  @brief State and scheduling of differential-inclusion evolution.
 */

#ifndef ARIADNE_EVOLUTION_INCLUSION_EVOLVER_STATE_HPP
#define ARIADNE_EVOLUTION_INCLUSION_EVOLVER_STATE_HPP
#include <algorithm>
#include <limits>
#include <stdexcept>

#include "evolution/differential_inclusion_evolver.hpp"

namespace Ariadne {
namespace EvolutionDetail {
/** @brief Compute the requested evolution operation. */
inline Nat scheduled_step(Nat current, Nat delay) {
    if (delay >= std::numeric_limits<Nat>::digits)
        throw std::overflow_error("Inclusion schedule delay is not representable");
    Nat offset = Nat(1u) << delay;
    if (current > std::numeric_limits<Nat>::max() - offset)
        throw std::overflow_error("Inclusion schedule step overflow");
    return current + offset;
}
} // namespace EvolutionDetail

/** @brief ScheduledApproximator: state and scheduling of differential-inclusion evolution. */
struct ScheduledApproximator {
    Nat step;
    InclusionIntegrator approximator;

    ScheduledApproximator(Nat s, InclusionIntegrator a) : step(s), approximator(a) {}
};

/** @brief Perform the specified evolution operation. */
inline OutputStream& operator<<(OutputStream& os, const ScheduledApproximator& sa) {
    return os << "(" << sa.step << ":" << sa.approximator << ")";
}

/** @brief ScheduledApproximatorComparator: state and scheduling of differential-inclusion evolution. */
struct ScheduledApproximatorComparator {
    /** @brief Perform the specified evolution operation. */
    inline bool operator()(const ScheduledApproximator& sa1, const ScheduledApproximator& sa2) {
        return (sa1.step > sa2.step);
    }
};

/** @brief InclusionEvolverState: state and scheduling of differential-inclusion evolution. */
class InclusionEvolverState {
  private:
    Nat _step;
    List<ScheduledApproximator> _schedule;
    Map<InclusionIntegrator, Nat> _approximator_check_delay;
    Map<InclusionIntegrator, Nat> _approximator_global_optima_count;
    Map<InclusionIntegrator, Nat> _approximator_local_optima_count;

  public:
    InclusionEvolverState(const DifferentialInclusion& ivf, const List<InputApproximation>& approximations, const IntegratorInterface& integrator)
        : _step(0u) {
        /** @brief Perform the specified evolution operation. */
        InclusionIntegratorFactory factory(integrator);
        bool has_reach_rigorous_approximations = false;
        for (auto appro: approximations) {
            if (appro.is_reach_rigorous())
                has_reach_rigorous_approximations = true;
            InclusionIntegrator approximator = factory.create(ivf.function(), ivf.inputs(), appro);
            _schedule.push_back(ScheduledApproximator(0u, approximator));
            _approximator_global_optima_count[approximator] = 0;
            _approximator_local_optima_count[approximator] = 0;
            _approximator_check_delay[approximator] = 0;
        }
        if (not has_reach_rigorous_approximations) {
            InclusionIntegrator approximator = factory.create(ivf.function(), ivf.inputs(), ConstantApproximation());
            _schedule.push_back(ScheduledApproximator(0u, approximator));
            _approximator_global_optima_count[approximator] = 0;
            _approximator_local_optima_count[approximator] = 0;
            _approximator_check_delay[approximator] = 0;
        }
    }

    Nat step() const { return _step; }
    const List<ScheduledApproximator>& schedule() const { return _schedule; }
    const Map<InclusionIntegrator, Nat>& check_delay() const { return _approximator_check_delay; }
    const Map<InclusionIntegrator, Nat>& global_optima_count() const { return _approximator_global_optima_count; }
    const Map<InclusionIntegrator, Nat>& local_optima_count() const { return _approximator_local_optima_count; }

    /** @brief Perform the specified evolution operation. */
    Void reset_local_optima_count() {
        for (auto entry: _approximator_local_optima_count)
            _approximator_local_optima_count[entry.first] = 0;
    }

    /** @brief Perform the specified evolution operation. */
    Void update_with_best(const InclusionIntegrator& best) {
        Nat max_reach_rigorous_index = 0;
        for (auto appro: approximators_to_use()) {
            if (appro.is_reach_rigorous())
                max_reach_rigorous_index = std::max(max_reach_rigorous_index, appro.index());
        }
        for (auto appro: approximators_to_use()) {
            if (best == appro)
                _reset_check_delay(appro);
            else {
                if (not best.is_reach_rigorous() and appro.index() == max_reach_rigorous_index)
                    /** @brief Perform the specified evolution operation. */
                    _reset_check_delay(appro);
                /** @brief Perform the specified evolution operation. */
                else
                    _increase_check_delay(appro);
            }

            _schedule.push_back(ScheduledApproximator(
                EvolutionDetail::scheduled_step(_step, _approximator_check_delay[appro]), appro));
        }
        std::sort(_schedule.begin(), _schedule.end(), ScheduledApproximatorComparator());

        /** @brief Perform the specified evolution operation. */
        _increase_optima_count(best);
    }

    /** @brief Perform the specified evolution operation. */
    List<InclusionIntegrator> approximators_to_use() const {
        List<InclusionIntegrator> result;
        for (auto i = _schedule.size(); i > 0; --i) {
            auto entry = _schedule[i - 1];
            if (entry.step == _step) {
                result.push_back(entry.approximator);
            } else {
                // No past entries survive next_step().
                break;
            }
        }
        return result;
    }

    /** @brief Perform the specified evolution operation. */
    Void next_step() {
        // Called after update_with_best: every current entry has a
        // strictly future successor, so retirement cannot empty the list.
        while (_schedule.back().step == _step)
            _schedule.pop_back();
        _step++;
    }

  private:
    Void _reset_check_delay(const InclusionIntegrator& approximator) { _approximator_check_delay[approximator] = 0; }
    Void _increase_check_delay(const InclusionIntegrator& approximator) { _approximator_check_delay[approximator]++; }
    /** @brief Perform the specified evolution operation. */
    Void _increase_optima_count(const InclusionIntegrator& approximator) {
        _approximator_global_optima_count[approximator]++;
        _approximator_local_optima_count[approximator]++;
    }
};

} // namespace Ariadne
#endif
