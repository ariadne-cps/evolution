/***************************************************************************
 *            discrete_time_approximate_evolver.cpp
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

/** @file src/discrete_time_approximate_evolver.cpp
 *  @brief Approximate discrete map iteration with paving-independent initial sampling.
 */

#include "evolution/discrete_time_approximate_evolver.hpp"

#include <stdexcept>

#include "evolution/initial_point_sampling.hpp"
#include "function/functional.hpp"
#include "function/function.hpp"

namespace Ariadne {
namespace {

using PointType = DiscreteTimeEvolver<ApproximateTag>::ApproximatePointType;

Point<FloatDPApproximation> canonical_state(const PointType& given,
                                            const RealSpace& state_space) {
    Point<FloatDPApproximation> state(state_space.dimension(), given.zero_element());
    for (SizeType i = 0u; i != state_space.dimension(); ++i) {
        if (!given.space().contains(state_space[i]))
            throw std::invalid_argument("Discrete approximate initial point is missing a state variable");
        state[i] = given[given.space()[state_space[i]]];
    }
    return state;
}

PointType with_auxiliary(const DifferenceAlgebraicEquation& system,
                         const Point<FloatDPApproximation>& state) {
    const auto states = system.state_space();
    const auto auxiliaries = system.auxiliary_space();
    const auto output_space = system.state_auxiliary_space();
    const auto values = evaluate(system.auxiliary_function(), Vector<FloatDPApproximation>(state));
    Point<FloatDPApproximation> output(output_space.dimension(), state.zero_element());
    for (SizeType i = 0u; i != output_space.dimension(); ++i) {
        const auto variable = output_space[i];
        output[i] = states.contains(variable) ? state[states[variable]]
                                             : values[auxiliaries[variable]];
    }
    return PointType(output_space, output);
}

} // namespace

DiscreteTimeEvolver<ApproximateTag>::DiscreteTimeEvolver(const SystemType& system)
    : _system(system.clone()), _configuration(std::make_shared<ConfigurationType>()) {}

DiscreteTimeEvolver<ApproximateTag>::DiscreteTimeEvolver(const DiscreteTimeEvolver& other)
    : _system(other._system->clone()),
      _configuration(std::make_shared<ConfigurationType>(*other._configuration)) {}

auto DiscreteTimeEvolver<ApproximateTag>::orbit(
    const UpperBoxType& initial_box, const TerminationType& stop) const -> OrbitType {
    if (initial_box.dimension() != system().dimension())
        throw std::invalid_argument("Discrete approximate initial box dimension differs from the map state space");
    return orbit(EvolutionDetail::sample_vector_field_initial_points(
        initial_box, system().state_space(), _configuration->sampling_counts()), stop);
}

auto DiscreteTimeEvolver<ApproximateTag>::orbit(
    const ExactBoxType& initial_box, const TerminationType& stop) const -> OrbitType {
    return orbit(UpperBoxType(initial_box), stop);
}

auto DiscreteTimeEvolver<ApproximateTag>::orbit(
    const InitialType& initial, const TerminationType& stop) const -> OrbitType {
    if (stop < Integer(0))
        throw std::invalid_argument("Discrete approximate horizon must be nonnegative");

    const auto state_space = system().state_space();
    InitialType canonical;
    for (const auto& given: initial)
        canonical.push_back(with_auxiliary(system(), canonical_state(given, state_space)));

    OrbitType result(canonical);
    for (SizeType trajectory = 0u; trajectory != canonical.size(); ++trajectory) {
        auto state = canonical_state(canonical[trajectory], state_space);
        FloatDP tick(0, dp);
        for (Integer step = 0; step < stop; ++step) {
            state = Point<FloatDPApproximation>(
                evaluate(system().update_function(), Vector<FloatDPApproximation>(state)));
            // Compute the tick from the exact integer index: raw FloatDP
            // addition yields Bounds<FloatDP>, not a FloatDP.
            const FloatDP next_tick(step + Integer(1), near, dp);
            if (next_tick == tick)
                throw std::overflow_error("Discrete approximate tick exceeds floating-point time precision");
            tick = next_tick;
            result.insert(tick, with_auxiliary(system(), state), trajectory);
        }
    }
    return result;
}

OutputStream& DiscreteTimeEvolverConfiguration<ApproximateTag>::_write(OutputStream& os) const {
    return os << "DiscreteTimeEvolverConfiguration<ApproximateTag>(sampling_counts="
              << sampling_counts() << ")";
}

} // namespace Ariadne
