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

#include <cmath>
#include <stdexcept>
#include <vector>

#include "evolution/initial_point_sampling.hpp"
#include "function/function.hpp"
#include "function/functional.hpp"

namespace Ariadne {
namespace {

using Evolver = DiscreteTimeEvolver<ApproximateTag>;
using PointType = Evolver::ApproximatePointType;
using ValueVector = Vector<FloatDPApproximation>;

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

ValueVector combine_arguments(const Point<FloatDPApproximation>& state,
                              const ValueVector& inputs) {
    ValueVector arguments(state.dimension() + inputs.size(), state.zero_element());
    for (SizeType i = 0u; i != state.dimension(); ++i)
        arguments[i] = state[i];
    for (SizeType i = 0u; i != inputs.size(); ++i)
        arguments[state.dimension() + i] = inputs[i];
    return arguments;
}

PointType with_auxiliary(const DiscreteDynamics& system,
                         const EffectiveVectorMultivariateFunction& auxiliary,
                         const Point<FloatDPApproximation>& state,
                         const ValueVector& inputs) {
    const auto states = system.state_space();
    const auto auxiliaries = system.auxiliary_space();
    const auto output_space = system.state_auxiliary_space();
    const ValueVector values = evaluate(auxiliary, combine_arguments(state, inputs));
    Point<FloatDPApproximation> output(output_space.dimension(), state.zero_element());
    for (SizeType i = 0u; i != output_space.dimension(); ++i) {
        const auto variable = output_space[i];
        output[i] = states.contains(variable) ? state[states[variable]]
                                              : values[auxiliaries[variable]];
    }
    return PointType(output_space, output);
}

ValueVector sampled_inputs(const DifferenceAlgebraicInclusion* inclusion,
                           const InputSamplingPolicy* policy,
                           SizeType trajectory, SizeType step, FloatDP time) {
    if (inclusion == nullptr)
        return ValueVector(0u, FloatDPApproximation(0, dp));
    if (policy == nullptr)
        throw std::invalid_argument("Approximate discrete inclusion requires an explicit input policy");

    const InputSamplingContext context{
        trajectory, step,
        FloatDPApproximation(ExactDouble(time.get_d()), dp),
        FloatDPApproximation(1, dp)};
    const ValueVector inputs = policy->sample(inclusion->inputs(), context);
    if (inputs.size() != inclusion->number_of_inputs())
        throw std::invalid_argument("Discrete input policy produced the wrong dimension");
    for (SizeType i = 0u; i != inputs.size(); ++i) {
        const double value = inputs[i].raw().get_d();
        if (!std::isfinite(value) ||
            value < inclusion->inputs()[i].lower_bound().get_d() ||
            value > inclusion->inputs()[i].upper_bound().get_d())
            throw std::invalid_argument("Discrete input sample lies outside the declared input box");
    }
    return inputs;
}

} // namespace

DiscreteTimeEvolver<ApproximateTag>::DiscreteTimeEvolver(
    const DifferenceAlgebraicEquation& system)
    : _input_sampling(nullptr), _system(std::make_shared<SystemType>(system)),
      _configuration(std::make_shared<ConfigurationType>()) {}

DiscreteTimeEvolver<ApproximateTag>::DiscreteTimeEvolver(
    const DifferenceAlgebraicInclusion& system, const InputSamplingPolicy& policy)
    : _input_sampling(policy.clone()), _system(std::make_shared<SystemType>(system)),
      _configuration(std::make_shared<ConfigurationType>()) {}

DiscreteTimeEvolver<ApproximateTag>::DiscreteTimeEvolver(const DiscreteTimeEvolver& other)
    : _input_sampling(other._input_sampling ? other._input_sampling->clone() : nullptr),
      _system(other._system->clone()),
      _configuration(std::make_shared<ConfigurationType>(*other._configuration)) {}

auto DiscreteTimeEvolver<ApproximateTag>::orbit(
    const UpperBoxType& initial_box, const TerminationType& stop) const -> OrbitType {
    if (initial_box.dimension() != system().dimension())
        throw std::invalid_argument("Discrete approximate initial box dimension differs from the state space");
    return orbit(EvolutionDetail::sample_vector_field_initial_points(
                     initial_box, system().state_space(), _configuration->sampling_counts()),
                 stop);
}

auto DiscreteTimeEvolver<ApproximateTag>::orbit(
    const ExactBoxType& initial_box, const TerminationType& stop) const -> OrbitType {
    return orbit(UpperBoxType(initial_box), stop);
}

auto DiscreteTimeEvolver<ApproximateTag>::orbit(
    const InitialType& initial, const TerminationType& stop) const -> OrbitType {
    if (stop < Integer(0))
        throw std::invalid_argument("Discrete approximate horizon must be nonnegative");

    const auto& dynamics = system();
    const auto* equation = dynamics.equation();
    const auto* inclusion = dynamics.inclusion();
    const auto& update = equation ? equation->update_function() : inclusion->update_function();
    const auto& auxiliary = equation ? equation->auxiliary_function() : inclusion->auxiliary_function();
    const auto state_space = dynamics.state_space();

    InitialType canonical;
    std::vector<ValueVector> initial_inputs;
    initial_inputs.reserve(initial.size());
    for (const auto& given: initial) {
        const SizeType trajectory = canonical.size();
        const ValueVector input = sampled_inputs(inclusion, _input_sampling.get(), trajectory, 0u, FloatDP(0, dp));
        canonical.push_back(with_auxiliary(dynamics, auxiliary, canonical_state(given, state_space), input));
        initial_inputs.push_back(input);
    }

    OrbitType result(canonical);
    for (SizeType trajectory = 0u; trajectory != canonical.size(); ++trajectory) {
        auto state = canonical_state(canonical[trajectory], state_space);
        ValueVector inputs = initial_inputs[trajectory];
        FloatDP time(0, dp);
        SizeType tick_index = 0u;
        for (Integer step = 0; step < stop; ++step) {
            state = Point<FloatDPApproximation>(
                evaluate(update, combine_arguments(state, inputs)));
            const FloatDP next_time(step + Integer(1), near, dp);
            if (next_time == time)
                throw std::overflow_error("Discrete approximate tick exceeds floating-point time precision");
            time = next_time;
            ++tick_index;
            inputs = sampled_inputs(inclusion, _input_sampling.get(), trajectory, tick_index, time);
            result.insert(time, with_auxiliary(dynamics, auxiliary, state, inputs), trajectory);
        }
    }
    return result;
}

OutputStream& DiscreteTimeEvolverConfiguration<ApproximateTag>::_write(OutputStream& os) const {
    return os << "DiscreteTimeEvolverConfiguration<ApproximateTag>(sampling_counts="
              << sampling_counts() << ")";
}

} // namespace Ariadne
