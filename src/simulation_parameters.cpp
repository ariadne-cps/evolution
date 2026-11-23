/***************************************************************************
 *            simulation_parameters.cpp
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

/** @file src/simulation_parameters.cpp
 *  @brief Prepare explicit parameter values for approximate evolution.
 */

#include "simulation_parameters.hpp"

#include <cmath>
#include <stdexcept>

namespace Ariadne::EvolutionDetail {
Void validate_simulation_parameters(const RealSpace& space, const ExactBoxType& domain,
                                    const RealValuation& values, bool exact_names) {
    ARIADNE_PRECONDITION(space.dimension() == domain.dimension());
    if (exact_names && values.values().size() != space.dimension())
        throw std::invalid_argument("Approximate evolution requires exactly the declared constant parameter names");
    for (SizeType i = 0u; i != space.dimension(); ++i) {
        const auto variable = space[i];
        if (!values.values().has_key(variable.name()))
            throw std::invalid_argument("Approximate evolution requires an explicit value for every model parameter");
        const FloatDPBounds bounds(values.get(variable), dp);
        if (!std::isfinite(bounds.lower_raw().get_d()) || !std::isfinite(bounds.upper_raw().get_d()))
            throw std::invalid_argument("Approximate evolution requires finite constant parameter values");
        const ExactIntervalType enclosure(bounds.lower_raw(), bounds.upper_raw());
        if (!definitely(subset(enclosure, domain[i])))
            throw std::invalid_argument("Cannot certify the selected constant parameter value within its declared domain");
    }
}
Void validate_simulation_inclusion(const OrdinaryDifferentialAlgebraicInclusion& inclusion, const InputSamplingPolicy* policy) {
    if (policy == nullptr)
        throw std::invalid_argument("Approximate ODAI evolution requires an explicit input sampling policy");
    if (!inclusion.has_input_function())
        throw std::invalid_argument("Approximate ODAI evolution requires a structured effective input function");
}
EffectiveVectorMultivariateFunction prepare_simulation_inclusion(const OrdinaryDifferentialAlgebraicInclusion& inclusion,
                                                                 const RealValuation& values, bool timed) {
    auto field = inclusion.function();
    const SizeType n = inclusion.dimension(), m = inclusion.number_of_inputs();
    ARIADNE_PRECONDITION(timed || !inclusion.has_explicit_time());
    if (timed && !inclusion.has_explicit_time()) {
        List<EffectiveScalarMultivariateFunction> arguments;
        for (SizeType i = 0u; i != field.argument_size(); ++i)
            arguments.append(EffectiveScalarMultivariateFunction::coordinate(field.argument_size() + 1u, i < n ? i : i + 1u));
        field = compose(field, EffectiveVectorMultivariateFunction(arguments));
    }
    return bind_simulation_parameters(field, n + (timed ? 1u : 0u) + m, inclusion.parameter_space(), values, m);
}
EffectiveVectorMultivariateFunction sample_simulation_field(const EffectiveVectorMultivariateFunction& field,
                                                            const OrdinaryDifferentialAlgebraicInclusion* inclusion, const InputSamplingPolicy* policy, const InputSamplingContext& context) {
    if (inclusion == nullptr)
        return field;
    validate_simulation_inclusion(*inclusion, policy);
    const auto sample = policy->sample(inclusion->inputs(), context);
    const SizeType m = inclusion->number_of_inputs();
    if (sample.size() != m)
        throw std::invalid_argument("Input sampling policy returned the wrong input dimension");
    ARIADNE_PRECONDITION(field.argument_size() >= m);
    const SizeType retained = field.argument_size() - m;
    List<EffectiveScalarMultivariateFunction> arguments;
    for (SizeType i = 0u; i != retained; ++i)
        arguments.append(EffectiveScalarMultivariateFunction::coordinate(retained, i));
    for (SizeType i = 0u; i != m; ++i) {
        const double value = sample[i].raw().get_d();
        if (!std::isfinite(value) || value < inclusion->inputs()[i].lower_bound().get_d() || value > inclusion->inputs()[i].upper_bound().get_d())
            throw std::invalid_argument("Input sampling policy returned a value outside its numerical input box");
        arguments.append(EffectiveScalarMultivariateFunction::constant(retained, Real(ExactDouble(value))));
    }
    return compose(field, EffectiveVectorMultivariateFunction(arguments));
}
} // namespace Ariadne::EvolutionDetail
