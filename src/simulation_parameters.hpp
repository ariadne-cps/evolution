/***************************************************************************
 *            simulation_parameters.hpp
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

/** @file src/simulation_parameters.hpp
 *  @brief Prepare explicit parameter values for approximate evolution.
 */

#ifndef ARIADNE_EVOLUTION_SIMULATION_PARAMETERS_HPP
#define ARIADNE_EVOLUTION_SIMULATION_PARAMETERS_HPP

#include "evolution/input_sampling.hpp"
#include "function/function.hpp"
#include "function/functional.hpp"
#include "geometry/box.hpp"
#include "symbolic/space.hpp"
#include "symbolic/valuation.hpp"
#include "system/odai.hpp"

namespace Ariadne::EvolutionDetail {
Void validate_simulation_parameters(const RealSpace&, const ExactBoxType&,
                                    const RealValuation&, bool exact_names);

template<class F>
F bind_simulation_parameters(const F& function, SizeType retained_arguments,
                             const RealSpace& space, const RealValuation& values, SizeType input_count = 0u) {
    ARIADNE_PRECONDITION(function.argument_size() == retained_arguments + space.dimension());
    ARIADNE_PRECONDITION(input_count <= retained_arguments);
    if (space.dimension() == 0u)
        return function;
    List<EffectiveScalarMultivariateFunction> arguments;
    for (SizeType i = 0u; i != retained_arguments - input_count; ++i)
        arguments.append(EffectiveScalarMultivariateFunction::coordinate(retained_arguments, i));
    for (const auto& variable: space.variables())
        arguments.append(EffectiveScalarMultivariateFunction::constant(retained_arguments, values.get(variable)));
    for (SizeType i = retained_arguments - input_count; i != retained_arguments; ++i)
        arguments.append(EffectiveScalarMultivariateFunction::coordinate(retained_arguments, i));
    return compose(function, EffectiveVectorMultivariateFunction(arguments));
}
Void validate_simulation_inclusion(const OrdinaryDifferentialAlgebraicInclusion&, const InputSamplingPolicy*);
EffectiveVectorMultivariateFunction prepare_simulation_inclusion(const OrdinaryDifferentialAlgebraicInclusion&,
                                                                 const RealValuation&, bool timed);
EffectiveVectorMultivariateFunction sample_simulation_field(const EffectiveVectorMultivariateFunction&,
                                                            const OrdinaryDifferentialAlgebraicInclusion*, const InputSamplingPolicy*, const InputSamplingContext&);
} // namespace Ariadne::EvolutionDetail
#endif
