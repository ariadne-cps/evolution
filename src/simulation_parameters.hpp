/***************************************************************************
 *            simulation_parameters.hpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#ifndef ARIADNE_EVOLUTION_SIMULATION_PARAMETERS_HPP
#define ARIADNE_EVOLUTION_SIMULATION_PARAMETERS_HPP

#include "function/function.hpp"
#include "function/functional.hpp"
#include "geometry/box.hpp"
#include "symbolic/space.hpp"
#include "symbolic/valuation.hpp"

namespace Ariadne::EvolutionDetail {
Void validate_simulation_parameters(const RealSpace&, const ExactBoxType&,
                                    const RealValuation&, bool exact_names);

template<class F>
F bind_simulation_parameters(const F& function, SizeType retained_arguments,
                             const RealSpace& space, const RealValuation& values) {
    ARIADNE_PRECONDITION(function.argument_size() == retained_arguments + space.dimension());
    if (space.dimension() == 0u)
        return function;
    List<EffectiveScalarMultivariateFunction> arguments;
    for (SizeType i = 0u; i != retained_arguments; ++i)
        arguments.append(EffectiveScalarMultivariateFunction::coordinate(retained_arguments, i));
    for (const auto& variable: space.variables())
        arguments.append(EffectiveScalarMultivariateFunction::constant(retained_arguments, values.get(variable)));
    return compose(function, EffectiveVectorMultivariateFunction(arguments));
}
} // namespace Ariadne::EvolutionDetail
#endif
