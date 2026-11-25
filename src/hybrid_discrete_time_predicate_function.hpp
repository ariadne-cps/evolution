/***************************************************************************
 *            src/hybrid_discrete_time_predicate_function.hpp
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

/** @file src/hybrid_discrete_time_predicate_function.hpp
 *  @brief Private discrete indicator construction preserving algebraic correlations.
 */

#ifndef ARIADNE_EVOLUTION_HYBRID_DISCRETE_TIME_PREDICATE_FUNCTION_HPP
#define ARIADNE_EVOLUTION_HYBRID_DISCRETE_TIME_PREDICATE_FUNCTION_HPP

#include <stdexcept>
#include "evolution/hybrid_enclosure.hpp"
#include "hybrid_discrete_time_inclusion_capability.hpp"
#include "symbolic/assignment.hpp"
#include "system/discrete_time_hybrid_automaton_interface.hpp"
#include "system/system_function.hpp"

namespace Ariadne::EvolutionDetail {

//! \brief Indicator together with the enclosure argument layout.
struct DiscreteSourceIndicator {
    EffectiveScalarMultivariateFunction function;
    EnclosureFunctionArguments arguments;
};

//! \brief Build a positive-side indicator of a discrete pre-state predicate.
//! \details Inclusion auxiliaries remain correlated with their current inputs.
inline DiscreteSourceIndicator discrete_source_indicator(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source, const ContinuousPredicate& predicate) {
    const DiscreteLocation location = source.location();
    if (!automaton.has_mode(location))
        throw std::invalid_argument("Discrete hybrid predicate source mode does not exist");
    if (source.state_space() != automaton.state_space(location))
        throw std::invalid_argument("Discrete hybrid predicate source state space mismatch");
    if (automaton.has_inclusion(location)) {
        validate_discrete_inclusion_source(automaton, source);
        const auto& inclusion = automaton.inclusion(location);
        if (inclusion.auxiliary_space().dimension() != 0u) {
            // Keep a as an argument rather than substituting a(x,u) into
            // an indicator over x alone. The enclosure supplies the exact
            // correlated (x,a(x,u_current)) validated parameter map.
            return {make_constraint_function(inclusion.state_auxiliary_space(),
                                             List<RealAssignment>{},
                                             predicate, Sign::POSITIVE),
                    EnclosureFunctionArguments::STATE_AUXILIARY};
        }
        return {make_constraint_function(inclusion.state_space(),
                                         List<RealAssignment>{},
                                         predicate, Sign::POSITIVE),
                EnclosureFunctionArguments::STATE};
    }
    if (!automaton.has_equation(location))
        throw std::invalid_argument("Unknown discrete hybrid predicate dynamics");
    const auto& equation = automaton.equation(location);
    return {make_constraint_function(equation.state_space(),
                                     algebraic_sort(equation.auxiliary_assignments()),
                                     predicate, Sign::POSITIVE),
            EnclosureFunctionArguments::STATE};
}

} // namespace Ariadne::EvolutionDetail

#endif // ARIADNE_EVOLUTION_HYBRID_DISCRETE_TIME_PREDICATE_FUNCTION_HPP
