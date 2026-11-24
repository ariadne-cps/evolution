/***************************************************************************
 *            hybrid_discrete_time_predicate_evaluator.cpp
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

/** @file src/hybrid_discrete_time_predicate_evaluator.cpp
 *  @brief Interval-certified hybrid predicate evidence on constrained enclosures.
 */

#include "evolution/hybrid_discrete_time_predicate_evaluator.hpp"

#include <stdexcept>

#include "symbolic/assignment.hpp"
#include "system/system_function.hpp"

namespace Ariadne {

DiscreteTimePredicateEvidence HybridDiscreteTimePredicateEvaluator::classify(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source,
    const ContinuousPredicate& predicate) {
    const DiscreteLocation location = source.location();
    if (!automaton.has_mode(location))
        throw std::invalid_argument("Discrete hybrid predicate source mode does not exist");
    if (!automaton.has_equation(location))
        throw std::invalid_argument("Certified discrete inclusion predicate evaluation is not implemented");

    const DifferenceAlgebraicEquation& equation = automaton.equation(location);
    if (source.state_space() != equation.state_space())
        throw std::invalid_argument("Discrete hybrid predicate source state space mismatch");

    // Build a scalar function of the physical state only. Substituting the
    // sorted algebraic assignments prevents treating outputs such as a=x
    // as independent interval variables when classifying a-x+1 >= 0.
    const EffectiveScalarMultivariateFunction indicator_function =
        make_constraint_function(equation.state_space(),
                                 algebraic_sort(equation.auxiliary_assignments()),
                                 predicate, Sign::POSITIVE);

    // Composing with the enclosure's state map preserves correlations in
    // its parameterization and retains the constraints on its source domain.
    const UpperIntervalType range =
        source.range_of(indicator_function, EnclosureFunctionArguments::STATE);

    if (definitely(range.lower_bound() > 0))
        return DiscreteTimePredicateEvidence::TRUE;
    if (definitely(range.upper_bound() < 0))
        return DiscreteTimePredicateEvidence::FALSE;
    return DiscreteTimePredicateEvidence::UNKNOWN;
}

} // namespace Ariadne
