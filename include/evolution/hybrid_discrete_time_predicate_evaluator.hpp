/***************************************************************************
 *            hybrid_discrete_time_predicate_evaluator.hpp
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

/** @file evolution/hybrid_discrete_time_predicate_evaluator.hpp
 *  @brief Certified three-valued evidence for deterministic hybrid guards.
 */

#ifndef ARIADNE_HYBRID_DISCRETE_TIME_PREDICATE_EVALUATOR_HPP
#define ARIADNE_HYBRID_DISCRETE_TIME_PREDICATE_EVALUATOR_HPP

#include "evolution/hybrid_discrete_time_step_planner.hpp"
#include "evolution/hybrid_enclosure.hpp"

namespace Ariadne {

/** @brief Certify a symbolic predicate over an entire hybrid state enclosure.
 *  @details The source mode must contain a DifferenceAlgebraicEquation.
 *  Explicit algebraic auxiliaries are substituted into the predicate, so
 *  dependent outputs remain correlated with the state. The indicator's
 *  validated range is computed over the entire constrained source enclosure,
 *  rather than over a sampled point or an independent auxiliary box.
 *
 *  The shared source-guard classifier, also used for continuous hybrid
 *  initial activation, certifies TRUE on a nonempty source with a strictly
 *  positive indicator; FALSE if the closed active-side candidate is provably
 *  empty; and UNKNOWN otherwise. This remains conservative for strict guards
 *  because the closed restriction includes the boundary.
 *  No LOWER image guarantee follows from the returned evidence.
 *
 *  Bounded-input inclusion modes support predicates over state and
 *  algebraic auxiliary outputs. The current independent input box is
 *  represented by correlated INPUT coordinates, so auxiliary guards
 *  preserve state/input dependencies in TRUE/FALSE/UNKNOWN evidence.
 */
class HybridDiscreteTimePredicateEvaluator {
  public:
    static DiscreteTimePredicateEvidence classify(
        const DiscreteTimeHybridAutomatonInterface& automaton,
        const HybridEnclosure& source,
        const ContinuousPredicate& predicate);
};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_DISCRETE_TIME_PREDICATE_EVALUATOR_HPP
