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
#include "hybrid_discrete_time_predicate_function.hpp"
#include "hybrid_guard_predicate.hpp"

#include <stdexcept>

namespace Ariadne {

DiscreteTimePredicateEvidence HybridDiscreteTimePredicateEvaluator::classify(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source,
    const ContinuousPredicate& predicate) {
    // Public direct classification also needs a full current input box.
    const HybridEnclosure prepared = automaton.has_inclusion(source.location())
        ? EvolutionDetail::prepare_discrete_inclusion_source(automaton, source)
        : source;
    const auto indicator =
        EvolutionDetail::discrete_source_indicator(automaton, prepared, predicate);
    switch (EvolutionDetail::classify_source_guard(
        prepared, indicator.function, indicator.arguments)) {
        case EvolutionDetail::HybridSourceGuardStatus::FULL_SOURCE:
            return DiscreteTimePredicateEvidence::TRUE;
        case EvolutionDetail::HybridSourceGuardStatus::EXCLUDED:
            return DiscreteTimePredicateEvidence::FALSE;
        case EvolutionDetail::HybridSourceGuardStatus::UNRESOLVED:
            return DiscreteTimePredicateEvidence::UNKNOWN;
        default:
            throw std::logic_error("Unhandled hybrid source predicate status");
    }
}

} // namespace Ariadne
