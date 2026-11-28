/***************************************************************************
 *            hybrid_discrete_time_step_planner.cpp
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

/** @file src/hybrid_discrete_time_step_planner.cpp
 *  @brief Structural one-step branching under three-valued predicate evidence.
 */

#include "evolution/hybrid_discrete_time_step_planner.hpp"

#include <stdexcept>

namespace Ariadne {
namespace {

DiscreteTimePredicateEvidence checked_evidence(DiscreteTimePredicateEvidence value) {
    switch (value) {
        case DiscreteTimePredicateEvidence::TRUE:
        case DiscreteTimePredicateEvidence::FALSE:
        case DiscreteTimePredicateEvidence::UNKNOWN:
            return value;
        default:
            throw std::invalid_argument("Unknown predicate evidence");
    }
}

DiscreteTimePredicateEvidence complement(DiscreteTimePredicateEvidence value) {
    switch (value) {
        case DiscreteTimePredicateEvidence::TRUE:
            return DiscreteTimePredicateEvidence::FALSE;
        case DiscreteTimePredicateEvidence::FALSE:
            return DiscreteTimePredicateEvidence::TRUE;
        case DiscreteTimePredicateEvidence::UNKNOWN:
            return DiscreteTimePredicateEvidence::UNKNOWN;
        default:
            throw std::invalid_argument("Unknown predicate evidence");
    }
}

} // namespace

HybridDiscreteTimeStepPlan HybridDiscreteTimeStepPlanner::plan(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    DiscreteLocation source, Integer index,
    const EvidenceEvaluator& evaluate) {
    if (!automaton.has_mode(source))
        throw std::invalid_argument("Discrete hybrid step source mode does not exist");
    if (!evaluate)
        throw std::invalid_argument("Discrete hybrid step requires predicate evidence evaluator");

    HybridDiscreteTimeStepPlan result{index, index + 1, {}};
    Bool progression_possible = true;
    Bool progression_certain = true;

    for (const DiscreteEvent& event: automaton.events(source)) {
        const auto guard = checked_evidence(evaluate(source, event, DiscreteTimePredicateRole::GUARD,
                                                     automaton.guard_predicate(source, event)));
        DiscreteTimePredicateEvidence invariant;
        if (automaton.transition_kind(source, event) == DiscreteTimeTransitionKind::URGENT) {
            // Derive from guard classification: the two predicates are
            // mathematical complements, not independently sampled facts.
            invariant = complement(guard);
        } else if (automaton.transition_kind(source, event) == DiscreteTimeTransitionKind::PERMISSIVE) {
            invariant = checked_evidence(evaluate(source, event, DiscreteTimePredicateRole::PROGRESSION_INVARIANT,
                                                  automaton.invariant_predicate(source, event)));
        } else {
            throw std::invalid_argument("Unknown discrete hybrid transition kind");
        }

        progression_possible = progression_possible && invariant != DiscreteTimePredicateEvidence::FALSE;
        progression_certain = progression_certain && invariant == DiscreteTimePredicateEvidence::TRUE;
        if (guard != DiscreteTimePredicateEvidence::FALSE)
            result.branches.push_back({HybridDiscreteTimeBranchKind::TRANSITION,
                                       automaton.target(source, event), event,
                                       guard == DiscreteTimePredicateEvidence::TRUE});
    }

    if (progression_possible)
        result.branches.push_back({HybridDiscreteTimeBranchKind::PROGRESSION, source,
                                   std::nullopt, progression_certain});
    return result;
}

} // namespace Ariadne
