/***************************************************************************
 *            src/hybrid_discrete_time_step_executor.cpp
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

/** @file src/hybrid_discrete_time_step_executor.cpp
 *  @brief One-tick deterministic hybrid UPPER fan-out with event provenance.
 */

#include "evolution/hybrid_discrete_time_step_executor.hpp"
#include "evolution/hybrid_discrete_time_predicate_evaluator.hpp"

#include <stdexcept>

#include "system/system_function.hpp"
#include "symbolic/assignment.hpp"

namespace Ariadne {
namespace {

// The indicator returned by System is positive at points satisfying the
// predicate. Its closed nonnegative half-space is an OUTER approximation,
// including the boundary even if the original comparison is strict.
// Always compose with the same parameterized source state: correlated
// algebraic auxiliaries must never become independent box coordinates.
Void require_prestate_predicate(const DiscreteTimeHybridAutomatonInterface& automaton,
                                HybridEnclosure& source, DiscreteEvent event,
                                const ContinuousPredicate& predicate) {
    const DifferenceAlgebraicEquation& equation = automaton.equation(source.location());
    const EffectiveScalarMultivariateFunction indicator_function =
        make_constraint_function(equation.state_space(),
                                 algebraic_sort(equation.auxiliary_assignments()),
                                 predicate, Sign::POSITIVE);
    const auto parameter_indicator = compose(indicator_function, source.state_function());
    source.new_parameter_constraint(event, parameter_indicator >= ExactNumber(0));
}

// In the certified path only, limit each alternative to the CLOSURE of the
// subset satisfying its necessary pre-state conditions. The supplied callback
// path below keeps its previous contract (no implicit numerical restriction).
HybridEnclosure constrained_branch_source(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source, const HybridDiscreteTimePlannedBranch& branch) {
    HybridEnclosure restricted = source;
    if (branch.kind == HybridDiscreteTimeBranchKind::TRANSITION) {
        if (!branch.event)
            throw std::invalid_argument("Discrete hybrid transition has no event");
        require_prestate_predicate(automaton, restricted, *branch.event,
                                   automaton.guard_predicate(source.location(), *branch.event));
    } else if (branch.kind == HybridDiscreteTimeBranchKind::PROGRESSION) {
        for (const DiscreteEvent& event: automaton.events(source.location())) {
            const auto kind = automaton.transition_kind(source.location(), event);
            if (kind != DiscreteTimeTransitionKind::URGENT &&
                kind != DiscreteTimeTransitionKind::PERMISSIVE)
                throw std::invalid_argument("Unknown discrete hybrid transition kind");
            // For URGENT events the System invariant is the logical
            // complement of the guard; PERMISSIVE uses the stored invariant.
            require_prestate_predicate(automaton, restricted, event,
                                       automaton.invariant_predicate(source.location(), event));
        }
    } else {
        throw std::invalid_argument("Unknown discrete hybrid branch kind");
    }
    return restricted;
}

HybridDiscreteTimeStepExecution execute_impl(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source, Integer source_index,
    const HybridDiscreteTimeStepPlanner::EvidenceEvaluator& evaluate,
    Bool restrict_prestate) {
    if (!automaton.has_mode(source.location()))
        throw std::invalid_argument("Discrete hybrid step source mode does not exist");
    if (!automaton.has_equation(source.location()))
        throw std::invalid_argument("Validated discrete inclusion images are not implemented");

    const HybridDiscreteTimeStepPlan plan =
        HybridDiscreteTimeStepPlanner::plan(automaton, source.location(), source_index, evaluate);

    HybridDiscreteTimeStepExecution result{plan.source_index, plan.successor_index, {}};
    result.successors.reserve(plan.branches.size());
    for (const HybridDiscreteTimePlannedBranch& branch: plan.branches) {
        // Every successor independently images the same original source,
        // restricted only by its own necessary pre-state predicates.
        // A failure aborts the entire operation; there is no partial UPPER.
        const HybridEnclosure branch_source = restrict_prestate
            ? constrained_branch_source(automaton, source, branch) : source;
        result.successors.push_back({
            branch, HybridDiscreteTimeStepImage::apply(automaton, branch_source, branch)});
    }
    return result;
}

} // namespace

HybridDiscreteTimeStepExecution HybridDiscreteTimeStepExecutor::execute_certified(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source, Integer source_index) {
    const auto evaluate = [&](DiscreteLocation, DiscreteEvent,
                              DiscreteTimePredicateRole, const ContinuousPredicate& predicate) {
        return HybridDiscreteTimePredicateEvaluator::classify(automaton, source, predicate);
    };
    return execute_impl(automaton, source, source_index, evaluate, true);
}

HybridDiscreteTimeStepExecution HybridDiscreteTimeStepExecutor::execute(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source, Integer source_index,
    const HybridDiscreteTimeStepPlanner::EvidenceEvaluator& evaluate) {
    return execute_impl(automaton, source, source_index, evaluate, false);
}

} // namespace Ariadne
