/***************************************************************************
 *            hybrid_discrete_time_step_image.cpp
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

/** @file src/hybrid_discrete_time_step_image.cpp
 *  @brief Validated UPPER images of discrete hybrid equations and bounded-input inclusions.
 */

#include "evolution/hybrid_discrete_time_step_image.hpp"

#include <stdexcept>

#include "hybrid_discrete_time_inclusion_capability.hpp"
#include "system/system_function.hpp"

namespace Ariadne {
namespace {

// The System model checks completeness and unique reset targets, but stores
// assignments in user order. Enclosure reset maps must use the target mode's
// canonical state order, irrespective of that user order.
List<PrimedRealAssignment> ordered_reset(const List<PrimedRealAssignment>& reset,
                                         const RealSpace& target_space) {
    List<PrimedRealAssignment> ordered;
    for (const RealVariable& variable: target_space.variables()) {
        Bool found = false;
        for (const PrimedRealAssignment& assignment: reset) {
            if (assignment.lhs.base() == variable) {
                ordered.append(assignment);
                found = true;
                break;
            }
        }
        if (!found)
            throw std::invalid_argument("Discrete hybrid reset does not assign a target state");
    }
    return ordered;
}

// A reset from either dynamics family depends on the pre-state and consumes
// one tick. The reset itself never samples or consumes a destination input;
// an input-dependent auxiliary at the destination instead needs a new full
// input box at the arrival index, independent of the source's input history.
HybridEnclosure apply_tick_reset(const DiscreteTimeHybridAutomatonInterface& automaton,
                                 const HybridEnclosure& source,
                                 const HybridDiscreteTimePlannedBranch& branch,
                                 const RealSpace& source_space,
                                 const List<RealAssignment>& source_auxiliaries) {
    if (branch.kind != HybridDiscreteTimeBranchKind::TRANSITION || !branch.event)
        throw std::invalid_argument("Invalid discrete hybrid transition branch metadata");
    const DiscreteLocation location = source.location();
    const DiscreteEvent event = *branch.event;
    if (!automaton.has_transition(location, event) ||
        automaton.target(location, event) != branch.target)
        throw std::invalid_argument("Discrete hybrid transition branch is not in the system");
    const RealSpace destination_space = automaton.state_space(branch.target);
    const auto reset = ordered_reset(automaton.reset_assignments(location, event),
                                     destination_space);
    const Bool correlated_auxiliary =
        automaton.has_inclusion(location) &&
        automaton.inclusion(location).auxiliary_space().dimension() != 0u;
    // For an inclusion, use the existing correlated map (x,a(x,u[n]))
    // directly, not an a(x,u[n]) substitution into a state-only function.
    // The source's current input is distinct from the arrival input.
    const RealSpace reset_arguments = correlated_auxiliary
                                        ? join(source_space, source.auxiliary_space())
                                        : source_space;
    const auto reset_function = make_reset_function(
        reset_arguments,
        correlated_auxiliary ? List<RealAssignment>{} : source_auxiliaries,
        reset);
    const EnclosureFunctionArguments source_layout = correlated_auxiliary
                                                       ? EnclosureFunctionArguments::STATE_AUXILIARY
                                                       : EnclosureFunctionArguments::STATE;
    HybridEnclosure next = source;
    if (automaton.has_equation(branch.target)) {
        const auto& destination = automaton.equation(branch.target);
        next.apply_reset(event, branch.target, destination_space,
                         reset_function, destination.auxiliary_space(),
                         destination.auxiliary_function(), source_layout);
    } else if (automaton.has_inclusion(branch.target)) {
        const auto& destination = automaton.inclusion(branch.target);
        // Reset is evaluated on the pre-tick source only, with no target
        // input. Drop the source auxiliary labels before setting the target
        // algebraic outputs (which refer to the *arrival* input).
        const EffectiveVectorMultivariateFunction empty_auxiliary =
            make_auxiliary_function(destination_space, List<RealAssignment>{});
        next.apply_reset(event, branch.target, destination_space,
                         reset_function, RealSpace{}, empty_auxiliary, source_layout);
        if (destination.auxiliary_space().dimension() != 0u) {
            // At arrival n+1 a fresh whole-box u[n+1] is needed for
            // a[n+1]=g(x[n+1],u[n+1]). This input did not participate in
            // the reset and cannot change its pre-state correlations.
            next = EvolutionDetail::prepare_discrete_inclusion_source(automaton, next);
        }
    } else {
        throw std::invalid_argument("Discrete hybrid transition has unsupported target dynamics");
    }
    next.set_time_function(next.time_function() + StepSizeType(1, 0u));
    return next;
}

} // namespace

HybridEnclosure HybridDiscreteTimeStepImage::apply(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source,
    const HybridDiscreteTimePlannedBranch& branch) {
    const DiscreteLocation location = source.location();
    if (!automaton.has_mode(location))
        throw std::invalid_argument("Discrete hybrid image source mode does not exist");
    if (source.state_space() != automaton.state_space(location))
        throw std::invalid_argument("Discrete hybrid image source state space mismatch");
    if (automaton.has_inclusion(location)) {
        const DifferenceAlgebraicInclusion& inclusion = automaton.inclusion(location);
        const HybridEnclosure prepared =
            EvolutionDetail::prepare_discrete_inclusion_source(automaton, source);
        if (branch.kind == HybridDiscreteTimeBranchKind::TRANSITION)
            return apply_tick_reset(automaton, prepared, branch, inclusion.state_space(),
                                    inclusion.auxiliary_assignments());
        if (branch.kind != HybridDiscreteTimeBranchKind::PROGRESSION ||
            branch.target != location || branch.event)
            throw std::invalid_argument("Invalid discrete hybrid inclusion progression metadata");

        // For auxiliary-bearing models a[n]=g(x[n],u[n]) is already
        // attached to the full current u[n] box. For models with no
        // auxiliaries, the input box is introduced on progression as before.
        Enclosure evolved = prepared.continuous_set();
        const Bool has_auxiliary = inclusion.auxiliary_space().dimension() != 0u;
        if (!has_auxiliary)
            evolved.adjoin_parameters(inclusion.inputs(), EnclosureVariableKind::INPUT);
        const auto inputs = evolved.input_parameter_function(inclusion.number_of_inputs());
        const auto updated_state = compose(inclusion.update_function(),
                                           join(evolved.state_function(), inputs));
        evolved.apply_parameterised_fixed_evolve_step(updated_state, StepSizeType(1, 0u));
        // Allocate u[n+1] independently, so the output at the next index
        // uses the new current input, not the just-consumed u[n].
        if (has_auxiliary)
            evolved.adjoin_parameters(inclusion.inputs(), EnclosureVariableKind::INPUT);
        HybridEnclosure next = prepared;
        next.replace_continuous_set(evolved);
        return next;
    }
    if (!automaton.has_equation(location))
        throw std::invalid_argument("Unknown discrete hybrid mode dynamics");

    const DifferenceAlgebraicEquation& equation = automaton.equation(location);
    if (source.auxiliary_space().size() != 0u && source.auxiliary_space() != equation.auxiliary_space())
        throw std::invalid_argument("Discrete hybrid image source auxiliary space mismatch");

    // Rebind equation auxiliaries on the source state before either operation.
    // The symbolic equation already substitutes auxiliaries in update and reset
    // functions. No independent auxiliary values are sampled.
    HybridEnclosure next = source;
    next.set_auxiliary(equation.auxiliary_space().variables(), equation.auxiliary_function());

    if (branch.kind == HybridDiscreteTimeBranchKind::PROGRESSION) {
        if (branch.event || branch.target != location)
            throw std::invalid_argument("Invalid discrete hybrid progression branch metadata");
        Enclosure evolved = next.continuous_set();
        evolved.apply_discrete_time_map_step(equation.update_function());
        next.replace_continuous_set(evolved);
        return next;
    }
    return apply_tick_reset(automaton, next, branch, equation.state_space(),
                            equation.auxiliary_assignments());
}

} // namespace Ariadne
