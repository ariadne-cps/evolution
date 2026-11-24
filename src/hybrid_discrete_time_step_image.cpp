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
 *  @brief Validated deterministic hybrid images without a new numerical engine.
 */

#include "evolution/hybrid_discrete_time_step_image.hpp"
#include "hybrid_discrete_time_inclusion_capability.hpp"

#include <stdexcept>

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
// one tick. Inclusion inputs are introduced only for a progression branch.
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
    if (!automaton.has_equation(branch.target))
        throw std::invalid_argument("Validated discrete inclusion target outputs are not implemented");

    const auto& destination = automaton.equation(branch.target);
    const auto reset = ordered_reset(automaton.reset_assignments(location, event),
                                     destination.state_space());
    const auto reset_function = make_reset_function(source_space, source_auxiliaries, reset);
    HybridEnclosure next = source;
    next.apply_reset(event, branch.target, destination.state_space(),
                     reset_function, destination.auxiliary_space(),
                     destination.auxiliary_function());
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
        EvolutionDetail::validate_discrete_inclusion_source(automaton, source);
        if (branch.kind == HybridDiscreteTimeBranchKind::TRANSITION)
            return apply_tick_reset(automaton, source, branch, inclusion.state_space(),
                                    inclusion.auxiliary_assignments());
        if (branch.kind != HybridDiscreteTimeBranchKind::PROGRESSION ||
            branch.target != location || branch.event)
            throw std::invalid_argument("Invalid discrete hybrid inclusion progression metadata");

        // Only progression chooses the full independent input box. Guard
        // restrictions were placed on the pre-state by the executor.
        Enclosure evolved = source.continuous_set();
        const SizeType input_offset = evolved.number_of_parameters();
        evolved.adjoin_parameters(inclusion.inputs(), EnclosureVariableKind::INPUT);
        const ExactBoxType domain = evolved.parameter_domain();
        const auto& factory = evolved.configuration().function_factory();
        auto inputs = factory.create_zeros(inclusion.number_of_inputs(), domain);
        for (SizeType j = 0u; j != inclusion.number_of_inputs(); ++j)
            inputs[j] = factory.create_coordinate(domain, input_offset + j);

        // Compose f(x,u) with the current validated state and the freshly
        // appended input coordinates; preserve the previous tick's inputs.
        const auto updated_state = compose(inclusion.update_function(),
                                           join(evolved.state_function(), inputs));
        evolved.apply_parameterised_fixed_evolve_step(updated_state, StepSizeType(1, 0u));
        HybridEnclosure next = source;
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
