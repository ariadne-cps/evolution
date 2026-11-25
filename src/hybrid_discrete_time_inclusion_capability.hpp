/***************************************************************************
 *            hybrid_discrete_time_inclusion_capability.hpp
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

/** @file src/hybrid_discrete_time_inclusion_capability.hpp
 *  @brief Shared supported-mode validation for bounded-input discrete hybrids.
 */

#ifndef ARIADNE_EVOLUTION_HYBRID_DISCRETE_TIME_INCLUSION_CAPABILITY_HPP
#define ARIADNE_EVOLUTION_HYBRID_DISCRETE_TIME_INCLUSION_CAPABILITY_HPP

#include <stdexcept>
#include "evolution/hybrid_enclosure.hpp"
#include "system/discrete_time_hybrid_automaton_interface.hpp"

namespace Ariadne::EvolutionDetail {

//! \brief Check the common capability boundary for a validated inclusion mode.
//! \details An auxiliary-bearing inclusion may have outgoing events only
//! when their guard, invariant and reset expressions depend on state alone.
//! Predicates or resets that reference input-dependent auxiliary outputs
//! require a separate certified correlation-aware implementation.
//! Inclusion targets with auxiliary outputs remain excluded.
//! Apply this gate even at zero horizon and at direct numerical entry points.
inline Void validate_discrete_inclusion_source(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source) {
    const DiscreteLocation location = source.location();
    if (!automaton.has_inclusion(location))
        throw std::invalid_argument("Expected a discrete inclusion source mode");
    const DifferenceAlgebraicInclusion& inclusion = automaton.inclusion(location);
    const Bool has_auxiliary = inclusion.auxiliary_space().dimension() != 0u;
    if (source.auxiliary_space().dimension() != 0u &&
        source.auxiliary_space() != inclusion.auxiliary_space())
        throw std::invalid_argument("Validated discrete inclusion source auxiliary-space mismatch");
    // A guard/invariant/reset depending only on state does not consume or
    // constrain the independent current input, even when the source mode
    // has auxiliary outputs a(x,u). Reject only the event expressions that
    // actually depend on auxiliary variables: substituting those into a
    // state-only indicator or reset map would be unsound.
    const Set<UntypedVariable> state_variables(make_set(inclusion.state_space().variables()));
    for (const DiscreteEvent& event: automaton.events(location)) {
        if (has_auxiliary) {
            if (!subset(automaton.guard_predicate(location, event).arguments(), state_variables) ||
                !subset(automaton.invariant_predicate(location, event).arguments(), state_variables))
                throw std::invalid_argument(
                    "Validated inclusion event guard/invariant with input-dependent auxiliary is unsupported");
            for (const PrimedRealAssignment& assignment: automaton.reset_assignments(location, event))
                if (!subset(assignment.rhs.arguments(), state_variables))
                    throw std::invalid_argument(
                        "Validated inclusion event reset with input-dependent auxiliary is unsupported");
        }
        const DiscreteLocation destination = automaton.target(location, event);
        if (!automaton.has_equation(destination) && !automaton.has_inclusion(destination))
            throw std::invalid_argument("Validated discrete inclusion event target has unsupported dynamics");
        if (automaton.has_inclusion(destination) &&
            automaton.inclusion(destination).auxiliary_space().dimension() != 0u)
            throw std::invalid_argument("Validated discrete inclusion target auxiliaries are not supported");
    }
}

//! \brief Attach the full current input box to an auxiliary-bearing
//! inclusion source, including at zero horizon and before state-only events.
//! Earlier INPUT axes remain independent and correlated with the state.
inline HybridEnclosure prepare_discrete_inclusion_source(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source) {
    validate_discrete_inclusion_source(automaton, source);
    const auto& model = automaton.inclusion(source.location());
    if (model.auxiliary_space().dimension() == 0u)
        return source;
    HybridEnclosure prepared = source;
    if (source.auxiliary_space().dimension() == 0u) {
        Enclosure extended = source.continuous_set();
        extended.adjoin_parameters(model.inputs(), EnclosureVariableKind::INPUT);
        prepared.replace_continuous_set(extended);
    } else {
        const Enclosure& set = source.continuous_set();
        if (set.auxiliary_arguments() != EnclosureFunctionArguments::STATE_INPUTS ||
            set.auxiliary_parameter_count() != model.number_of_inputs())
            throw std::invalid_argument("Hybrid inclusion current-input binding is missing");
        SizeType remaining = model.number_of_inputs();
        const auto domain = set.parameter_domain();
        const auto& kinds = set.variable_kinds();
        for (SizeType i = kinds.size(); i != 0u && remaining != 0u; --i) {
            if (kinds[i - 1u] == EnclosureVariableKind::INPUT) {
                --remaining;
                if (!definitely(subset(domain[i - 1u], model.inputs()[remaining])))
                    throw std::invalid_argument("Hybrid inclusion current input is outside declared box");
            }
        }
        if (remaining != 0u)
            throw std::invalid_argument("Hybrid inclusion current input coordinates are missing");
    }
    prepared.set_auxiliary(model.auxiliary_space().variables(), model.auxiliary_function(),
                           EnclosureFunctionArguments::STATE_INPUTS, model.number_of_inputs());
    return prepared;
}

} // namespace Ariadne::EvolutionDetail

#endif // ARIADNE_EVOLUTION_HYBRID_DISCRETE_TIME_INCLUSION_CAPABILITY_HPP
