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
//! \details Input-dependent auxiliary outputs are supported in event-free
//! inclusion source modes only. Guard, invariant and reset evaluation using
//! input-dependent auxiliary outputs requires a separate certified contract.
//! Inclusion targets with auxiliary outputs remain excluded. This gate is
//! applied even at zero horizon.
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
    if (has_auxiliary && !automaton.events(location).empty())
        throw std::invalid_argument(
            "Validated discrete inclusion input-dependent auxiliaries with events are not yet supported");
    for (const DiscreteEvent& event: automaton.events(location)) {
        const DiscreteLocation destination = automaton.target(location, event);
        if (!automaton.has_equation(destination) && !automaton.has_inclusion(destination))
            throw std::invalid_argument("Validated discrete inclusion event target has unsupported dynamics");
        if (automaton.has_inclusion(destination) &&
            automaton.inclusion(destination).auxiliary_space().dimension() != 0u)
            throw std::invalid_argument("Validated discrete inclusion target auxiliaries are not supported");
    }
}

//! \brief Attach the full current input box to an event-free inclusion
//! source whose output depends on current input, including at zero horizon.
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
