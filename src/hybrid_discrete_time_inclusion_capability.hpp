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
//! \details An auxiliary-bearing inclusion may have outgoing events referring to
//! state or algebraic outputs. They are evaluated against the correlated
//! pre-state/auxiliary enclosure with current input coordinates.
//! Inclusion targets with auxiliary outputs are prepared at reset arrival
//! with an independent full target input box, after evaluating the reset.
//! Apply this gate even at zero horizon and at direct numerical entry points.
inline Void validate_discrete_inclusion_source(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source) {
    const DiscreteLocation location = source.location();
    if (!automaton.has_inclusion(location))
        throw std::invalid_argument("Expected a discrete inclusion source mode");
    const DifferenceAlgebraicInclusion& inclusion = automaton.inclusion(location);
    if (source.auxiliary_space().dimension() != 0u &&
        source.auxiliary_space() != inclusion.auxiliary_space())
        throw std::invalid_argument("Validated discrete inclusion source auxiliary-space mismatch");
    // System already checks every event expression against the state and
    // auxiliary symbols. Both are represented by the correlated enclosure.
    for (const DiscreteEvent& event: automaton.events(location)) {
        const DiscreteLocation destination = automaton.target(location, event);
        if (!automaton.has_equation(destination) && !automaton.has_inclusion(destination))
            throw std::invalid_argument("Validated discrete inclusion event target has unsupported dynamics");
        // A target inclusion with auxiliary outputs is admitted: its
        // fresh arrival input is bound *after* the source reset image.
    }
}

//! \brief Attach the full current input box to an auxiliary-bearing
//! inclusion source, including at zero horizon and before auxiliary guards.
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
