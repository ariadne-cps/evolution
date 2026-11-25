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
//! \details System's event contract forbids free input variables in guards,
//! invariants and reset expressions. Here there are no auxiliary outputs,
//! so every predicate depends only on the pre-state. Progression universally
//! encloses fresh inputs; event resets do not consume a fresh input.
//! Targets may be deterministic or bounded-input inclusion modes without
//! auxiliary outputs. Resets consume no new input; a subsequent progression
//! introduces its own fresh independent input parameters. Call this gate at
//! zero horizon and from direct numerical entry points.
inline Void validate_discrete_inclusion_source(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source) {
    const DiscreteLocation location = source.location();
    if (!automaton.has_inclusion(location))
        throw std::invalid_argument("Expected a discrete inclusion source mode");
    const DifferenceAlgebraicInclusion& inclusion = automaton.inclusion(location);
    if (inclusion.auxiliary_space().dimension() != 0u ||
        source.auxiliary_space().dimension() != 0u)
        throw std::invalid_argument(
            "Validated discrete inclusion mode cannot have auxiliary outputs");
    for (const DiscreteEvent& event: automaton.events(location)) {
        const DiscreteLocation destination = automaton.target(location, event);
        if (!automaton.has_equation(destination) && !automaton.has_inclusion(destination))
            throw std::invalid_argument("Validated discrete inclusion event target has unsupported dynamics");
        if (automaton.has_inclusion(destination) &&
            automaton.inclusion(destination).auxiliary_space().dimension() != 0u)
            throw std::invalid_argument("Validated discrete inclusion target auxiliaries are not supported");
    }
}

} // namespace Ariadne::EvolutionDetail

#endif // ARIADNE_EVOLUTION_HYBRID_DISCRETE_TIME_INCLUSION_CAPABILITY_HPP
