/***************************************************************************
 *            enclosure_augmented_paving.hpp
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

// Outer-paving adapters for continuous and hybrid evolution enclosures.
// These definitions belong to evolution, not to the reasoning analyser.
/** @file evolution/enclosure_augmented_paving.hpp
 *  @brief Outer-paving adapters for evolution enclosures.
 */

#ifndef ARIADNE_EVOLUTION_ENCLOSURE_AUGMENTED_PAVING_HPP
#define ARIADNE_EVOLUTION_ENCLOSURE_AUGMENTED_PAVING_HPP

#include "evolution/enclosure.hpp"
#include "evolution/hybrid_enclosure.hpp"
#include "hybrid/hybrid_augmented_paving.hpp"
#include "paving/augmented_paving.hpp"
#include "paving/paver.hpp"

namespace Ariadne {

/** @brief Add a validated outer approximation of a continuous enclosure to an augmented paving.
 *  The caller explicitly supplies the paver.
 */
Void adjoin_outer_approximation_to(const Enclosure& enclosure, AugmentedPaving& paving, Nat fineness, const Paver& paver);

/** @brief Construct an augmented paving carrying the enclosure's auxiliary mapping
 *  and a validated outer approximation of its state set.
 */
AugmentedPaving outer_approximation(const Enclosure& enclosure, const Grid& grid, Nat fineness, const Paver& paver);

/** @brief Add the location-aware outer approximation of a hybrid enclosure
 *  without coupling to an automaton or reasoning-layer objects.
 */
Void adjoin_outer_approximation_to(const HybridEnclosure& enclosure, HybridAugmentedPaving& paving, Nat fineness, const Paver& paver);

} // namespace Ariadne

#endif // ARIADNE_EVOLUTION_ENCLOSURE_AUGMENTED_PAVING_HPP
