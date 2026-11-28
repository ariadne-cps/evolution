/***************************************************************************
 *            src/hybrid_guard_predicate.hpp
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

/** @file src/hybrid_guard_predicate.hpp
 *  @brief Shared validated hybrid source guard classification and restriction.
 */

#ifndef ARIADNE_EVOLUTION_HYBRID_GUARD_PREDICATE_HPP
#define ARIADNE_EVOLUTION_HYBRID_GUARD_PREDICATE_HPP

#include <cstdint>

#include "evolution/hybrid_enclosure.hpp"

namespace Ariadne::EvolutionDetail {

//! \brief Validated status for one entire source, not an inner reach certificate.
enum class HybridSourceGuardStatus : std::uint8_t {
    EXCLUDED,   //!< Closed active-side candidate is provably empty.
    UNRESOLVED, //!< Partial activation, boundary or insufficient evidence.
    FULL_SOURCE //!< Nonempty source has a strictly positive guard range.
};

//! \brief Restrict a source to the closed positive side of a guard indicator.
//! \details The closed cut is an UPPER relaxation of strict predicates.
Void restrict_source_guard(HybridEnclosure& source, DiscreteEvent event,
                           const EffectiveScalarMultivariateFunction& positive_indicator,
                           EnclosureFunctionArguments arguments = EnclosureFunctionArguments::STATE);

//! \brief Classify entire-source activation through validated range and emptiness.
//! \details FULL_SOURCE is never inferred from a restricted candidate.
//! EXCLUDED needs an empty closed candidate; otherwise return UNRESOLVED.
HybridSourceGuardStatus classify_source_guard(
    const HybridEnclosure& source,
    const EffectiveScalarMultivariateFunction& positive_indicator,
    EnclosureFunctionArguments arguments = EnclosureFunctionArguments::STATE);

} // namespace Ariadne::EvolutionDetail

#endif // ARIADNE_EVOLUTION_HYBRID_GUARD_PREDICATE_HPP
