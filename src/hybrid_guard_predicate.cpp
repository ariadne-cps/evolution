/***************************************************************************
 *            src/hybrid_guard_predicate.cpp
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

/** @file src/hybrid_guard_predicate.cpp
 *  @brief Common validated source guard classification and restriction.
 */

#include "hybrid_guard_predicate.hpp"

namespace Ariadne::EvolutionDetail {

Void restrict_source_guard(HybridEnclosure& source, DiscreteEvent event,
                           const EffectiveScalarMultivariateFunction& positive_indicator,
                           EnclosureFunctionArguments arguments) {
    source.new_activation(event, positive_indicator, arguments);
}

HybridSourceGuardStatus classify_source_guard(
    const HybridEnclosure& source,
    const EffectiveScalarMultivariateFunction& positive_indicator,
    EnclosureFunctionArguments arguments) {
    if (definitely(source.is_empty()))
        return HybridSourceGuardStatus::EXCLUDED;
    if (definitely(source.range_of(positive_indicator, arguments).lower_bound() > 0))
        return HybridSourceGuardStatus::FULL_SOURCE;

    HybridEnclosure candidate = source;
    restrict_source_guard(candidate, DiscreteEvent("_guard_test"), positive_indicator, arguments);
    return definitely(candidate.is_empty())
             ? HybridSourceGuardStatus::EXCLUDED
             : HybridSourceGuardStatus::UNRESOLVED;
}

} // namespace Ariadne::EvolutionDetail
