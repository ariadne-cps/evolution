/***************************************************************************
 *            enclosure_augmented_paving.cpp
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

/** @file src/enclosure_augmented_paving.cpp
 *  @brief Implement rendering-independent outer-paving adapters for enclosures.
 */

// Adapt evolved enclosures into drawing-independent augmented pavings.
#include "evolution/enclosure_augmented_paving.hpp"

#include "evolution/hybrid_enclosure.hpp"
#include "function/constraint.hpp"
#include "function/functional.hpp"
#include "hybrid/hybrid_augmented_paving.hpp"
#include "labelled/figure.hpp"
#include "logging/logging.hpp"
#include "paving/binary_tree_paving.hpp"
#include "set/function_set.hpp"
#include "symbolic/space.hpp"

namespace Ariadne {
namespace {
void require_state_auxiliary_paving(const Enclosure& enclosure) {
    if (enclosure.auxiliary_arguments() == EnclosureFunctionArguments::STATE_TIME && enclosure.auxiliary_mapping().result_size() != 0u)
        throw std::invalid_argument("A state-only augmented paving cannot represent time-dependent auxiliary outputs");
}
} // namespace

Void adjoin_outer_approximation_to(const Enclosure& enclosure, AugmentedPaving& paving, Nat fineness, const Paver& paver) {
    require_state_auxiliary_paving(enclosure);
    if (enclosure.auxiliary_mapping().result_size() != 0 &&
        enclosure.auxiliary_mapping().managed_pointer() == paving.auxiliary_mapping().managed_pointer()) {
        ARIADNE_WARN("enclosure=" << enclosure << ", "
                                                  "enclosure.auxiliary_mapping()="
                                  << enclosure.auxiliary_mapping() << ", "
                                                                      "paving.auxiliary_mapping()="
                                  << paving.auxiliary_mapping());
    }
    paver.adjoin_outer_approximation(paving.state_paving(), enclosure.state_set(), fineness);
}

AugmentedPaving outer_approximation(const Enclosure& enclosure, const Grid& grid, Nat fineness, const Paver& paver) {
    require_state_auxiliary_paving(enclosure);
    const auto auxiliary = enclosure.auxiliary_mapping().result_size() == 0u
                             ? EffectiveVectorMultivariateFunction(0u, EuclideanDomain(enclosure.state_dimension()))
                             : enclosure.auxiliary_mapping();
    AugmentedPaving result(grid, auxiliary);
    adjoin_outer_approximation_to(enclosure, result, fineness, paver);
    return result;
}

Void adjoin_outer_approximation_to(const HybridEnclosure& enclosure, HybridAugmentedPaving& paving, Nat fineness, const Paver& paver) {
    enclosure.adjoin_outer_approximation_to(paving.state_paving(), fineness, paver);
}

} // namespace Ariadne
