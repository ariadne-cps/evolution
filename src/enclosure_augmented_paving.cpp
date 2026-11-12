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
#include "function/functional.hpp"
#include "function/constraint.hpp"
#include "set/function_set.hpp"
#include "evolution/enclosure_augmented_paving.hpp"
#include "hybrid/hybrid_augmented_paving.hpp"
#include "evolution/hybrid_enclosure.hpp"
#include "paving/binary_tree_paving.hpp"
#include "symbolic/space.hpp"
#include "labelled/figure.hpp"
#include "logging/logging.hpp"

namespace Ariadne {

Void adjoin_outer_approximation_to(Enclosure const& enclosure, AugmentedPaving& paving, Nat fineness)
{
    if (enclosure.auxiliary_mapping().result_size()!=0 &&
            enclosure.auxiliary_mapping().managed_pointer()==paving.auxiliary_mapping().managed_pointer())
    {
        ARIADNE_WARN("enclosure="<<enclosure<<", "
                    "enclosure.auxiliary_mapping()="<<enclosure.auxiliary_mapping()<<", "
                    "paving.auxiliary_mapping()="<<paving.auxiliary_mapping());
    }
    enclosure.configuration().paver().adjoin_outer_approximation(paving.state_paving(),enclosure.state_set(),fineness);
}


AugmentedPaving outer_approximation(Enclosure const& enclosure, Grid const& grid, Nat fineness)
{
    AugmentedPaving result(grid,enclosure.auxiliary_mapping());
    adjoin_outer_approximation_to(enclosure,result,fineness);
    return result;
}

Void adjoin_outer_approximation_to(HybridEnclosure const& enclosure, HybridAugmentedPaving& paving, Nat fineness)
{
    enclosure.adjoin_outer_approximation_to(paving.state_paving(),fineness);
}

} // namespace Ariadne
