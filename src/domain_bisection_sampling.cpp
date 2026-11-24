/***************************************************************************
 *            src/domain_bisection_sampling.cpp
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

/** @file src/domain_bisection_sampling.cpp
 *  @brief Paving-independent recursive domain bisection.
 */

#include "domain_bisection_sampling.hpp"

#include <limits>
#include <utility>

#include "utility/macros.hpp"

namespace Ariadne {
namespace {

Void bisect_axis(const UpperBoxType& box, SizeType axis, Nat leaves,
                 std::vector<UpperBoxType>& result) {
    if (leaves == 1u) {
        result.push_back(box);
        return;
    }
    const auto halves = box.split(axis);
    const Nat lower = leaves / 2u;
    bisect_axis(halves.first, axis, lower, result);
    bisect_axis(halves.second, axis, leaves - lower, result);
}

} // namespace

namespace EvolutionDetail {

Nat bisection_count(Nat depth, Bool coarse) {
    if (coarse && depth != 0u)
        --depth;
    Nat count = 1u;
    for (Nat i = 0u; i < depth; ++i) {
        ARIADNE_PRECONDITION(count <= std::numeric_limits<Nat>::max() / 2u);
        count *= 2u;
    }
    return count;
}

std::vector<UpperBoxType> bisect_domain(
    const UpperBoxType& domain, const std::vector<Nat>& counts) {
    ARIADNE_PRECONDITION(counts.size() == domain.dimension());
    std::vector<UpperBoxType> result{domain};
    for (SizeType axis = 0u; axis < domain.dimension(); ++axis) {
        const Nat leaves = counts[axis];
        ARIADNE_PRECONDITION(leaves > 0u);
        ARIADNE_PRECONDITION(
            result.size() <= std::numeric_limits<std::size_t>::max() / leaves);
        std::vector<UpperBoxType> next;
        next.reserve(result.size() * leaves);
        for (const auto& box: result)
            bisect_axis(box, axis, leaves, next);
        result = std::move(next);
    }
    return result;
}

} // namespace EvolutionDetail
} // namespace Ariadne
