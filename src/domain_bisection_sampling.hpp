/***************************************************************************
 *            src/domain_bisection_sampling.hpp
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

/** @file src/domain_bisection_sampling.hpp
 *  @brief Paving-independent subdivision of box domains for sampling.
 */

#ifndef ARIADNE_EVOLUTION_DOMAIN_BISECTION_SAMPLING_HPP
#define ARIADNE_EVOLUTION_DOMAIN_BISECTION_SAMPLING_HPP

#include <vector>

#include "geometry/box.hpp"

namespace Ariadne::EvolutionDetail {

/** @brief Bisect the domain along each coordinate into the given leaf counts.
 *  @details Each axis is processed in order; leaves are traversed lower-first.
 *  For nonpowers of two, the lower half receives floor(count/2) leaves,
 *  and the upper half receives the rest. This is a sampling heuristic.
 *  @pre counts.size()==domain.dimension() and every count is positive.
 */
std::vector<UpperBoxType> bisect_domain(
    const UpperBoxType& domain, const std::vector<Nat>& counts);

/** @brief Number of intervals per non-degenerate coordinate at the given depth.
 *  Recombine deliberately uses one fewer bisection than Mince: it is a
 *  lower-resolution heuristic, not a reconstruction of paving recombination.
 */
Nat bisection_count(Nat depth, Bool coarse);

} // namespace Ariadne::EvolutionDetail

#endif // ARIADNE_EVOLUTION_DOMAIN_BISECTION_SAMPLING_HPP
