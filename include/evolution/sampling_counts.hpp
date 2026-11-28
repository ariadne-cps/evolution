/***************************************************************************
 *            sampling_counts.hpp
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

/** @file evolution/sampling_counts.hpp
 *  @brief Sampling subinterval counts by named coordinate.
 */

#ifndef ARIADNE_EVOLUTION_SAMPLING_COUNTS_HPP
#define ARIADNE_EVOLUTION_SAMPLING_COUNTS_HPP

#include "symbolic/space.hpp"
#include "utility/container.hpp"

namespace Ariadne {

/** @brief Number of sampling subintervals per named variable.
 *  @details Unspecified axes use one subinterval. Counts must be positive.
 *  Counts are numbers of subintervals, not numbers of bisection steps.
 */
using SamplingCounts = Map<RealVariable, Nat>;

} // namespace Ariadne

#endif // ARIADNE_EVOLUTION_SAMPLING_COUNTS_HPP
