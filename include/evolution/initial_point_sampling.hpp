/***************************************************************************
 *            evolution/initial_point_sampling.hpp
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

/** @file evolution/initial_point_sampling.hpp
 *  @brief Private sampling API for continuous and hybrid simulations.
 */

#ifndef ARIADNE_EVOLUTION_INITIAL_POINT_SAMPLING_HPP
#define ARIADNE_EVOLUTION_INITIAL_POINT_SAMPLING_HPP

#include "evolution/hybrid_continuous_time_approximate_evolver.hpp"
#include "evolution/continuous_time_approximate_evolver.hpp"

namespace Ariadne::EvolutionDetail {

/** @brief Sample at centres of axis-wise subdivided initial boxes.
 * Unspecified or flat axes contribute one subinterval.
 */
ContinuousTimeEvolver<ApproximateTag>::ApproximateListPointType sample_vector_field_initial_points(
    const UpperBoxType& domain, const RealSpace& space,
    const SamplingCounts& counts);

/** @brief Hybrid midpoint sampling with the same per-axis counts. */
HybridContinuousTimeEvolver<ApproximateTag>::HybridApproximateListPointType sample_hybrid_initial_points(
    const UpperBoxType& domain, DiscreteLocation location,
    const RealSpace& space, const SamplingCounts& counts);

} // namespace Ariadne::EvolutionDetail

#endif // ARIADNE_EVOLUTION_INITIAL_POINT_SAMPLING_HPP
