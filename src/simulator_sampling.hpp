/***************************************************************************
 *            src/simulator_sampling.hpp
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

/** @file src/simulator_sampling.hpp
 *  @brief Private sampling API for continuous and hybrid simulations.
 */

#ifndef ARIADNE_EVOLUTION_SIMULATOR_SAMPLING_HPP
#define ARIADNE_EVOLUTION_SIMULATOR_SAMPLING_HPP

#include "evolution/hybrid_simulator.hpp"
#include "evolution/vector_field_simulator.hpp"

namespace Ariadne::EvolutionDetail {

/** @brief Sample a continuous initial domain by deterministic bisection.
 * Explicit subspace counts control their respective coordinate; otherwise
 * Mince resolves every non-flat axis at subdivisions, Recombine at one
 * lower resolution, and None takes the domain centre.
 */
VectorFieldSimulator::ApproximateListPointType sample_vector_field_initial_points(
    const UpperBoxType& initial_box, const RealSpace& space,
    const subspace& selected_subspace, DiscretisationType discretisation,
    Nat subdivisions);

/** @brief Sample a hybrid initial domain directly without grid widening.
 * HMince uses fineness per non-degenerate axis; HRecombine uses one
 * fewer bisection. Every resulting point stays in the original box.
 */
HybridSimulator::HybridApproximateListPointType sample_hybrid_initial_points(
    const UpperBoxType& domain, DiscreteLocation location,
    const RealSpace& space, DiscretizationHybridType discretisation,
    Nat fineness);

} // namespace Ariadne::EvolutionDetail

#endif // ARIADNE_EVOLUTION_SIMULATOR_SAMPLING_HPP
