/***************************************************************************
 *            hybrid_timing.hpp
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

/** @file src/hybrid_timing.hpp
 *  @brief Compute hybrid evolution and finishing time maps.
 */

#ifndef ARIADNE_EVOLUTION_HYBRID_TIMING_HPP
#define ARIADNE_EVOLUTION_HYBRID_TIMING_HPP

#include "evolution/hybrid_continuous_time_evolver.hpp"

namespace Ariadne {
namespace EvolutionDetail {
/** @brief Select the hybrid horizon, guard creep and clock unwinding.
 *  Only a certified creep adjustment can retire a crossing from the map.
 *  No result publication, reset, or workload scheduling takes place here.
 */
TimingData estimate_hybrid_timing(
    Real final_time, const FlowFunctionModel& flow,
    Map<DiscreteEvent, CrossingData>& crossings,
    const Map<DiscreteEvent, TransitionData>& transitions,
    const HybridEnclosure& initial, const ContinuousTimeHybridAutomatonInterface& system,
    const ValidatedFunctionPatchFactory& patch_factory,
    const EquationSolverInterface& solver, Bool allow_creep, Bool allow_unwind);
} // namespace EvolutionDetail
} // namespace Ariadne
#endif
