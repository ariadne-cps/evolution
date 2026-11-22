/***************************************************************************
 *            hybrid_timing.hpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
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
    const HybridEnclosure& initial, const HybridAutomatonInterface& system,
    const ValidatedFunctionPatchFactory& patch_factory,
    const EquationSolverInterface& solver, Bool allow_creep, Bool allow_unwind);
} // namespace EvolutionDetail
} // namespace Ariadne
#endif
