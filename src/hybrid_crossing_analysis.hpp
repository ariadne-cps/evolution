/***************************************************************************
 *            hybrid_crossing_analysis.hpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#ifndef ARIADNE_EVOLUTION_HYBRID_CROSSING_ANALYSIS_HPP
#define ARIADNE_EVOLUTION_HYBRID_CROSSING_ANALYSIS_HPP

#include "evolution/hybrid_crossing.hpp"
#include "evolution/hybrid_enclosure.hpp"

namespace Ariadne {
class EquationSolverInterface;
namespace EvolutionDetail {
/** @brief Equation-flow crossing analysis, independent of events and scheduling.
 *  Callers supply a state-only autonomous flow and the configured solver.
 */
Map<DiscreteEvent, CrossingData> analyze_crossings(
    const Set<DiscreteEvent>& active_events,
    const EffectiveVectorMultivariateFunction& dynamic,
    const Map<DiscreteEvent, EffectiveScalarMultivariateFunction>& guards,
    const ValidatedVectorMultivariateFunctionPatch& flow,
    const HybridEnclosure& initial, const EquationSolverInterface& solver);
} // namespace EvolutionDetail
} // namespace Ariadne
#endif
