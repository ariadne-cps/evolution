/***************************************************************************
 *            hybrid_crossing_analysis.hpp
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

/** @file src/hybrid_crossing_analysis.hpp
 *  @brief Analyze equation crossing geometry and numerical evidence.
 */

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
