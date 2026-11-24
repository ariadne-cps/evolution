/***************************************************************************
 *            hybrid_discrete_time_step_image.hpp
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

/** @file evolution/hybrid_discrete_time_step_image.hpp
 *  @brief Reuse validated enclosure operations for deterministic hybrid ticks.
 */

#ifndef ARIADNE_HYBRID_DISCRETE_TIME_STEP_IMAGE_HPP
#define ARIADNE_HYBRID_DISCRETE_TIME_STEP_IMAGE_HPP

#include "evolution/hybrid_discrete_time_step_planner.hpp"
#include "evolution/hybrid_enclosure.hpp"

namespace Ariadne {

/** @brief Validated outer image of a structurally selected, deterministic hybrid tick.
 *  @details This adapter reuses Enclosure::apply_discrete_time_map_step
 *  for progression, HybridEnclosure::apply_reset for jumps and the shared
 *  System symbolic reset builder. Predicates are classified separately by
 *  HybridDiscreteTimeStepPlanner; applying an image to the full source enclosure
 *  is an UPPER over-approximation of a possibly enabled branch, not LOWER proof.
 *
 *  The current kernel requires deterministic source and destination modes.
 *  It rejects bounded-input inclusion modes rather than selecting/fixing
 *  inputs or falsely claiming a validated set-valued image.
 *  It performs no scheduling and does not assemble an Orbit.
 */
class HybridDiscreteTimeStepImage {
  public:
    /** @brief Apply one planner branch to an enclosure at its source location.
     *  @throws std::invalid_argument for incompatible source state, invalid
     *  branch metadata, unknown modes or unsupported inclusion dynamics.
     */
    static HybridEnclosure apply(const DiscreteTimeHybridAutomatonInterface& automaton,
                                 const HybridEnclosure& source,
                                 const HybridDiscreteTimePlannedBranch& branch);
};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_DISCRETE_TIME_STEP_IMAGE_HPP
