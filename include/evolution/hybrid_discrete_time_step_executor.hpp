/***************************************************************************
 *            include/evolution/hybrid_discrete_time_step_executor.hpp
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

/** @file include/evolution/hybrid_discrete_time_step_executor.hpp
 *  @brief Execute all possible deterministic one-tick hybrid UPPER branches.
 */

#ifndef ARIADNE_HYBRID_DISCRETE_TIME_STEP_EXECUTOR_HPP
#define ARIADNE_HYBRID_DISCRETE_TIME_STEP_EXECUTOR_HPP

#include <vector>

#include "evolution/hybrid_discrete_time_step_image.hpp"
#include "evolution/hybrid_discrete_time_step_planner.hpp"

namespace Ariadne {

//! \brief One possible successor, paired with its event and predicate evidence.
//! \details A certain precondition is not a LOWER reachability certificate.
struct HybridDiscreteTimeStepSuccessor {
    HybridDiscreteTimePlannedBranch branch;
    HybridEnclosure enclosure;
};

//! \brief All possibly admissible successors of a single integer-time step.
//! \details Each successor is an alternative at successor_index, never a
//! sequence of zero-time events. Empty successors mean no branch was admitted.
struct HybridDiscreteTimeStepExecution {
    Integer source_index;
    Integer successor_index;
    std::vector<HybridDiscreteTimeStepSuccessor> successors;
};

//! \brief Compose the existing hybrid tick planner with validated images.
//! \details The evidence callback must certify predicates over the same
//! complete source enclosure and its correlated auxiliaries. UNKNOWN admits
//! both alternatives. execute_certified restricts each branch to a
//! closed, validated over-approximation of the pre-state satisfying the
//! corresponding guard or all progression invariants. This may include
//! predicate boundaries and does not establish LOWER reachability.
//! execute() with an external evidence callback retains its original
//! full-source imaging contract: the callback alone cannot authorize a
//! numerical restriction. LOWER evolution is not supported.
//! Deterministic equations and event-free bounded-input inclusion
//! progression without auxiliary outputs are supported. Other inclusion
//! operations throw instead of sampling or returning partial results.
//! The caller must
//! maintain consistency between source_index and the enclosure's global
//! time function. This is a single-step adapter, not an orbit scheduler.
class HybridDiscreteTimeStepExecutor {
  public:
    /** @brief Run a deterministic UPPER tick using certified enclosure evidence. */
    static HybridDiscreteTimeStepExecution execute_certified(
        const DiscreteTimeHybridAutomatonInterface& automaton,
        const HybridEnclosure& source, Integer source_index);

    static HybridDiscreteTimeStepExecution execute(
        const DiscreteTimeHybridAutomatonInterface& automaton,
        const HybridEnclosure& source, Integer source_index,
        const HybridDiscreteTimeStepPlanner::EvidenceEvaluator& evaluate);
};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_DISCRETE_TIME_STEP_EXECUTOR_HPP
