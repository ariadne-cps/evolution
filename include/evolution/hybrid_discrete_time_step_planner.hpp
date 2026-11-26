/***************************************************************************
 *            hybrid_discrete_time_step_planner.hpp
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

/** @file evolution/hybrid_discrete_time_step_planner.hpp
 *  @brief Sound predicate-evidence scheduling for one integer-time hybrid step.
 */

#ifndef ARIADNE_HYBRID_DISCRETE_TIME_STEP_PLANNER_HPP
#define ARIADNE_HYBRID_DISCRETE_TIME_STEP_PLANNER_HPP

#include <cstdint>
#include <functional>
#include <optional>
#include <vector>

#include "system/discrete_time_hybrid_automaton_interface.hpp"

namespace Ariadne {

//! \brief Truth of a predicate over the *entire* source set.
//! \details UNKNOWN admits both truth values: this planner never resolves it
//! by sampling or by assuming that a branch is disabled.
enum class DiscreteTimePredicateEvidence : std::uint8_t {
    FALSE,   //!< Certified false throughout the source set.
    UNKNOWN, //!< Neither uniform truth nor uniform falsity certified.
    TRUE     //!< Certified true throughout the source set.
};

//! \brief Whether the predicate is a guard or an explicit advancement invariant.
enum class DiscreteTimePredicateRole : std::uint8_t { GUARD,
                                                      PROGRESSION_INVARIANT };

//! \brief Type of a one-tick branch; each branch advances the index by exactly one.
enum class HybridDiscreteTimeBranchKind : std::uint8_t { PROGRESSION,
                                                         TRANSITION };

//! \brief A possible next-index branch and what has been proved about its precondition.
//! \details certain means that its *predicate precondition* holds over the
//! whole input set; it is not proof of a LOWER numerical image.
struct HybridDiscreteTimePlannedBranch {
    HybridDiscreteTimeBranchKind kind;
    DiscreteLocation target;
    std::optional<DiscreteEvent> event;
    Bool certain;
};

//! \brief A complete one-step plan; all branches are alternatives, not a sequence.
//! \details Source state and all predicates are evaluated at source_index.
//! Each successor, including a reset, belongs to source_index + 1.
struct HybridDiscreteTimeStepPlan {
    Integer source_index;
    Integer successor_index;
    std::vector<HybridDiscreteTimePlannedBranch> branches;
};

//! \brief Schedule all possibly enabled branches without computing numerical images.
//! \details The evaluator must return certified three-valued predicate evidence
//! for the same source enclosure and consistently correlated auxiliary data.
//! For urgent events the progression invariant is inferred from the guard
//! evidence, never evaluated independently. Permissive progression invariants
//! are evaluated explicitly. No event priority or zero-time chaining occurs.
//! The planner does not certify numerical reachability or LOWER semantics.
class HybridDiscreteTimeStepPlanner {
  public:
    using EvidenceEvaluator = std::function<DiscreteTimePredicateEvidence(
        DiscreteLocation, DiscreteEvent, DiscreteTimePredicateRole, const ContinuousPredicate&)>;

    static HybridDiscreteTimeStepPlan plan(const DiscreteTimeHybridAutomatonInterface& automaton,
                                           DiscreteLocation source, Integer source_index,
                                           const EvidenceEvaluator& evaluate);
};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_DISCRETE_TIME_STEP_PLANNER_HPP
