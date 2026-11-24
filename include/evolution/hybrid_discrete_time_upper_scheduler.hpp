/***************************************************************************
 *            hybrid_discrete_time_upper_scheduler.hpp
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

/** @file evolution/hybrid_discrete_time_upper_scheduler.hpp
 *  @brief Finite-horizon UPPER exploration for deterministic hybrid discrete dynamics.
 */

#ifndef ARIADNE_HYBRID_DISCRETE_TIME_UPPER_SCHEDULER_HPP
#define ARIADNE_HYBRID_DISCRETE_TIME_UPPER_SCHEDULER_HPP

#include <cstdint>
#include <optional>
#include <vector>

#include "evolution/hybrid_discrete_time_step_executor.hpp"
#include "evolution/hybrid_orbit.hpp"

namespace Ariadne {

//! \brief One reached enclosure, its integer index and incoming branch.
//! \details A parent refers to an earlier entry in the same exploration.
//! A certain precondition is not a LOWER image or reachability certificate.
struct HybridDiscreteTimeUpperNode {
    Integer index;
    HybridEnclosure enclosure;
    std::optional<SizeType> parent;
    std::optional<HybridDiscreteTimePlannedBranch> incoming_branch;
};

//! \brief Why a leaf of a finite-horizon exploration has no children.
enum class HybridDiscreteTimeUpperTerminalKind : std::uint8_t {
    HORIZON, //!< The requested terminal index was reached.
    BLOCKED  //!< Neither progression nor any transition was possible.
};

//! \brief Reference to a leaf node and its termination reason.
struct HybridDiscreteTimeUpperTerminal {
    SizeType node_index;
    HybridDiscreteTimeUpperTerminalKind kind;
};

//! \brief All explored branches, preserving event labels and parentage.
//! \details The initial enclosure is nodes[0]. Every later node has an
//! incoming branch and an earlier parent; every terminal references a node.
//! These are UPPER enclosures, not a certificate of LOWER reachability.
struct HybridDiscreteTimeUpperExploration {
    Integer initial_index;
    Integer termination_index;
    std::vector<HybridDiscreteTimeUpperNode> nodes;
    std::vector<HybridDiscreteTimeUpperTerminal> terminals;
};

//! \brief A projection of a discrete exploration into the existing hybrid orbit.
//! \details reach contains each non-root reached node; intermediate contains
//! non-horizon reached nodes (including blocked leaves); final contains only
//! leaves at the requested horizon. Blocked leaves are also reported separately
//! so they are never mistaken for horizon results. Parentage and event history
//! remain available in the original exploration, not in Orbit.
struct HybridDiscreteTimeUpperOrbitProjection {
    Orbit<HybridEnclosure> orbit;
    ListSet<HybridEnclosure> blocked;
};

//! \brief Explore all possible deterministic hybrid successors up to an index.
//! \details Each branch consumes one tick. Uncertain guard and invariant
//! evidence always retains every possibly admissible successor. The complete
//! source enclosure is imaged separately for each branch (UPPER semantics).
//! No zero-time transitions, implicit event priorities, or LOWER claims.
//!
//! The initial enclosure's entire certified time range must be the singleton
//! initial_index; a range merely containing that index is insufficient.
//! The same invariant is checked for all generated nodes.
//! Unsupported difference inclusions throw rather
//! than being sampled. A hard maximum_nodes budget aborts with an exception
//! instead of dropping branches or returning an incomplete exploration.
//! This is a single-purpose scheduler, not yet a public Evolver specialization.
class HybridDiscreteTimeUpperScheduler {
  public:
    static HybridDiscreteTimeUpperExploration evolve(
        const DiscreteTimeHybridAutomatonInterface& automaton,
        const HybridEnclosure& initial,
        Integer initial_index,
        Integer termination_index,
        SizeType maximum_nodes);

    //! \brief Project the explored tree without merging independent branches.
    //! \details Only HORIZON terminals contribute to orbit.final().
    //! BLOCKED terminals are returned separately, including a blocked root.
    static HybridDiscreteTimeUpperOrbitProjection project_orbit(
        const HybridDiscreteTimeUpperExploration& exploration);
};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_DISCRETE_TIME_UPPER_SCHEDULER_HPP
