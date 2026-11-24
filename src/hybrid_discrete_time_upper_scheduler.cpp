/***************************************************************************
 *            hybrid_discrete_time_upper_scheduler.cpp
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

/** @file src/hybrid_discrete_time_upper_scheduler.cpp
 *  @brief Breadth-first UPPER fan-out for integer-indexed hybrid systems.
 */

#include "evolution/hybrid_discrete_time_upper_scheduler.hpp"

#include <stdexcept>

namespace Ariadne {

HybridDiscreteTimeUpperExploration HybridDiscreteTimeUpperScheduler::evolve(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& initial,
    Integer initial_index,
    Integer termination_index,
    SizeType maximum_nodes) {
    if (termination_index < initial_index)
        throw std::invalid_argument("Discrete hybrid termination index precedes initial index");
    if (maximum_nodes == 0u)
        throw std::invalid_argument("Discrete hybrid maximum node budget must be positive");
    if (!automaton.has_mode(initial.location()))
        throw std::invalid_argument("Discrete hybrid initial mode does not exist");
    if (!automaton.has_equation(initial.location()))
        throw std::invalid_argument("Discrete hybrid UPPER scheduler does not support inclusion modes");
    if (initial.state_space() != automaton.state_space(initial.location()))
        throw std::invalid_argument("Discrete hybrid initial state space mismatch");

    HybridDiscreteTimeUpperExploration result{initial_index, termination_index, {}, {}};
    result.nodes.push_back({initial_index, initial, std::nullopt, std::nullopt});

    // Parent references are stable indices, not pointers or references
    // that vector growth could invalidate.
    for (SizeType position = 0u; position < result.nodes.size(); ++position) {
        const Integer current_index = result.nodes[position].index;
        if (current_index == termination_index) {
            result.terminals.push_back({position, HybridDiscreteTimeUpperTerminalKind::HORIZON});
            continue;
        }

        // Do not retain a reference to result.nodes[position]: push_back
        // may move the underlying vector. Every branch is evaluated against
        // an independent copy of the same original source enclosure.
        const HybridEnclosure source = result.nodes[position].enclosure;
        const HybridDiscreteTimeStepExecution step =
            HybridDiscreteTimeStepExecutor::execute_certified(
                automaton, source, current_index);

        if (step.source_index != current_index || step.successor_index != current_index + 1)
            throw std::logic_error("Discrete hybrid step produced inconsistent sample indices");
        if (step.successors.size() > maximum_nodes - result.nodes.size())
            throw std::length_error("Discrete hybrid UPPER exploration exceeds node budget");

        if (step.successors.empty()) {
            result.terminals.push_back({position, HybridDiscreteTimeUpperTerminalKind::BLOCKED});
            continue;
        }

        for (const HybridDiscreteTimeStepSuccessor& successor: step.successors)
            result.nodes.push_back({step.successor_index, successor.enclosure,
                                    position, successor.branch});
    }
    return result;
}

} // namespace Ariadne
