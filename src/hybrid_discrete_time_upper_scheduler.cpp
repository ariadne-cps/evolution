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

#include "hybrid_discrete_time_inclusion_capability.hpp"

namespace Ariadne {
namespace {

// The discrete index must describe the whole validated time set, not merely
// intersect it. In particular, ambiguous or non-singleton clocks are rejected.
Void check_time_index(const HybridEnclosure& enclosure, const Integer& index) {
    const Real expected_time(index);
    const UpperIntervalType range = enclosure.time_range();
    if (!definitely(range.lower_bound() >= expected_time) ||
        !definitely(range.upper_bound() <= expected_time))
        throw std::invalid_argument("Discrete hybrid enclosure global time does not match sample index");
}

} // namespace

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
    if (initial.state_space() != automaton.state_space(initial.location()))
        throw std::invalid_argument("Discrete hybrid initial state space mismatch");
    HybridEnclosure prepared = initial;
    if (automaton.has_inclusion(initial.location()))
        prepared = EvolutionDetail::prepare_discrete_inclusion_source(automaton, initial);
    else if (!automaton.has_equation(initial.location()))
        throw std::invalid_argument("Discrete hybrid initial mode has unsupported dynamics");
    check_time_index(prepared, initial_index);

    HybridDiscreteTimeUpperExploration result{initial_index, termination_index, {}, {}};
    result.nodes.push_back({initial_index, prepared, std::nullopt, std::nullopt});

    // Parent references are stable indices, not pointers or references
    // that vector growth could invalidate.
    for (SizeType position = 0u; position < result.nodes.size(); ++position) {
        const Integer current_index = result.nodes[position].index;
        check_time_index(result.nodes[position].enclosure, current_index);
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

        for (const HybridDiscreteTimeStepSuccessor& successor: step.successors) {
            check_time_index(successor.enclosure, step.successor_index);
            result.nodes.push_back({step.successor_index, successor.enclosure,
                                    position, successor.branch});
        }
    }
    return result;
}

// Orbit is a set-oriented view. The exploration retains the complete tree,
// including parent indices and event evidence, without merging alternatives.
HybridDiscreteTimeUpperOrbitProjection HybridDiscreteTimeUpperScheduler::project_orbit(
    const HybridDiscreteTimeUpperExploration& exploration) {
    if (exploration.nodes.empty())
        throw std::invalid_argument("Cannot project an empty discrete hybrid exploration");

    HybridDiscreteTimeUpperOrbitProjection projection{
        Orbit<HybridEnclosure>(exploration.nodes.front().enclosure), {}};

    // The initial set is stored only in orbit.initial(), except when it is
    // itself a terminal (zero horizon or blocked at the first index).
    for (SizeType position = 1u; position < exploration.nodes.size(); ++position) {
        const HybridDiscreteTimeUpperNode& node = exploration.nodes[position];
        projection.orbit.adjoin_reach(node.enclosure);
        if (node.index < exploration.termination_index)
            projection.orbit.adjoin_intermediate(node.enclosure);
    }
    for (const HybridDiscreteTimeUpperTerminal& terminal: exploration.terminals) {
        if (terminal.node_index >= exploration.nodes.size())
            throw std::invalid_argument("Discrete hybrid terminal references an invalid node");
        const HybridDiscreteTimeUpperNode& node = exploration.nodes[terminal.node_index];
        switch (terminal.kind) {
            case HybridDiscreteTimeUpperTerminalKind::HORIZON:
                if (node.index != exploration.termination_index)
                    throw std::invalid_argument("Discrete hybrid horizon terminal has wrong index");
                projection.orbit.adjoin_final(node.enclosure);
                break;
            case HybridDiscreteTimeUpperTerminalKind::BLOCKED:
                if (node.index >= exploration.termination_index)
                    throw std::invalid_argument("Discrete hybrid blocked terminal is not before horizon");
                projection.blocked.adjoin(node.enclosure);
                break;
            default:
                throw std::invalid_argument("Unknown discrete hybrid terminal kind");
        }
    }
    return projection;
}

} // namespace Ariadne
