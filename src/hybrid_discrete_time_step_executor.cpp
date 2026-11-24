/***************************************************************************
 *            src/hybrid_discrete_time_step_executor.cpp
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

/** @file src/hybrid_discrete_time_step_executor.cpp
 *  @brief One-tick deterministic hybrid UPPER fan-out with event provenance.
 */

#include "evolution/hybrid_discrete_time_step_executor.hpp"

#include <stdexcept>

namespace Ariadne {

HybridDiscreteTimeStepExecution HybridDiscreteTimeStepExecutor::execute(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source, Integer source_index,
    const HybridDiscreteTimeStepPlanner::EvidenceEvaluator& evaluate) {
    if (!automaton.has_mode(source.location()))
        throw std::invalid_argument("Discrete hybrid step source mode does not exist");
    if (!automaton.has_equation(source.location()))
        throw std::invalid_argument("Validated discrete inclusion images are not implemented");

    const HybridDiscreteTimeStepPlan plan =
        HybridDiscreteTimeStepPlanner::plan(automaton, source.location(), source_index, evaluate);

    HybridDiscreteTimeStepExecution result{plan.source_index, plan.successor_index, {}};
    result.successors.reserve(plan.branches.size());
    for (const HybridDiscreteTimePlannedBranch& branch: plan.branches) {
        // Each image starts from the original enclosure: no event ordering,
        // zero-time transition chain or implicit guard-priority selection.
        // If any image is unsupported the entire execution fails rather than
        // returning an incomplete set of successors.
        result.successors.push_back({
            branch, HybridDiscreteTimeStepImage::apply(automaton, source, branch)});
    }
    return result;
}

} // namespace Ariadne
