/***************************************************************************
 *            hybrid_outer_flow.cpp
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

/** @file src/hybrid_outer_flow.cpp
 *  @brief Validated hybrid outer flow and certified transition publication.
 */

#include "hybrid_outer_flow.hpp"

#include <stdexcept>

#include "continuous_integration.hpp"
#include "hybrid_parameters.hpp"
#include "hybrid_transition_handlers.hpp"
#include "system/continuous_time_hybrid_automaton_interface.hpp"

namespace Ariadne::EvolutionDetail {

bool try_publish_lower_urgent_crossing(
    const HybridEnclosure& initial, const AppliedIntegrationStep& applied,
    const Map<DiscreteEvent, TransitionData>& transitions,
    const Real& final_time, double processed_time, HybridEvolutionStepResult& result) {
    // A single certified tube gives a constructive LOWER first-crossing
    // witness for an URGENT transition without solving for an exact
    // crossing time. A real trajectory starts strictly before the guard,
    // and its hypothetical uninterrupted endpoint is strictly beyond it.
    // Continuity gives a first zero; excluding every other blocker on
    // the *whole tube* makes that first zero an admissible hybrid jump.
    // The guard-equality slice of the tube contains this real witness,
    // even though it need not be an INNER set.
    if (applied.reach.size() == 1u && !transitions.empty()) {
        HybridEnclosure tube = initial;
        tube.replace_continuous_set(applied.reach[0u]);
        HybridEnclosure endpoint = initial;
        endpoint.replace_continuous_set(applied.evolved);
        // With uncertain arrival times a step bounded from the earliest
        // source time can end after the requested horizon for the actual
        // LOWER subset. Do not pick earlier outer times as witnesses.
        if (!definitely(endpoint.time_range().upper_bound() <= final_time))
            return false;
        for (const auto& entry: transitions) {
            const auto& crossing = entry.second;
            if (crossing.event_kind != EventKind::URGENT)
                continue;
            const auto at_start = initial.range_of(crossing.guard_function, crossing.source_arguments);
            const auto at_finish = endpoint.range_of(crossing.guard_function, crossing.source_arguments);
            // A certified nonnegative endpoint also guarantees a real
            // first zero (possibly at the step endpoint) by continuity.
            // Strict positivity would unnecessarily exclude closed
            // boundary contacts where the guard is exactly zero.
            if (!definitely(at_start.upper_bound() < 0) ||
                !definitely(at_finish.lower_bound() >= 0))
                continue;
            bool no_competing_blocker = true;
            for (const auto& other_entry: transitions) {
                if (other_entry.first == entry.first)
                    continue;
                const auto& other = other_entry.second;
                if (!transition_policy(other.event_kind).blocks_flow)
                    continue;
                const auto range = tube.range_of(other.guard_function, other.source_arguments);
                const bool safe = other.event_kind == EventKind::INVARIANT
                                    ? definitely(range.upper_bound() <= 0)
                                    : definitely(range.upper_bound() < 0);
                if (!safe) {
                    no_competing_blocker = false;
                    break;
                }
            }
            if (!no_competing_blocker)
                continue;

            HybridEnclosure witnessed_crossing = tube;
            witnessed_crossing.new_guard(crossing.event,
                                         crossing.guard_function,
                                         crossing.source_arguments);
            if (definitely(witnessed_crossing.is_empty()))
                continue;
            result.record_reach(witnessed_crossing);
            apply_transition_reset(witnessed_crossing, crossing);
            result.record_intermediate(witnessed_crossing);
            result.transition_to(witnessed_crossing, crossing.event, crossing.event_kind);
            // One certified integration step furnished this crossing.
            // The precise first-entry time is deliberately not guessed.
            result.continuous_evolution_time = processed_time;
            return true;
        }
    }
    return false;
}

Void validate_outer_event_kinds(const ContinuousTimeHybridAutomatonInterface& system, DiscreteLocation location) {
    for (const auto& event: system.events(location)) {
        const auto kind = system.event_kind(location, event);
        if (kind != EventKind::INVARIANT && kind != EventKind::PROGRESS &&
            kind != EventKind::URGENT && kind != EventKind::PERMISSIVE &&
            kind != EventKind::IMPACT)
            throw std::invalid_argument(
                "Hybrid UPPER outer flow supports state invariants, PROGRESS blockers and "
                "conservative URGENT/PERMISSIVE/IMPACT jumps");
    }
}

Map<DiscreteEvent, TransitionData> outer_flow_transitions(
    const ContinuousTimeHybridAutomatonInterface& system, DiscreteLocation location) {
    validate_outer_event_kinds(system, location);
    Map<DiscreteEvent, TransitionData> transitions;
    for (const auto& event: system.events(location)) {
        TransitionData transition;
        transition.event = event;
        transition.event_kind = system.event_kind(location, event);
        transition.source_arguments = hybrid_function_arguments(system, location);
        if (transition.event_kind == EventKind::INVARIANT) {
            transition.guard_function = system.invariant_function(location, event);
        } else {
            transition.guard_function = system.guard_function(location, event);
            // A PROGRESS guard blocks continuous flow; it is not a reset.
            if (transition_policy(transition.event_kind).activates_reset) {
                transition.target = system.target(location, event);
                transition.target_space = system.continuous_state_space(transition.target);
                transition.reset_function = system.reset_function(location, event);
                transition.target_auxiliary_space = system.continuous_auxiliary_space(transition.target);
                transition.target_auxiliary_function = system.auxiliary_function(transition.target);
                configure_transition_target(system, transition);
            }
        }
        transitions.insert(event, transition);
    }
    return transitions;
}

Void apply_state_invariants(HybridEnclosure& set,
                            const Map<DiscreteEvent, TransitionData>& transitions) {
    for (const auto& item: transitions) {
        const DiscreteEvent event = item.first;
        const TransitionData& transition = item.second;
        if (transition.event_kind == EventKind::INVARIANT &&
            possibly(set.range_of(transition.guard_function, transition.source_arguments).upper_bound() >= 0)) {
            set.new_invariant(event, transition.guard_function, transition.source_arguments);
        }
    }
}

Void apply_outer_flow_blockers(
    HybridEnclosure& set, const Map<DiscreteEvent, TransitionData>& transitions) {
    apply_state_invariants(set, transitions);
    for (const auto& item: transitions) {
        const auto& transition = item.second;
        // This must match the ODAE event policy: a validated UPPER
        // continuous branch cannot lie strictly beyond an URGENT or
        // PROGRESS/IMPACT boundary before taking a transition. PERMISSIVE
        // never blocks flow. Guard intersections remain outer filters.
        if (transition.event_kind != EventKind::INVARIANT &&
            transition_policy(transition.event_kind).blocks_flow &&
            possibly(set.range_of(transition.guard_function, transition.source_arguments).upper_bound() >= 0)) {
            set.new_invariant(transition.event, transition.guard_function, transition.source_arguments);
        }
    }
}

Void publish_outer_jump_candidates(
    const HybridEnclosure& reach, const Map<DiscreteEvent, TransitionData>& transitions,
    const Real& final_time, HybridEvolutionStepResult& result) {
    // The validated outer reach contains each feasible trajectory state/time
    // at which a guard can activate. Its intersection with g>=0 can include
    // spurious candidates, which is allowed for UPPER. Neither a first entry
    // nor a full-source urgent transition is inferred from such an intersection.
    for (const auto& item: transitions) {
        const auto& transition = item.second;
        if (!transition_policy(transition.event_kind).activates_reset)
            continue;
        HybridEnclosure candidate = reach;
        if (transition.event_kind == EventKind::IMPACT) {
            // Every admissible impact, including a zero-derivative contact,
            // lies on g=0. Retain the entire boundary intersection: arbitrary
            // measurable inputs need not supply a unique directional value.
            // Inward or repeated contacts may remain spurious UPPER jumps.
            candidate.new_guard(transition.event, transition.guard_function, transition.source_arguments);
        } else {
            candidate.new_activation(transition.event, transition.guard_function, transition.source_arguments);
        }
        if (possibly(candidate.time_range().upper_bound() > final_time))
            candidate.bound_time(final_time);
        if (definitely(candidate.is_empty()))
            continue;
        apply_transition_reset(candidate, transition);
        result.transition_to(candidate, transition.event, transition.event_kind);
    }
}

HybridEvolutionStepResult publish_outer_step(const HybridEnclosure& initial,
                                             const AppliedIntegrationStep& applied, const Map<DiscreteEvent, TransitionData>& invariants,
                                             const Real& final_time, double processed_time) {
    HybridEvolutionStepResult result;
    HybridEnclosure endpoint = initial;
    endpoint.replace_continuous_set(applied.evolved);
    // Necessary state restrictions include invariants, PROGRESS, URGENT and IMPACT
    // blockers. They do not certify an inclusion crossing time.
    EvolutionDetail::apply_outer_flow_blockers(endpoint, invariants);
    const auto endpoint_time = endpoint.time_range();
    const bool at_horizon = definitely(endpoint_time.lower_bound() >= final_time) &&
                            definitely(endpoint_time.upper_bound() <= final_time);
    HybridEvolutionStepResult completed;
    for (const auto& continuous: applied.reach) {
        auto reach = initial;
        reach.replace_continuous_set(continuous);
        // Restrict the certified outer tube by invariant and every
        // blocking guard, including competing URGENT events. The guard
        // boundary is retained for conservative event candidates.
        EvolutionDetail::apply_outer_flow_blockers(reach, invariants);
        // Publish the pre-reset guard intersection, retaining the full
        // source-time chart. The reset is a separate successor; UPPER may
        // include infeasible (spurious) candidate transitions.
        EvolutionDetail::publish_outer_jump_candidates(reach, invariants, final_time, completed);
        if (!at_horizon && possibly(reach.time_range().upper_bound() >= final_time)) {
            auto final_slice = reach;
            final_slice.set_time(final_time);
            if (!definitely(final_slice.is_empty())) {
                // Keep the old time equality as a constraint on the state and
                // dwell maps, then expose the certified stopping time itself.
                final_slice.set_time_function(initial.configuration().function_factory().create_constant(
                    final_slice.parameter_domain(), final_time));
                completed.record_final(final_slice);
            }
        }
        if (possibly(reach.time_range().upper_bound() > final_time))
            reach.bound_time(final_time);
        if (!definitely(reach.is_empty()))
            completed.record_reach(reach);
    }
    if (at_horizon) {
        endpoint.set_time(final_time);
        if (!definitely(endpoint.is_empty()))
            completed.record_final(endpoint);
    } else if (possibly(endpoint_time.lower_bound() < final_time)) {
        if (possibly(endpoint_time.upper_bound() > final_time))
            endpoint.bound_time(final_time);
        if (!definitely(endpoint.is_empty())) {
            completed.record_intermediate(endpoint);
            completed.continue_with(endpoint);
        }
    }
    completed.continuous_evolution_time = processed_time;
    if (completed.successors.empty())
        completed.stop_reason = HybridEvolutionStopReason::TIME_LIMIT;
    result = std::move(completed);
    return result;
}

} // namespace Ariadne::EvolutionDetail
