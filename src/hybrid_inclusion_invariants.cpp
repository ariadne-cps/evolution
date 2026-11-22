/***************************************************************************
 *            hybrid_inclusion_invariants.cpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#include "hybrid_inclusion_invariants.hpp"

#include <stdexcept>

#include "system/hybrid_automaton_interface.hpp"

namespace Ariadne::EvolutionDetail {

Void validate_inclusion_event_kinds(const HybridAutomatonInterface& system, DiscreteLocation location) {
    for (const auto& event: system.events(location)) {
        const auto kind = system.event_kind(location, event);
        if (kind != EventKind::INVARIANT && kind != EventKind::PROGRESS &&
            kind != EventKind::URGENT && kind != EventKind::PERMISSIVE)
            throw std::invalid_argument(
                "Hybrid inclusion UPPER currently supports state invariants, PROGRESS blockers and "
                "conservative URGENT/PERMISSIVE jumps; IMPACT requires dedicated event certificates");
    }
}

Map<DiscreteEvent, TransitionData> inclusion_upper_transitions(
    const HybridAutomatonInterface& system, DiscreteLocation location) {
    validate_inclusion_event_kinds(system, location);
    Map<DiscreteEvent, TransitionData> transitions;
    for (const auto& event: system.events(location)) {
        TransitionData transition;
        transition.event = event;
        transition.event_kind = system.event_kind(location, event);
        if (transition.event_kind == EventKind::INVARIANT) {
            transition.guard_function = system.invariant_function(location, event);
        } else {
            transition.guard_function = system.guard_function(location, event);
            // A PROGRESS guard blocks continuous flow; it is not a reset.
            if (transition.event_kind == EventKind::URGENT ||
                transition.event_kind == EventKind::PERMISSIVE) {
                transition.target = system.target(location, event);
                transition.target_space = system.continuous_state_space(transition.target);
                transition.reset_function = system.reset_function(location, event);
                transition.target_auxiliary_space = system.continuous_auxiliary_space(transition.target);
                transition.target_auxiliary_function = system.auxiliary_function(transition.target);
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
            possibly(set.satisfies(transition.guard_function >= 0))) {
            set.new_invariant(event, transition.guard_function);
        }
    }
}

Void apply_inclusion_flow_blockers(
    HybridEnclosure& set, const Map<DiscreteEvent, TransitionData>& transitions) {
    apply_state_invariants(set, transitions);
    for (const auto& item: transitions) {
        const auto& transition = item.second;
        if (transition.event_kind == EventKind::PROGRESS &&
            possibly(set.satisfies(transition.guard_function >= 0))) {
            set.new_invariant(transition.event, transition.guard_function);
        }
    }
}

Void publish_inclusion_jump_candidates(
    const HybridEnclosure& reach, const Map<DiscreteEvent, TransitionData>& transitions,
    const Real& final_time, HybridEvolutionStepResult& result) {
    // The validated outer reach contains each feasible trajectory state/time
    // at which a guard can activate. Its intersection with g>=0 can include
    // spurious candidates, which is allowed for UPPER. Neither a first entry
    // nor a full-source urgent transition is inferred from such an intersection.
    for (const auto& item: transitions) {
        const auto& transition = item.second;
        if (transition.event_kind != EventKind::URGENT &&
            transition.event_kind != EventKind::PERMISSIVE)
            continue;
        HybridEnclosure candidate = reach;
        candidate.new_activation(transition.event, transition.guard_function);
        if (possibly(candidate.time_range().upper_bound() > final_time))
            candidate.bound_time(final_time);
        if (definitely(candidate.is_empty()))
            continue;
        candidate.apply_reset(transition.event, transition.target, transition.target_space,
                              transition.reset_function, transition.target_auxiliary_space,
                              transition.target_auxiliary_function);
        result.transition_to(candidate, transition.event, transition.event_kind);
    }
}

} // namespace Ariadne::EvolutionDetail
