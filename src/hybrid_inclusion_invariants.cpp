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
        if (system.event_kind(location, event) != EventKind::INVARIANT)
            throw std::invalid_argument(
                "Hybrid inclusion modes support only state invariants in UPPER semantics; "
                "transitions, progress predicates and impact events require certified inclusion event analysis");
    }
}

Map<DiscreteEvent, TransitionData> inclusion_invariant_transitions(
    const HybridAutomatonInterface& system, DiscreteLocation location) {
    validate_inclusion_event_kinds(system, location);
    Map<DiscreteEvent, TransitionData> transitions;
    for (const auto& event: system.events(location)) {
        TransitionData transition;
        transition.event = event;
        transition.event_kind = EventKind::INVARIANT;
        transition.guard_function = system.invariant_function(location, event);
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

} // namespace Ariadne::EvolutionDetail
