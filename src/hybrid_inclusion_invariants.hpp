/***************************************************************************
 *            hybrid_inclusion_invariants.hpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#ifndef ARIADNE_EVOLUTION_HYBRID_INCLUSION_INVARIANTS_HPP
#define ARIADNE_EVOLUTION_HYBRID_INCLUSION_INVARIANTS_HPP

#include "evolution/hybrid_continuous_time_evolver.hpp"

namespace Ariadne::EvolutionDetail {

/** @brief Admit only state invariants for an autonomous UPPER inclusion mode.
 *  Arbitrary measurable inputs offer no certified first-entry, urgent,
 *  progress or impact timing in the current hybrid integration bridge.
 */
Void validate_inclusion_event_kinds(const HybridAutomatonInterface&, DiscreteLocation);

/** @brief Build invariant-only data without evaluating a nonexistent ODAI vector field.
 *  The unused Lie-derivative field and reset data remain empty.
 */
Map<DiscreteEvent, TransitionData> inclusion_invariant_transitions(
    const HybridAutomatonInterface&, DiscreteLocation);

/** @brief Intersect a state enclosure with necessary invariant conditions.
 *  Applied to both reach and endpoint of an UPPER inclusion step. This is
 *  a sound outer filter, not a certification of pathwise invariant satisfaction
 *  and not a first-contact or event-time calculation.
 */
Void apply_state_invariants(HybridEnclosure&, const Map<DiscreteEvent, TransitionData>&);

} // namespace Ariadne::EvolutionDetail
#endif
