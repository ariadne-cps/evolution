/***************************************************************************
 *            hybrid_inclusion_invariants.hpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#ifndef ARIADNE_EVOLUTION_HYBRID_INCLUSION_INVARIANTS_HPP
#define ARIADNE_EVOLUTION_HYBRID_INCLUSION_INVARIANTS_HPP

#include "evolution/hybrid_continuous_time_evolver.hpp"

namespace Ariadne::EvolutionDetail {

/** @brief Admit state invariants, PROGRESS blockers and conservative
 *  URGENT/PERMISSIVE candidates under UPPER only. Neither first entry nor impact direction is
 *  certified for arbitrary measurable inputs.
 */
Void validate_inclusion_event_kinds(const HybridAutomatonInterface&, DiscreteLocation);

/** @brief Extract state guards and reset data without evaluating an ODAI Lie
 *  derivative. This intentionally does not claim first contact or crossing.
 */
Map<DiscreteEvent, TransitionData> inclusion_upper_transitions(
    const HybridAutomatonInterface&, DiscreteLocation);

/** @brief Intersect a state enclosure with necessary invariant conditions.
 *  Applied to both reach and endpoint of an UPPER inclusion step. This is
 *  a sound outer filter, not a certification of pathwise invariant satisfaction
 *  and not a first-contact or event-time calculation.
 */
Void apply_state_invariants(HybridEnclosure&, const Map<DiscreteEvent, TransitionData>&);

/** @brief Constrain an UPPER reach/endpoint to the admissible side of
 *  invariants and PROGRESS boundaries. Unlike invariants, PROGRESS is
 *  not imposed on already-terminal zero-horizon initial states.
 */
Void apply_inclusion_flow_blockers(HybridEnclosure&, const Map<DiscreteEvent, TransitionData>&);

/** @brief Publish all possible UPPER reset candidates from one validated reach
 *  enclosure. Candidates retain the source time/dwell map until reset;
 *  no candidate is claimed to cover all source trajectories.
 */
Void publish_inclusion_jump_candidates(
    const HybridEnclosure& reach, const Map<DiscreteEvent, TransitionData>&,
    const Real& final_time, HybridEvolutionStepResult&);

} // namespace Ariadne::EvolutionDetail
#endif
