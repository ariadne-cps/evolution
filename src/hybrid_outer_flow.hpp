/***************************************************************************
 *            hybrid_outer_flow.hpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#ifndef ARIADNE_EVOLUTION_HYBRID_OUTER_FLOW_HPP
#define ARIADNE_EVOLUTION_HYBRID_OUTER_FLOW_HPP

#include "evolution/hybrid_continuous_time_evolver.hpp"

namespace Ariadne::EvolutionDetail {
struct AppliedIntegrationStep;
/** @brief Common UPPER guard filtering/publication for ordinary outer steps.
 *  Equation and inclusion producers retain their own numerical contracts.
 */
HybridEvolutionStepResult publish_outer_step(const HybridEnclosure&, const AppliedIntegrationStep&,
                                             const Map<DiscreteEvent, TransitionData>&, const Real& final_time, double processed_time);

/** @brief Admit state invariants, PROGRESS blockers and conservative
 *  URGENT/PERMISSIVE/IMPACT candidates under UPPER only. Neither first entry nor impact direction is
 *  certified for arbitrary measurable inputs.
 */
Void validate_outer_event_kinds(const HybridAutomatonInterface&, DiscreteLocation);

/** @brief Extract state/time guards and reset data without evaluating a field Lie
 *  derivative. This intentionally does not claim first contact or crossing.
 */
Map<DiscreteEvent, TransitionData> outer_flow_transitions(
    const HybridAutomatonInterface&, DiscreteLocation);

/** @brief Intersect a state enclosure with necessary invariant conditions.
 *  Applied to both reach and endpoint of an UPPER outer step. This is
 *  a sound outer filter, not a certification of pathwise invariant satisfaction
 *  and not a first-contact or event-time calculation.
 */
Void apply_state_invariants(HybridEnclosure&, const Map<DiscreteEvent, TransitionData>&);

/** @brief Constrain UPPER reach/endpoint by all blocking event guards:
 *  invariants, PROGRESS, URGENT and IMPACT, but never PERMISSIVE.
 *  This is a necessary outer restriction, not certified first-entry timing.
 *  Unlike invariants, transition and PROGRESS guards are not imposed
 *  retroactively on an already-terminal zero-horizon initial state.
 */
Void apply_outer_flow_blockers(HybridEnclosure&, const Map<DiscreteEvent, TransitionData>&);

/** @brief Publish all possible UPPER reset candidates from one validated reach
 *  enclosure. Candidates retain the source time/dwell map until reset;
 *  no candidate is claimed to cover all source trajectories. IMPACT uses the
 *  closed guard boundary without directional exclusions, retaining tangencies
 *  and possibly spurious contacts under UPPER.
 */
Void publish_outer_jump_candidates(
    const HybridEnclosure& reach, const Map<DiscreteEvent, TransitionData>&,
    const Real& final_time, HybridEvolutionStepResult&);

} // namespace Ariadne::EvolutionDetail
#endif
