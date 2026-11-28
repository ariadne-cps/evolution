/***************************************************************************
 *            hybrid_outer_flow.hpp
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

/** @file src/hybrid_outer_flow.hpp
 *  @brief Validated hybrid outer flow and certified transition publication.
 */

#ifndef ARIADNE_EVOLUTION_HYBRID_OUTER_FLOW_HPP
#define ARIADNE_EVOLUTION_HYBRID_OUTER_FLOW_HPP

#include "evolution/hybrid_continuous_time_evolver.hpp"

namespace Ariadne::EvolutionDetail {
struct AppliedIntegrationStep;
/** @brief Sufficient whole-enclosure certificate for an uninterrupted LOWER flow.
 *  Invariants require g<=0; URGENT/PROGRESS/IMPACT require g<0.
 *  Optional PERMISSIVE resets never block the no-jump behavior.
 */
bool lower_flow_predicate_certified(const HybridEnclosure&, const TransitionData&);

/** @brief Apply the common LOWER flow certificate to every complete reach tube. */
bool lower_flow_step_certified(
    const HybridEnclosure&, const AppliedIntegrationStep&, const Map<DiscreteEvent, TransitionData>&);

/** @brief Append all certified PERMISSIVE endpoint resets, retaining no-jump outputs.
 *  The caller must first certify the entire uninterrupted flow prefix.
 *  Only FULL_SOURCE endpoint guards and whole-image target invariants admit
 *  a reset, and all endpoint times must lie within the final horizon.
 *  No source/target restriction selects a reached subset; every reset copies
 *  the complete endpoint. The shared continuous prefix is credited once.
 */
Void publish_lower_permissive_endpoint_resets(
    const ContinuousTimeHybridAutomatonInterface&, const HybridEnclosure&, const Enclosure& evolved,
    const Map<DiscreteEvent, TransitionData>&, const Real& final_time, HybridEvolutionStepResult&);

/** @brief Publish a LOWER URGENT crossing only with universal certificates.
 *  Requires one complete structured tube, a strictly inactive source guard,
 *  a nonnegative whole endpoint guard, an endpoint entirely within the
 *  final horizon and exclusion of every competing blocker over the tube.
 *  The common reset retains the nonempty
 *  reached first-crossing subset, including with measurable inputs.
 *  False leaves result unchanged; IMPACT direction is not inferred.
 */
bool try_publish_lower_urgent_crossing(
    const HybridEnclosure&, const AppliedIntegrationStep&,
    const Map<DiscreteEvent, TransitionData>&, const Real& final_time, double processed_time,
    HybridEvolutionStepResult&);

/** @brief Common outer time filtering/publication for ordinary numerical steps.
 *  LOWER callers certify the flow independently and pass no transition
 *  candidates; certified optional resets are appended separately.
 *  Equation and inclusion producers retain their own numerical contracts.
 */
HybridEvolutionStepResult publish_outer_step(const HybridEnclosure&, const AppliedIntegrationStep&,
                                             const Map<DiscreteEvent, TransitionData>&, const Real& final_time, double processed_time);

/** @brief Admit event kinds handled by the outer-flow bridge.
 *  UPPER may produce conservative candidates. LOWER requires independent
 *  flow/reset certificates; event-kind admission alone proves no crossing
 *  or impact direction for arbitrary measurable inputs.
 */
Void validate_outer_event_kinds(const ContinuousTimeHybridAutomatonInterface&, DiscreteLocation);

/** @brief Extract state/time guards and reset data without evaluating a field Lie
 *  derivative. This intentionally does not claim first contact or crossing.
 */
Map<DiscreteEvent, TransitionData> outer_flow_transitions(
    const ContinuousTimeHybridAutomatonInterface&, DiscreteLocation);

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
