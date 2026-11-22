/***************************************************************************
 *            hybrid_transition_handlers.hpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#ifndef ARIADNE_EVOLUTION_HYBRID_TRANSITION_HANDLERS_HPP
#define ARIADNE_EVOLUTION_HYBRID_TRANSITION_HANDLERS_HPP

#include <cstdint>

#include "evolution/hybrid_continuous_time_evolver.hpp"

namespace Ariadne::EvolutionDetail {

/** @brief Event policy independent of crossing geometry and numerical solver. */
struct HybridTransitionPolicy {
    bool blocks_flow;
    bool activates_reset;
    bool inspect_initial_activation;
    bool resets_initially;
    bool requires_impact_direction;
};

HybridTransitionPolicy transition_policy(EventKind);

/** @brief Validated initial guard activation on one entire invariant-valid source.
 *  UNRESOLVED can mean partial activation or insufficient numerical evidence;
 *  it never certifies a LOWER transition. This is not a collective tree check.
 */
enum class HybridInitialGuardActivation : std::uint8_t {
    EXCLUDED,   //!< Intersection with the active side is definitely empty.
    UNRESOLVED, //!< A candidate may activate, but not the whole source is certified.
    FULL_SOURCE //!< The entire represented source satisfies the guard.
};

/** @brief Use guard intersection for possibility and universal satisfaction for coverage. */
HybridInitialGuardActivation classify_initial_guard(
    const HybridEnclosure& invariant_valid_source, const TransitionData& transition);

/** @brief UPPER admits possible jumps; LOWER only certifiably whole-source jumps.
 *  For split trajectories, a future collective check must cover every piece.
 */
bool permits_initial_reset(HybridInitialGuardActivation, Semantics);

/** @brief Private flow-restriction implementations, selected by actual evidence. */
enum class HybridFlowGuardMethod : std::uint8_t {
    KEEP,
    FIRST_ENTRY_BOUND,
    GUARD_SIDE,
    TWO_SAMPLES,
    MAXIMUM_SPLIT,
    DISCARD
};

/** @brief Private pre-reset jump implementations, separate from flow guards. */
enum class HybridJumpGuardMethod : std::uint8_t {
    SOURCE_ACTIVATION,
    TUBE_ACTIVATION,
    FIRST_ENTRY,
    BOUNDARY,
    DIRECTED_BOUNDARY,
    DISCARD
};

/** @brief The numerical evidence required by the chosen flow operation.
 *  All calls from the evolver require a validated autonomous local state-flow
 *  patch. NONE is not a claim that the source geometry is certified.
 */
struct HybridFlowGuardCapability {
    HybridFlowGuardMethod method;
    CrossingEvidenceKind evidence;
    bool uses_timed_flow;
};

/** @brief Evidence and additional analytic operations required by a jump.
 *  In particular, a directed entry or an impact uses a Lie derivative; the
 *  latter condition is independent of the chosen entry operation.
 */
struct HybridJumpGuardCapability {
    HybridJumpGuardMethod method;
    CrossingEvidenceKind evidence;
    bool uses_timed_flow;
    bool requires_lie_derivative;
};

/** @brief Classify actual flow strategy; semantics do not imply inner guarantees. */
HybridFlowGuardCapability flow_guard_capability(const CrossingData&, Semantics);

/** @brief Classify actual jump strategy and its analytic requirements. */
HybridJumpGuardCapability jump_guard_capability(EventKind, const CrossingData&, Semantics);

/** @brief Preserve the continuation copy before these endpoint restrictions. */
bool defer_endpoint_guard(const CrossingData&);

/** @brief Borrowed immutable numerical data; no workload or backend is owned. */
struct HybridGuardContext {
    const HybridEnclosure& initial;
    const ValidatedVectorMultivariateFunctionPatch& flow;
    const TransitionData& transition;
    const CrossingData& crossing;
    Semantics semantics;
    ValidatedFunctionPatchFactory factory;
};

/** @brief A jump candidate before reset, paired with its actual elapsed map. */
struct HybridJumpCandidate {
    HybridEnclosure enclosure;
    ValidatedScalarMultivariateFunctionPatch elapsed_time;
};

HybridJumpCandidate make_jump_candidate(const HybridGuardContext&,
                                        const EffectiveVectorMultivariateFunction& dynamic,
                                        const TimingData&);

/** @brief Apply a blocking guard with stable primary/additional branch order. */
Void apply_flow_guard(List<HybridEnclosure>&, const HybridGuardContext&,
                      const ValidatedScalarMultivariateFunctionPatch& elapsed_time);

} // namespace Ariadne::EvolutionDetail
#endif
