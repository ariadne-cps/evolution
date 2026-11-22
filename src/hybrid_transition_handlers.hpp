/***************************************************************************
 *            hybrid_transition_handlers.hpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#ifndef ARIADNE_EVOLUTION_HYBRID_TRANSITION_HANDLERS_HPP
#define ARIADNE_EVOLUTION_HYBRID_TRANSITION_HANDLERS_HPP

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

/** @brief Preserve the continuation copy before these endpoint restrictions. */
bool defer_endpoint_guard(CrossingKind);

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
