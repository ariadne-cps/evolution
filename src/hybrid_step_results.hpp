/***************************************************************************
 *            hybrid_step_results.hpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#ifndef ARIADNE_EVOLUTION_HYBRID_STEP_RESULTS_HPP
#define ARIADNE_EVOLUTION_HYBRID_STEP_RESULTS_HPP

#include <vector>

#include "evolution/hybrid_continuous_time_evolver.hpp"
#include "hybrid_transition_handlers.hpp"

namespace Ariadne::EvolutionDetail {

/** @brief An endpoint's reported constraints and its continuation checkpoint.
 *  Deferred blocking guards change only reported, while continuation retains
 *  the enclosure before those restrictions. A geometric split copies the
 *  original checkpoint into every resulting reported branch.
 */
struct HybridEndpointBranch {
    HybridEnclosure reported;
    HybridEnclosure continuation;
};

using HybridEndpointBranches = std::vector<HybridEndpointBranch>;

/** @brief Associate each endpoint with the same pre-deferred checkpoint. */
HybridEndpointBranches pair_endpoint_branches(const List<HybridEnclosure>&);

/** @brief Restrict reported endpoints, retaining origin checkpoints and branch order. */
Void restrict_reported_endpoints(HybridEndpointBranches&, const HybridGuardContext&,
                                 const ValidatedScalarMultivariateFunctionPatch& elapsed_time);

/** @brief Local publication flags, independent of event/reset processing. */
struct HybridStepPublication {
    bool progress = false;
    bool finishing = false;
};

/** @brief Publish reach, reported endpoints, and continuation checkpoints.
 *  Preserve the closed horizon, final-versus-intermediate rules and output order.
 *  The caller publishes jumps only after this phase.
 */
HybridStepPublication publish_step_contributions(
    List<HybridEnclosure>& reach_sets, HybridEndpointBranches& branches,
    const TimingData&, HybridEvolutionStepResult&);

} // namespace Ariadne::EvolutionDetail
#endif
