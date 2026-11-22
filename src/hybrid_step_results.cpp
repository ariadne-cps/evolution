/***************************************************************************
 *            hybrid_step_results.cpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#include "hybrid_step_results.hpp"

#include "logging/logging.hpp"
#include "numeric/numeric.hpp"

namespace Ariadne::EvolutionDetail {

HybridEndpointBranches pair_endpoint_branches(const List<HybridEnclosure>& endpoints) {
    HybridEndpointBranches branches;
    branches.reserve(endpoints.size());
    for (const auto& endpoint: endpoints)
        branches.push_back({endpoint, endpoint});
    return branches;
}

Void restrict_reported_endpoints(HybridEndpointBranches& branches, const HybridGuardContext& context,
                                 const ValidatedScalarMultivariateFunctionPatch& elapsed_time) {
    const SizeType original_size = branches.size();
    for (SizeType i = 0u; i != original_size; ++i) {
        List<HybridEnclosure> restricted{branches[i].reported};
        apply_flow_guard(restricted, context, elapsed_time);
        branches[i].reported = restricted.front();
        // Snapshot before append: the vector may reallocate when a guard splits.
        const HybridEnclosure checkpoint = branches[i].continuation;
        for (SizeType j = 1u; j != restricted.size(); ++j)
            branches.push_back({restricted[j], checkpoint});
    }
}

HybridStepPublication publish_step_contributions(
    List<HybridEnclosure>& reach_sets, HybridEndpointBranches& branches,
    const TimingData& timing_data, HybridEvolutionStepResult& result) {
    LOGGING_SCOPE_CREATE;
    HybridStepPublication publication;
    if (timing_data.finishing_kind == FinishingKind::AT_FINAL_TIME) {
        for (const auto& branch: branches) {
            const HybridEnclosure& evolve_set = branch.reported;
            if (!definitely(evolve_set.is_empty())) {
                LOGGING_PRINTLN_AT(1, "final_set=" << evolve_set);
                result.record_final(evolve_set);
                publication.finishing = true;
            }
        }
        for (List<HybridEnclosure>::ConstIterator reach_set_iter = reach_sets.begin();
             reach_set_iter != reach_sets.end(); ++reach_set_iter) {
            const HybridEnclosure& reach_set = *reach_set_iter;
            result.record_reach(reach_set);
        }
    } else { // (timing_data.finishing_kind!=AT_FINAL_TIME)
        for (auto& branch: branches) {
            HybridEnclosure& evolve_set = branch.reported;
            HybridEnclosure& next_working_set = branch.continuation;
            UpperIntervalType evolve_set_time_range = evolve_set.time_range();
            if (definitely(evolve_set_time_range.lower_bound() > timing_data.final_time)) {
                // Do nothing, since evolve set is past final time is definitely empty
            } else if (definitely(evolve_set_time_range.upper_bound() <= timing_data.final_time)) {
                // No need to introduce timing constraints
                if (!definitely(evolve_set.is_empty())) {
                    LOGGING_PRINTLN_AT(1, "evolve_set=" << evolve_set);
                    LOGGING_PRINTLN_AT(1, "next_working_set=" << next_working_set);
                    result.continue_with(next_working_set);
                    result.record_intermediate(evolve_set);
                    publication.progress = true;
                }
            } else {
                // Bound time if necessary
                // This arm is the negation of the preceding definite <=
                // test, hence exceeding final_time is already possible.
                evolve_set.bound_time(timing_data.final_time);
                // Only continue evolution if the time-bounded evolve set is nonempty;
                // However, continue evolution without adding time constraint,
                // since this will be introduced in the next step
                if (!definitely(evolve_set.is_empty())) {
                    LOGGING_PRINTLN_AT(1, "evolve_set=" << evolve_set);
                    LOGGING_PRINTLN_AT(1, "next_working_set=" << next_working_set);
                    result.continue_with(next_working_set);
                    result.record_intermediate(evolve_set);
                    publication.progress = true;
                }
            }
        }
        for (List<HybridEnclosure>::Iterator reach_set_iter = reach_sets.begin();
             reach_set_iter != reach_sets.end(); ++reach_set_iter) {
            HybridEnclosure& reach_set = *reach_set_iter;
            // reach_set.continuous_set().reduce();
            //  The closed reach interval also contains the final time at its upper endpoint.
            if (possibly(reach_set.time_range().upper_bound() >= FloatDPBounds(timing_data.final_time, dp))) {
                HybridEnclosure final_set = reach_set;
                final_set.set_time(timing_data.final_time);
                if (timing_data.finishing_kind != FinishingKind::BEFORE_FINAL_TIME && !definitely(final_set.is_empty())) {
                    // if(!definitely(final_set.is_empty())) {
                    LOGGING_PRINTLN_AT(1, "final_set=" << final_set);
                    result.record_final(final_set);
                    publication.finishing = true;
                }
            }
            reach_set.bound_time(timing_data.final_time);
            result.record_reach(reach_set);
        }
    }

    return publication;
}

} // namespace Ariadne::EvolutionDetail
