/***************************************************************************
 *            hybrid_transition_handlers.cpp
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

/** @file src/hybrid_transition_handlers.cpp
 *  @brief Select transition policies and geometric handlers.
 */

#include "hybrid_transition_handlers.hpp"

#include "function/function.hpp"
#include "hybrid_guard_predicate.hpp"
#include "logging/logging.hpp"

namespace Ariadne::EvolutionDetail {
namespace {
const FloatDP handler_zero(0, dp);
const DiscreteEvent guard_step_event("_h_");

void validate_semantics(Semantics semantics) {
    switch (semantics) {
        case Semantics::UPPER:
        case Semantics::LOWER:
            return;
        default:
            ARIADNE_FAIL_MSG("Unhandled semantics.");
    }
}

struct FlowContext {
    const HybridGuardContext& guard;
    const ValidatedScalarMultivariateFunctionPatch& elapsed;
    ValidatedVectorMultivariateFunctionPatch initial_state;
};
using ExtraBranches = List<HybridEnclosure>;
using FlowHandler = ExtraBranches (*)(HybridEnclosure&, const FlowContext&);

ExtraBranches keep_flow(HybridEnclosure&, const FlowContext&) { return {}; }

ExtraBranches restrict_crossing_time(HybridEnclosure& set, const FlowContext& context) {
    ARIADNE_PRECONDITION(set.parameter_domain() == context.initial_state.domain());
    set.new_parameter_constraint(context.guard.transition.event,
                                 context.elapsed <= unchecked_compose(context.guard.crossing.crossing_time()->time(), context.initial_state));
    return {};
}

ExtraBranches restrict_guard_side(HybridEnclosure& set, const FlowContext& context) {
    set.new_invariant(context.guard.transition.event, context.guard.transition.guard_function);
    return {};
}

ExtraBranches discard_flow(HybridEnclosure& set, const FlowContext& context) {
    set.new_parameter_constraint(context.guard.transition.event,
                                 context.guard.factory.create_constant(set.parameter_domain(), 1) <= handler_zero);
    return {};
}

ExtraBranches restrict_sampled_flow(HybridEnclosure& set, const FlowContext& context) {
    // Necessary conditions at the same two sample times as the original UPPER rule.
    static const Nat subdivisions = 2;
    for (Nat i = 0; i != subdivisions; ++i) {
        const FloatDPBounds alpha = FloatDP(i + 1, dp) / subdivisions;
        const auto intermediate_guard = compose(context.guard.transition.guard_function,
                                                unchecked_compose(context.guard.flow, join(context.initial_state, alpha * context.elapsed)));
        set.new_parameter_constraint(context.guard.transition.event, intermediate_guard <= handler_zero);
    }
    return {};
}

ExtraBranches restrict_maximum_flow(HybridEnclosure& set, const FlowContext& context) {
    const auto event = context.guard.transition.event;
    const auto& guard_function = context.guard.transition.guard_function;
    const auto& guard_crossing_data = context.guard.crossing;
    const auto& starting_state_function = context.initial_state;
    const auto& elapsed_time_function = context.elapsed;
    const auto& flow = context.guard.flow;
    ExtraBranches additional;

    ValidatedScalarMultivariateFunctionPatch critical_time_function = unchecked_compose(guard_crossing_data.guard_maximum()->time(), starting_state_function);
    ValidatedScalarMultivariateFunctionPatch final_guard_function = compose(guard_function, unchecked_compose(flow, join(starting_state_function, elapsed_time_function)));
    ValidatedScalarMultivariateFunctionPatch maximal_guard_function = compose(guard_function, unchecked_compose(flow, join(starting_state_function, critical_time_function)));
    UpperIntervalType guard_range_at_critical_time = guard_crossing_data.guard_maximum()->guard_range();
    // Embedding above gives the starting map the elapsed-time
    // domain; composition and join preserve that argument domain.
    // The per-set domain check already establishes compatibility.

    // If no points in the set arise from trajectories which will later leave the progress set,
    // then we only need to look at the maximum value of the guard.
    HybridEnclosure eventually_hitting_set = set;
    // Test whether the trajectory can hit at any future time
    // from its starting state. An endpoint already beyond the
    // maximum may have crossed and returned; it is not a miss.
    eventually_hitting_set.new_parameter_constraint(event, handler_zero <= critical_time_function);
    eventually_hitting_set.new_parameter_constraint(event, maximal_guard_function >= handler_zero);
    LOGGING_PRINTLN_AT(1, "eventually_hitting_set.is_empty()=" << eventually_hitting_set.is_empty());
    if (definitely(eventually_hitting_set.is_empty())) {
        // NOTE: Constraint below should always be satisfied
        // set.new_parameter_constraint(event, maximal_guard_function <= zero);
        return {};
    }
    // If no points in the set arise from trajectories which leave the progress set and
    // later return, then we only need to look at the guard at the final value
    HybridEnclosure returning_set = set;
    returning_set.new_parameter_constraint(event, elapsed_time_function >= critical_time_function);
    returning_set.new_parameter_constraint(event, final_guard_function <= handler_zero);
    ValidatedLowerKleenean returning_set_empty = returning_set.is_empty();
    LOGGING_PRINTLN_AT(1, "returning_set.is_empty()=" << returning_set_empty);
    if (definitely(returning_set_empty)) {
        set.new_parameter_constraint(event, final_guard_function <= handler_zero);
        return {};
    }

    // Split the set into three components,
    //   hitting_set: points which hit the guard in the future
    //   missing_set: points whose trajectories miss the guard completely
    //   past_set: points which would have hit the guard in the past, but do not in the future
    HybridEnclosure& missing_set = set;
    HybridEnclosure hitting_set = set;
    HybridEnclosure past_set = set;
    missing_set.new_parameter_constraint(event, maximal_guard_function <= handler_zero);
    hitting_set.new_parameter_constraint(event, maximal_guard_function >= handler_zero);
    hitting_set.new_parameter_constraint(event, elapsed_time_function <= critical_time_function);
    hitting_set.new_parameter_constraint(event, final_guard_function <= handler_zero);
    past_set.new_parameter_constraint(event, maximal_guard_function >= handler_zero);
    past_set.new_parameter_constraint(event, critical_time_function <= handler_zero);
    hitting_set.reduce();
    missing_set.reduce();
    past_set.reduce();
    LOGGING_PRINTLN_AT(1, "final_guard_range=" << final_guard_function.range());
    LOGGING_PRINTLN_AT(1, "maximal_guard_function=" << maximal_guard_function);
    LOGGING_PRINTLN_AT(1, "final_guard_function=" << final_guard_function);
    LOGGING_PRINTLN_AT(1, "elapsed_time_function=" << elapsed_time_function);
    LOGGING_PRINTLN_AT(1, "critical_time_function=" << critical_time_function);
    LOGGING_PRINTLN_AT(1, "missing_set.is_empty()=" << missing_set.is_empty() << ", hitting_set.is_empty()=" << hitting_set.is_empty() << ", past_set.is_empty()=" << past_set.is_empty());
    LOGGING_PRINTLN_AT(1, "hitting_set=" << hitting_set);
    // FIXME: The guard range at the critical time may provide a stronger guarantee of emptiness than maximal_guard_function > zero, but we should ideally not rely on this
    if (definitely(guard_range_at_critical_time > 0 || missing_set.is_empty())) {
        // swap out missing set
        std::swap(hitting_set, missing_set);
    } else if (not definitely(hitting_set.is_empty())) {
        additional.append(hitting_set);
    }
    if (not definitely(past_set.is_empty())) {
        additional.append(past_set);
    }

    return additional;
}

FlowHandler select_flow_handler(const HybridFlowGuardCapability& capability) {
    switch (capability.method) {
        case HybridFlowGuardMethod::KEEP:
            return keep_flow;
        case HybridFlowGuardMethod::FIRST_ENTRY_BOUND:
            return restrict_crossing_time;
        case HybridFlowGuardMethod::GUARD_SIDE:
            return restrict_guard_side;
        case HybridFlowGuardMethod::TWO_SAMPLES:
            return restrict_sampled_flow;
        case HybridFlowGuardMethod::MAXIMUM_SPLIT:
            return restrict_maximum_flow;
        case HybridFlowGuardMethod::DISCARD:
            return discard_flow;
        default:
            ARIADNE_FAIL_MSG("Unhandled hybrid flow-guard method");
    }
}

HybridJumpCandidate initial_candidate(const HybridGuardContext& context) {
    return {context.initial, context.factory.create_zero(context.initial.parameter_domain())};
}

HybridJumpCandidate tube_candidate(const HybridGuardContext& context, const TimingData& timing) {
    auto result = initial_candidate(context);
    result.enclosure.apply_parameter_reach_step(context.flow, timing.parameter_dependent_evolution_time);
    result.elapsed_time = embed(context.initial.parameter_domain(), timing.evolution_time_coordinate);
    return result;
}

using JumpHandler = HybridJumpCandidate (*)(const HybridGuardContext&,
                                            const EffectiveVectorMultivariateFunction&, const TimingData&);

HybridJumpCandidate activate_source(const HybridGuardContext& context,
                                    const EffectiveVectorMultivariateFunction&, const TimingData&) {
    auto result = initial_candidate(context);
    restrict_source_guard(result.enclosure, context.transition.event, context.transition.guard_function,
                          context.transition.source_arguments);
    return result;
}

HybridJumpCandidate activate_tube(const HybridGuardContext& context,
                                  const EffectiveVectorMultivariateFunction&, const TimingData& timing) {
    auto result = tube_candidate(context, timing);
    restrict_source_guard(result.enclosure, context.transition.event, context.transition.guard_function,
                          context.transition.source_arguments);
    return result;
}

HybridJumpCandidate first_entry_at_time(const HybridGuardContext& context,
                                        const EffectiveVectorMultivariateFunction&, const TimingData& timing) {
    auto result = initial_candidate(context);
    result.elapsed_time = unchecked_compose(context.crossing.crossing_time()->time(), context.initial.state_function());
    if (timing.step_kind != StepKind::CONSTANT_EVOLUTION_TIME &&
        possibly((result.elapsed_time - timing.parameter_dependent_evolution_time).range().upper_bound() > handler_zero)) {
        result.enclosure.new_parameter_constraint(guard_step_event,
                                                  result.elapsed_time <= timing.parameter_dependent_evolution_time);
    }
    result.enclosure.apply_parameter_evolve_step(context.flow, result.elapsed_time);
    return result;
}

HybridJumpCandidate boundary_entry(const HybridGuardContext& context,
                                   const EffectiveVectorMultivariateFunction&, const TimingData& timing) {
    auto result = tube_candidate(context, timing);
    result.enclosure.new_guard(context.transition.event, context.transition.guard_function);
    return result;
}

HybridJumpCandidate directed_entry(const HybridGuardContext& context,
                                   const EffectiveVectorMultivariateFunction& dynamic, const TimingData& timing) {
    auto result = boundary_entry(context, dynamic, timing);
    result.enclosure.new_invariant(context.transition.event,
                                   -lie_derivative(context.transition.guard_function, dynamic));
    return result;
}

HybridJumpCandidate discard_jump(const HybridGuardContext& context,
                                 const EffectiveVectorMultivariateFunction&, const TimingData&) {
    auto result = initial_candidate(context);
    result.enclosure.new_invariant(context.transition.event, context.transition.guard_function * 0 + 1);
    return result;
}

JumpHandler select_entry_handler(const HybridJumpGuardCapability& capability) {
    switch (capability.method) {
        case HybridJumpGuardMethod::SOURCE_ACTIVATION:
            return activate_source;
        case HybridJumpGuardMethod::TUBE_ACTIVATION:
            return activate_tube;
        case HybridJumpGuardMethod::FIRST_ENTRY:
            return first_entry_at_time;
        case HybridJumpGuardMethod::BOUNDARY:
            return boundary_entry;
        case HybridJumpGuardMethod::DIRECTED_BOUNDARY:
            return directed_entry;
        case HybridJumpGuardMethod::DISCARD:
            return discard_jump;
        default:
            ARIADNE_FAIL_MSG("Unhandled hybrid jump-guard method");
    }
}

void apply_impact_direction(HybridEnclosure& set, const HybridGuardContext& context,
                            const EffectiveVectorMultivariateFunction& dynamic) {
    const auto derivative = lie_derivative(context.transition.guard_function, dynamic);
    const auto derivative_range = apply(derivative, set.state_bounding_box());
    if (possibly(derivative_range.lower_bound() <= 0)) {
        set.new_activation(context.transition.event, derivative);
    }
    // Preserve the original contraction before other blocking guards inspect
    // the impact candidate, even when its emptiness result is not yet used.
    (void)definitely(set.is_empty());
}
} // namespace

HybridTransitionPolicy transition_policy(EventKind event) {
    switch (event) {
        case EventKind::INVARIANT:
            return {true, false, false, false, false};
        case EventKind::PROGRESS:
            return {true, false, true, false, false};
        case EventKind::PERMISSIVE:
            return {false, true, true, true, false};
        case EventKind::URGENT:
            return {true, true, true, true, false};
        case EventKind::IMPACT:
            return {true, true, false, false, true};
        default:
            ARIADNE_FAIL_MSG("Unrecognised hybrid event kind: " << event);
    }
}

HybridFlowGuardCapability flow_guard_capability(const CrossingData& crossing, Semantics semantics) {
    validate_semantics(semantics);
    if (crossing.crossing_time() != nullptr)
        return {HybridFlowGuardMethod::FIRST_ENTRY_BOUND, CrossingEvidenceKind::CROSSING_TIME, true};
    if (crossing.guard_maximum() != nullptr)
        return {HybridFlowGuardMethod::MAXIMUM_SPLIT, CrossingEvidenceKind::GUARD_MAXIMUM, true};
    switch (crossing.geometry()) {
        case CrossingGeometry::INCREASING:
        case CrossingGeometry::CONVEX:
            return {HybridFlowGuardMethod::GUARD_SIDE, CrossingEvidenceKind::NONE, false};
        case CrossingGeometry::CONCAVE:
        case CrossingGeometry::UNRESOLVED:
            return semantics == Semantics::UPPER
                     ? HybridFlowGuardCapability{HybridFlowGuardMethod::TWO_SAMPLES, CrossingEvidenceKind::NONE, true}
                     : HybridFlowGuardCapability{HybridFlowGuardMethod::DISCARD, CrossingEvidenceKind::NONE, false};
        case CrossingGeometry::NEGATIVE:
        case CrossingGeometry::POSITIVE:
        case CrossingGeometry::DECREASING:
            return {HybridFlowGuardMethod::KEEP, CrossingEvidenceKind::NONE, false};
        default:
            ARIADNE_FAIL_MSG("Unhandled crossing geometry in flow restriction");
    }
}

HybridJumpGuardCapability jump_guard_capability(EventKind event, const CrossingData& crossing, Semantics semantics) {
    validate_semantics(semantics);
    const auto policy = transition_policy(event);
    if (!policy.activates_reset)
        ARIADNE_FAIL_MSG("Invalid event kind for jump construction: " << event);
    if (!policy.blocks_flow)
        return {HybridJumpGuardMethod::TUBE_ACTIVATION, CrossingEvidenceKind::NONE, true, false};

    const bool impact = policy.requires_impact_direction;
    if (crossing.crossing_time() != nullptr)
        return {HybridJumpGuardMethod::FIRST_ENTRY, CrossingEvidenceKind::CROSSING_TIME, true, impact};
    if (crossing.guard_maximum() != nullptr)
        return {HybridJumpGuardMethod::DIRECTED_BOUNDARY, CrossingEvidenceKind::GUARD_MAXIMUM, true, true};
    switch (crossing.geometry()) {
        case CrossingGeometry::INCREASING:
        case CrossingGeometry::CONVEX:
            return {HybridJumpGuardMethod::BOUNDARY, CrossingEvidenceKind::NONE, true, impact};
        case CrossingGeometry::CONCAVE:
            return {HybridJumpGuardMethod::DIRECTED_BOUNDARY, CrossingEvidenceKind::NONE, true, true};
        case CrossingGeometry::UNRESOLVED:
            return semantics == Semantics::UPPER
                     ? HybridJumpGuardCapability{HybridJumpGuardMethod::DIRECTED_BOUNDARY, CrossingEvidenceKind::NONE, true, true}
                     : HybridJumpGuardCapability{HybridJumpGuardMethod::DISCARD, CrossingEvidenceKind::NONE, false, impact};
        case CrossingGeometry::NEGATIVE:
        case CrossingGeometry::POSITIVE:
        case CrossingGeometry::DECREASING:
            return {HybridJumpGuardMethod::SOURCE_ACTIVATION, CrossingEvidenceKind::NONE, false, impact};
        default:
            ARIADNE_FAIL_MSG("Unhandled crossing geometry in jump construction");
    }
}

HybridInitialGuardActivation classify_initial_guard(
    const HybridEnclosure& invariant_valid_source, const TransitionData& transition) {
    switch (classify_source_guard(invariant_valid_source, transition.guard_function,
                                  transition.source_arguments)) {
        case HybridSourceGuardStatus::EXCLUDED:
            return HybridInitialGuardActivation::EXCLUDED;
        case HybridSourceGuardStatus::UNRESOLVED:
            return HybridInitialGuardActivation::UNRESOLVED;
        case HybridSourceGuardStatus::FULL_SOURCE:
            return HybridInitialGuardActivation::FULL_SOURCE;
        default:
            ARIADNE_FAIL_MSG("Unhandled hybrid source guard status");
    }
}

bool permits_initial_reset(HybridInitialGuardActivation activation, Semantics semantics) {
    switch (semantics) {
        case Semantics::UPPER:
            return activation != HybridInitialGuardActivation::EXCLUDED;
        case Semantics::LOWER:
            return activation == HybridInitialGuardActivation::FULL_SOURCE;
        default:
            ARIADNE_FAIL_MSG("Unhandled semantics in initial guard activation");
    }
}

bool defer_endpoint_guard(const CrossingData& crossing) {
    return crossing.crossing_time() != nullptr || crossing.geometry() == CrossingGeometry::INCREASING;
}

HybridJumpCandidate make_jump_candidate(const HybridGuardContext& context,
                                        const EffectiveVectorMultivariateFunction& dynamic, const TimingData& timing) {
    LOGGING_SCOPE_CREATE;
    const auto capability = jump_guard_capability(context.transition.event_kind, context.crossing, context.semantics);
    LOGGING_PRINTLN_AT(1, "jump event=" << context.transition.event << " kind=" << context.transition.event_kind
                                        << " semantics=" << (context.semantics == Semantics::UPPER ? "UPPER" : "LOWER"));
    const JumpHandler handler = select_entry_handler(capability);
    auto result = handler(context, dynamic, timing);
    if (transition_policy(context.transition.event_kind).requires_impact_direction)
        apply_impact_direction(result.enclosure, context, dynamic);
    ARIADNE_PRECONDITION(result.elapsed_time.domain() == result.enclosure.parameter_domain());
    return result;
}

Void apply_flow_guard(List<HybridEnclosure>& sets, const HybridGuardContext& context,
                      const ValidatedScalarMultivariateFunctionPatch& elapsed) {
    LOGGING_SCOPE_CREATE;
    validate_semantics(context.semantics);
    const FlowHandler handler = select_flow_handler(flow_guard_capability(context.crossing, context.semantics));
    LOGGING_PRINTLN_AT(1, "flow restriction event=" << context.transition.event
                                                    << " crossing=" << context.crossing.summary()
                                                    << " semantics=" << (context.semantics == Semantics::UPPER ? "UPPER" : "LOWER"));
    auto initial_state = context.initial.state_function();
    if (elapsed.domain().dimension() > context.initial.parameter_domain().dimension()) {
        const auto time_domain = elapsed.domain()[elapsed.argument_size() - 1u];
        initial_state = embed(initial_state, time_domain);
    }
    const FlowContext flow_context{context, elapsed, initial_state};
    const SizeType original_size = sets.size();
    for (SizeType i = 0u; i != original_size; ++i) {
        ARIADNE_PRECONDITION(sets[i].parameter_domain() == elapsed.domain());
        // Numerical handlers finish using the primary reference before any
        // append can invalidate it. Preserve primaries first, then extras.
        const auto additional = handler(sets[i], flow_context);
        for (const auto& extra: additional)
            sets.append(extra);
    }
}
} // namespace Ariadne::EvolutionDetail
