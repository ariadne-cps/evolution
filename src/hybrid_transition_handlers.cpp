/***************************************************************************
 *            hybrid_transition_handlers.cpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#include "hybrid_transition_handlers.hpp"

#include "function/function.hpp"
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
                                 context.elapsed <= unchecked_compose(context.guard.crossing.crossing_time, context.initial_state));
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

ExtraBranches restrict_grazing_flow(HybridEnclosure& set, const FlowContext& context) {
    const auto event = context.guard.transition.event;
    const auto& guard_function = context.guard.transition.guard_function;
    const auto& guard_crossing_data = context.guard.crossing;
    const auto& starting_state_function = context.initial_state;
    const auto& elapsed_time_function = context.elapsed;
    const auto& flow = context.guard.flow;
    ExtraBranches additional;

    ValidatedScalarMultivariateFunctionPatch critical_time_function = unchecked_compose(guard_crossing_data.critical_time, starting_state_function);
    ValidatedScalarMultivariateFunctionPatch final_guard_function = compose(guard_function, unchecked_compose(flow, join(starting_state_function, elapsed_time_function)));
    ValidatedScalarMultivariateFunctionPatch maximal_guard_function = compose(guard_function, unchecked_compose(flow, join(starting_state_function, critical_time_function)));
    UpperIntervalType guard_range_at_critical_time = guard_crossing_data.guard_range_at_critical_time;
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

FlowHandler select_flow_handler(CrossingKind crossing, Semantics semantics) {
    switch (crossing) {
        case CrossingKind::TRANSVERSE:
            return restrict_crossing_time;
        case CrossingKind::INCREASING:
        case CrossingKind::CONVEX:
            return restrict_guard_side;
        case CrossingKind::GRAZING:
            return restrict_grazing_flow;
        case CrossingKind::CONCAVE:
        case CrossingKind::DEGENERATE:
            return semantics == Semantics::UPPER ? restrict_sampled_flow : discard_flow;
        case CrossingKind::NEGATIVE:
        case CrossingKind::POSITIVE:
        case CrossingKind::DECREASING:
            return keep_flow;
        default:
            ARIADNE_FAIL_MSG("Unhandled crossing kind in flow restriction: " << crossing);
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
    result.enclosure.new_activation(context.transition.event, context.transition.guard_function);
    return result;
}

HybridJumpCandidate activate_tube(const HybridGuardContext& context,
                                  const EffectiveVectorMultivariateFunction&, const TimingData& timing) {
    auto result = tube_candidate(context, timing);
    result.enclosure.new_activation(context.transition.event, context.transition.guard_function);
    return result;
}

HybridJumpCandidate transverse_entry(const HybridGuardContext& context,
                                     const EffectiveVectorMultivariateFunction&, const TimingData& timing) {
    auto result = initial_candidate(context);
    result.elapsed_time = unchecked_compose(context.crossing.crossing_time, context.initial.state_function());
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

JumpHandler select_entry_handler(CrossingKind crossing, Semantics semantics) {
    switch (crossing) {
        case CrossingKind::TRANSVERSE:
            return transverse_entry;
        case CrossingKind::INCREASING:
        case CrossingKind::CONVEX:
            return boundary_entry;
        case CrossingKind::GRAZING:
        case CrossingKind::CONCAVE:
            return directed_entry;
        case CrossingKind::DEGENERATE:
            return semantics == Semantics::UPPER ? directed_entry : discard_jump;
        case CrossingKind::NEGATIVE:
        case CrossingKind::POSITIVE:
        case CrossingKind::DECREASING:
            return activate_source;
        default:
            ARIADNE_FAIL_MSG("Unhandled crossing kind in jump construction: " << crossing);
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

bool defer_endpoint_guard(CrossingKind crossing) {
    switch (crossing) {
        case CrossingKind::TRANSVERSE:
        case CrossingKind::INCREASING:
            return true;
        case CrossingKind::NEGATIVE:
        case CrossingKind::POSITIVE:
        case CrossingKind::DECREASING:
        case CrossingKind::CONVEX:
        case CrossingKind::GRAZING:
        case CrossingKind::CONCAVE:
        case CrossingKind::DEGENERATE:
            return false;
        default:
            ARIADNE_FAIL_MSG("Unhandled endpoint crossing kind: " << crossing);
    }
}

HybridJumpCandidate make_jump_candidate(const HybridGuardContext& context,
                                        const EffectiveVectorMultivariateFunction& dynamic, const TimingData& timing) {
    LOGGING_SCOPE_CREATE;
    validate_semantics(context.semantics);
    const auto policy = transition_policy(context.transition.event_kind);
    LOGGING_PRINTLN_AT(1, "jump event=" << context.transition.event << " kind=" << context.transition.event_kind
                                        << " semantics=" << (context.semantics == Semantics::UPPER ? "UPPER" : "LOWER"));
    if (!policy.activates_reset) {
        ARIADNE_FAIL_MSG("Invalid event kind for jump construction: " << context.transition.event_kind);
    }
    const JumpHandler handler = policy.blocks_flow
                                  ? select_entry_handler(context.crossing.crossing_kind, context.semantics)
                                  : activate_tube;
    auto result = handler(context, dynamic, timing);
    if (policy.requires_impact_direction)
        apply_impact_direction(result.enclosure, context, dynamic);
    ARIADNE_PRECONDITION(result.elapsed_time.domain() == result.enclosure.parameter_domain());
    return result;
}

Void apply_flow_guard(List<HybridEnclosure>& sets, const HybridGuardContext& context,
                      const ValidatedScalarMultivariateFunctionPatch& elapsed) {
    LOGGING_SCOPE_CREATE;
    validate_semantics(context.semantics);
    const FlowHandler handler = select_flow_handler(context.crossing.crossing_kind, context.semantics);
    LOGGING_PRINTLN_AT(1, "flow restriction event=" << context.transition.event
                                                    << " crossing=" << context.crossing.crossing_kind
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
