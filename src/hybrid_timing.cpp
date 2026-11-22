/***************************************************************************
 *            hybrid_timing.cpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#include "hybrid_timing.hpp"

#include <iomanip>

#include "algebra/vector.hpp"
#include "function/function.hpp"
#include "function/function_model.hpp"
#include "function/function_patch.hpp"
#include "logging/logging.hpp"
#include "numeric/numeric.hpp"
#include "system/hybrid_automaton_interface.hpp"

#include "equation-solving/equation_solver.hpp"

namespace Ariadne {
template<class PR>
inline auto operator*(FloatUpperBound<PR> x1, PositiveFloat<PR> x2) -> FloatUpperBound<PR> {
    return FloatUpperBound<PR>(mul(up, x1._u, x2._v));
}
template<class PR>
inline auto operator*(PositiveFloatUpperBound<PR> x1, PositiveFloat<PR> x2) -> PositiveFloatUpperBound<PR> {
    return PositiveFloatUpperBound<PR>(mul(up, x1._u, x2._v));
}

inline auto operator+(Int n, FloatDPBounds x) -> decltype(FloatDP(n, dp) + x) { return FloatDP(n, dp) + x; }
inline auto operator>(FloatDPLowerBound x, Real r) -> decltype(x > declval<FloatDPBounds>()) { return x > FloatDPBounds(r, dp); }
inline auto operator>=(FloatDPLowerBound x, Real r) -> decltype(x >= declval<FloatDPBounds>()) { return x >= FloatDPBounds(r, dp); }
inline auto operator>(FloatDPUpperBound x, Real r) -> decltype(x > declval<FloatDPBounds>()) { return x > FloatDPBounds(r, dp); }
inline auto operator<=(FloatDPUpperBound x, Real r) -> decltype(x <= declval<FloatDPBounds>()) { return x <= FloatDPBounds(r, dp); }

namespace EvolutionDetail {
namespace {
const FloatDP timing_zero = {0, dp};
using ExactIntervalVectorType = Vector<ExactIntervalType>;

/** @brief Choose the final-time strategy before consulting guard crossings. */
ValidatedScalarMultivariateFunctionPatch select_temporal_evolution(
    TimingData& result, const FloatDP& step_size, const FloatDPBounds& final_time_bounds,
    const UpperIntervalType& remaining_time_range,
    const ValidatedScalarMultivariateFunctionPatch& time_identity,
    const ValidatedFunctionPatchFactory& patch_factory, Bool allow_creep) {
    LOGGING_SCOPE_CREATE;
    // The time-dependent part of the evolution time
    ValidatedScalarMultivariateFunctionPatch temporal_evolution_time = patch_factory.create_zero(ExactIntervalVectorType(1u, time_domain));

    LOGGING_PRINTLN(std::fixed << "remaining_time_range=" << remaining_time_range);
    if (possibly(remaining_time_range.lower_bound() < timing_zero)) {
        // Some of the points may already have reached the final time.
        // Don't try anything fancy, just do a simple constant time step.
        if (definitely(remaining_time_range.upper_bound() <= step_size)) {
            result.step_kind = StepKind::CONSTANT_FINISHING_TIME;
            result.finishing_kind = FinishingKind::AT_FINAL_TIME;
            temporal_evolution_time = final_time_bounds - time_identity;
        } else {
            result.step_kind = StepKind::CONSTANT_EVOLUTION_TIME;
            result.finishing_kind = FinishingKind::STRADDLE_FINAL_TIME;
            temporal_evolution_time = FloatDPBounds(step_size);
        }
    } else if (definitely(remaining_time_range.upper_bound() <= result.step_size)) {
        // The rest of the evolution can be computed within a single time step.
        // The finishing kind is given as AT_FINAL_TIME so that the evolution algorithm
        // knows that the evolved set does not need to be evolved further.
        // This knowledge is required to be given combinarially, since
        // specifying the final time as a constant Function is not
        // exact if the final_time parameter is not exactly representable as
        // a FloatDP
        result.step_kind = StepKind::CONSTANT_FINISHING_TIME;
        result.finishing_kind = FinishingKind::AT_FINAL_TIME;
        temporal_evolution_time = final_time_bounds - time_identity;
    } else if (possibly(remaining_time_range.lower_bound() <= step_size) && allow_creep) {
        // Some of the evolved points can be evolved to the final time in a single step
        // The evolution is performed over a step size which moves points closer to the final time, but does not cross.

        // Using the final_time as a guide, set the finishing time to closer to the final time.
        // This method ensures that points do not pass the final time after the transition.
        result.step_kind = StepKind::SPACETIME_DEPENDENT_FINISHING_TIME;
        result.finishing_kind = FinishingKind::BEFORE_FINAL_TIME;
        PositiveFloatDP sf = {1u, dp};
        while (possibly(remaining_time_range.upper_bound() * sf > step_size)) {
            sf = hlf(sf);
        }
        temporal_evolution_time = FloatDPBounds(sf) * (final_time_bounds - time_identity);
    } else { // remaining_time_range.lower_bound()>step_size)
        // As far as timing goes, perform the evolution over a full time step
        result.step_kind = StepKind::CONSTANT_EVOLUTION_TIME;
        result.finishing_kind = FinishingKind::BEFORE_FINAL_TIME;
        temporal_evolution_time = result.step_size;
    }

    LOGGING_PRINTLN_AT(1, "finishing_kind=" << result.finishing_kind);
    LOGGING_PRINTLN_AT(1, "temporal_evolution_time=" << temporal_evolution_time);

    return temporal_evolution_time;
}

/** @brief Test whether the full-step endpoint admits a partial guard crossing. */
Bool needs_guard_creep(
    const FlowFunctionModel& flow, const HybridEnclosure& initial_set,
    const Map<DiscreteEvent, CrossingData>& crossings,
    const Map<DiscreteEvent, TransitionData>& transitions,
    const ValidatedFunctionPatchFactory& patch_factory,
    const FloatDP& step_size, const TimingData& result, Bool allow_creep) {
    LOGGING_SCOPE_CREATE;
    Bool creep = false;

    // Test for creep step
    if (allow_creep && !crossings.empty() && (result.finishing_kind == FinishingKind::BEFORE_FINAL_TIME || result.finishing_kind == FinishingKind::STRADDLE_FINAL_TIME)) {
        // If an event is possible, but only some points reach the guard set
        // after a full time step, then terminating the evolution here will
        // cause a splitting of the enclosure set. To prevent this, attempt
        // to "creep" up to the event guard boundary, so that in the next step,
        // all points can be made to cross.
        // Test to see if a creep step is needed
        LOGGING_PRINTLN("Possible creep step; h=" << step_size);

        HybridEnclosure evolve_set = initial_set;
        evolve_set.apply_fixed_evolve_step(flow, flow.step_size());

        for (Map<DiscreteEvent, CrossingData>::Iterator crossing_iter = crossings.begin();
             crossing_iter != crossings.end(); ++crossing_iter) {
            DiscreteEvent event = crossing_iter->first;
            EventKind event_kind = transitions[event].event_kind;
            const auto geometry = crossing_iter->second.geometry();
            LOGGING_PRINTLN_AT(1, "Event " << event << ": " << event_kind << ": " << geometry);
            if (event_kind != EventKind::PERMISSIVE) {
                evolve_set.new_invariant(event, transitions[event].guard_function);
            }
            // FIXME: When using permissive crossings, jumps in the step after the crossing time are lost.
            // A hack to fix this is to only creep on non-permissive events. Check that evolution is correct in this case.
            // FIXME: What should we do on increasing but non-transverse crossings?
            if ((crossing_iter->second.crossing_time() != nullptr) && event_kind != EventKind::PERMISSIVE) {
                LOGGING_PRINTLN_AT(1, "crossing_time_range=" << crossing_iter->second.crossing_time()->time().range());
                const ValidatedScalarMultivariateFunctionPatch& crossing_time = crossing_iter->second.crossing_time()->time();
                UpperIntervalType crossing_time_range = crossing_time.range();
                // Every supported non-permissive event blocks evolution;
                // the enclosing condition has already excluded PERMISSIVE.
                if (definitely(crossing_time_range.upper_bound() < step_size)) {
                    // NOTE: Use strict comparison here so that guard is fully crossed
                    // In principle this is not necessary, but this would involve testing
                    // to ensure that the evolved set is not propagated

                    // This event ensures that the evolve set is empty after a full step, so use this.
                    LOGGING_PRINTLN_AT(1, std::setprecision(18) << "crossing_time_range=" << crossing_time_range << ", crossing_time_range .upper_bound()=" << crossing_time_range.upper_bound() << ", step_size=" << step_size);
                    EffectiveScalarMultivariateFunction guard = transitions[event].guard_function;
                    ValidatedVectorMultivariateFunctionPatch identity = patch_factory.create_identity(crossing_time.domain());
                    // TODO: Remove use of cast
                    ValidatedScalarMultivariateFunctionPatch step_time = crossing_time * 0 + static_cast<Dyadic>(step_size);
                    LOGGING_PRINTLN_AT(1, "full flow=" << compose(flow, join(identity, step_time)));
                    LOGGING_PRINTLN_AT(1, "guard range at crossing time=" << compose(guard, compose(flow, join(initial_set.state_function(), compose(crossing_time, initial_set.state_function())))).range());
                    LOGGING_PRINTLN_AT(1, "guard range at crossing time=" << compose(guard, compose(flow, join(identity, crossing_time))).range());
                    LOGGING_PRINTLN_AT(1, "No creep; event " << crossing_iter->first << " completely taken");
                    creep = false;
                    break;
                } else if (possibly(crossing_time_range.lower_bound() <= timing_zero)) {
                    // This event is already partially active, so carry on with a full step
                    LOGGING_PRINTLN_AT(1, "No creep; event " << crossing_iter->first << " already partially active");
                    creep = false;
                    break;
                } else if (definitely(crossing_time_range.lower_bound() >= step_size)) {
                    LOGGING_PRINTLN_AT(1, "Event " << crossing_iter->first << " is not actually reached");
                    // crossings.erase(pending_erase_crossing_iter);
                } else {
                    LOGGING_PRINTLN_AT(1, "Event " << crossing_iter->first << " can be crept up to.");
                    creep = true;
                }
            }
        }
        // If the evolved set is definitely empty, no creeping occurs
        if (definitely(evolve_set.is_empty())) {
            LOGGING_PRINTLN_AT(1, "No creep; evolve set is empty");
            creep = false;
        }
    }

    if (allow_creep && creep == false) {
        LOGGING_PRINTLN("No creep");
    }

    return creep;
}

/** @brief Construct the creep time using the configured implicit solver.
 *  Preserve guard iteration order and erase only successfully solved crossings.
 */
Void apply_guard_creep(
    Bool creep, TimingData& result,
    ValidatedScalarMultivariateFunctionPatch& spacial_evolution_time,
    Map<DiscreteEvent, CrossingData>& crossings,
    const Map<DiscreteEvent, TransitionData>& transitions,
    const FlowFunctionModel& flow, const HybridEnclosure& initial_set,
    const HybridAutomatonInterface& system,
    const EquationSolverInterface& solver,
    const ValidatedFunctionPatchFactory& patch_factory,
    const FloatDP& step_size, SizeType n) {
    LOGGING_SCOPE_CREATE;
    if (creep) {
        // If an event is possible, but only some points reach the guard set
        // after a full time step, then terminating the evolution here will
        // cause a splitting of the enclosure set. To prevent this, attempt
        // to "creep" up to the event guard boundary, so that in the next step,
        // all points can be made to cross.
        // Test to see if a creep step is needed
        LOGGING_PRINTLN("Possible creep step; h=" << step_size);
        LOGGING_PRINTLN("crossings=" << crossings);

        // The configured solver is passed in; its fallback remains unchanged.

        HybridEnclosure evolve_set = initial_set;
        evolve_set.apply_fixed_evolve_step(flow, flow.step_size());
        EffectiveVectorMultivariateFunction dynamic = system.dynamic_function(initial_set.location());
        ExactBoxType flow_spacial_domain = project(flow.domain(), range(0, flow.argument_size() - 1u));
        ExactIntervalType flow_time_domain = flow.domain()[flow.argument_size() - 1u];
        ValidatedScalarMultivariateFunctionPatch zero_function = factory(flow).create_zero();
        ValidatedVectorMultivariateFunctionPatch identity_function = factory(flow).create_identity();
        // TODO: Remove use of cast
        ValidatedVectorMultivariateFunctionPatch space_projection = flow * static_cast<Dyadic>(timing_zero);
        for (SizeType i = 0; i != n; ++i) {
            space_projection[i] = space_projection[i] + identity_function[i];
        }

        // static const ExactDouble CREEP_MAXIMUM=1.0_X;
        static const ExactDouble CREEP_MAXIMUM = 0.9375_x;
        spacial_evolution_time = patch_factory.create_constant(flow.space_domain(), flow.step_size() * CREEP_MAXIMUM);

        for (Map<DiscreteEvent, CrossingData>::Iterator crossing_iter = crossings.begin();
             crossing_iter != crossings.end();) {
            DiscreteEvent event = crossing_iter->first;
            EventKind event_kind = transitions[event].event_kind;
            const auto geometry = crossing_iter->second.geometry();
            EffectiveScalarMultivariateFunction guard_function = transitions[event].guard_function;
            LOGGING_PRINTLN_AT(1, "Event " << event << ": " << event_kind << ": " << geometry);
            if (event_kind != EventKind::PERMISSIVE) {
                UpperIntervalType guard_range = compose(guard_function, flow).range();
                ARIADNE_ASSERT(decide(guard_range.lower_bound() < timing_zero));
                UpperIntervalType guard_derivative_range = compose(lie_derivative(guard_function, dynamic), flow).range();

                FloatDPBounds alpha_val = (1 + flow.step_size() * cast_exact(guard_derivative_range.lower_bound()) / cast_exact(guard_range.lower_bound()));
                FloatDP alpha = alpha_val.value();
                LOGGING_PRINTLN_AT(1, "step_size: " << flow.step_size() << ", guard_range: " << guard_range << ", guard_derivative_range: " << guard_derivative_range << ", alpha: " << alpha);
                if (alpha > 0 && alpha <= 1) {
                    ValidatedScalarMultivariateFunctionPatch guard_creep_time;
                    Bool successfully_computed_guard_creep_time = false;
                    try {
                        guard_creep_time = solver.implicit(compose(guard_function, flow) - static_cast<Dyadic>(alpha) * compose(guard_function, space_projection),
                                                           flow_spacial_domain, flow_time_domain);
                        LOGGING_PRINTLN_AT(1, "guard_creep_time= " << guard_creep_time);
                        LOGGING_PRINTLN_AT(1, "guard_creep_time.range()=" << guard_creep_time.range());
                        successfully_computed_guard_creep_time = true;
                        LOGGING_PRINTLN_AT(2, "sucessfully_computed_guard_creep_time=" << successfully_computed_guard_creep_time);
                    } catch (...) {
                        LOGGING_PRINTLN_AT(1, "Error in computing guard creep time");
                    }
                    if (successfully_computed_guard_creep_time) {
                        spacial_evolution_time = spacial_evolution_time * (guard_creep_time / flow.step_size());
                        LOGGING_PRINTLN_AT(2, "spacial_evolution_time=" << spacial_evolution_time);
                        LOGGING_PRINTLN_AT(2, "crossings before erasing=" << crossings);
                        crossings.erase(crossing_iter++);
                        LOGGING_PRINTLN_AT(2, "crossings after erasing=" << crossings);
                    } else {
                        ++crossing_iter;
                    }
                } else {
                    ++crossing_iter;
                }
            } else {
                ++crossing_iter;
            }
        }
        spacial_evolution_time.clobber();

        LOGGING_PRINTLN_AT(1, "Creep step: spacial_evolution_time=" << spacial_evolution_time);
        LOGGING_PRINTLN_AT(1, "spacial_evolution_time.range()=" << spacial_evolution_time.range());
        LOGGING_PRINTLN_AT(1, "remaining crossings=" << crossings);
        result.step_kind = StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME;
        result.finishing_kind = FinishingKind::BEFORE_FINAL_TIME;
    }

}

/** @brief Reduce a spread in clock values only in the absence of crossings. */
Void apply_clock_unwinding(
    TimingData& result, const Map<DiscreteEvent, CrossingData>& crossings,
    const HybridEnclosure& initial_set,
    const ValidatedScalarMultivariateFunctionPatch& starting_time_function,
    const UpperIntervalType& starting_time_range, const FloatDP& step_size,
    const ValidatedFunctionPatchFactory& patch_factory, Bool allow_unwind) {
    LOGGING_SCOPE_CREATE;
    // Test to see if it is possible to unwind crossings
    if (allow_unwind && crossings.empty() && (starting_time_range.lower_bound().raw() < starting_time_range.upper_bound().raw()) && result.step_kind == StepKind::CONSTANT_EVOLUTION_TIME && result.finishing_kind == FinishingKind::BEFORE_FINAL_TIME) {
        LOGGING_PRINTLN_AT(1, "Possible unwind step; starting_time_range=" << starting_time_range << ", step_size=" << step_size);
        // Try to unwind the evolution time to a constant
        result.step_kind = StepKind::PARAMETER_DEPENDENT_FINISHING_TIME;
        result.finishing_kind = FinishingKind::BEFORE_FINAL_TIME;
        if (decide(starting_time_range.width() * 2 < step_size)) {
            result.parameter_dependent_finishing_time = patch_factory.create_constant(initial_set.parameter_domain(), cast_exact(starting_time_range.lower_bound()) + step_size);
        } else {
            // Try to reduce the time interval by half the step size
            // Corresponds to setting omega(smin)=tau(smin)+h, omega(smax)=tau(smax)+h/2
            // Taking omega(s)=a tau(s) + b, we obtain
            //   a=1-h/2(tmax-tmin);  b=h(tmax-tmin/2)/(tmax-tmin) = (2tmax-tmin)a
            FloatDP h = {result.step_size, dp};
            FloatDP tmin = cast_exact(starting_time_range.lower_bound());
            FloatDP tmax = cast_exact(starting_time_range.upper_bound());
            FloatDPBounds a = 1 - (hlf(h) / (tmax - tmin));
            FloatDPBounds b = h * (tmax - hlf(tmin)) / (tmax - tmin);
            result.parameter_dependent_finishing_time = a * starting_time_function + b;
        }
        LOGGING_PRINTLN_AT(1, "Unwinding to time " << result.parameter_dependent_finishing_time);
        result.parameter_dependent_evolution_time = result.parameter_dependent_finishing_time - starting_time_function;
    }

}

} // namespace

TimingData estimate_hybrid_timing(
    Real final_time, const FlowFunctionModel& flow,
    Map<DiscreteEvent, CrossingData>& crossings,
    const Map<DiscreteEvent, TransitionData>& transitions,
    const HybridEnclosure& initial_set, const HybridAutomatonInterface& system,
    const ValidatedFunctionPatchFactory& patch_factory,
    const EquationSolverInterface& solver, Bool allow_creep, Bool allow_unwind) {
    // Compute the evolution time for the given step.
    LOGGING_SCOPE_CREATE;

    const SizeType n = flow.result_size();
    const FloatDP step_size = flow.domain()[flow.domain().size() - 1].upper_bound();
    const FloatDPBounds final_time_bounds(final_time, dp);

    TimingData result;
    result.step_size = flow.step_size();
    result.final_time = final_time;

    ExactBoxType state_domain = cast_exact_box(initial_set.state_bounding_box());
    ExactIntervalType time_domain = cast_exact_interval(initial_set.time_range() + ExactIntervalType(timing_zero, step_size));
    ExactBoxType statetime_domain = product(state_domain, time_domain);

    // ValidatedVectorMultivariateFunctionPatch space_coordinates=patch_factory.create_identity(space_domain);
    ValidatedScalarMultivariateFunctionPatch time_coordinate = patch_factory.create_coordinate(statetime_domain, n);
    ValidatedScalarMultivariateFunctionPatch time_identity = patch_factory.create_coordinate({time_domain}, 0);

    result.evolution_time_domain = ExactIntervalType(timing_zero, step_size);
    result.evolution_time_coordinate = patch_factory.create_coordinate({result.evolution_time_domain}, 0);

    ExactBoxType flow_state_domain = ExactBoxType(project(flow.domain(), range(0, n)));
    if (!subset(state_domain, flow_state_domain)) {
        ARIADNE_WARN(std::setprecision(17) << "Bounding box " << state_domain << " is not subset of the flow spacial domain " << flow_state_domain);
        state_domain = hull(state_domain, flow_state_domain);
    }

    // NOTE: The starting time function may be negative or greater than the final time
    // over part of the parameter domain.
    const ValidatedVectorMultivariateFunctionPatch& starting_state_function = initial_set.state_function();
    const ValidatedScalarMultivariateFunctionPatch& starting_time_function = initial_set.time_function();
    UpperIntervalType starting_time_range = initial_set.time_range();
    UpperIntervalType remaining_time_range = final_time_bounds - starting_time_range;

    LOGGING_PRINTLN(std::fixed << "starting_time_range=" << starting_time_range << " step_size=" << step_size << " final_time=" << final_time);

    const auto temporal_evolution_time = select_temporal_evolution(
        result, step_size, final_time_bounds, remaining_time_range,
        time_identity, patch_factory, allow_creep);

    ValidatedScalarMultivariateFunctionPatch spacial_evolution_time = patch_factory.create_constant(state_domain, FloatDP(step_size));

    const Bool creep = needs_guard_creep(
        flow, initial_set, crossings, transitions, patch_factory, step_size, result, allow_creep);
    apply_guard_creep(creep, result, spacial_evolution_time, crossings, transitions,
                      flow, initial_set, system, solver, patch_factory, step_size, n);

    // TODO: Remove use of cast
    ValidatedScalarMultivariateFunctionPatch evolution_time = embed(spacial_evolution_time, time_domain) * embed(state_domain, temporal_evolution_time / static_cast<Dyadic>(step_size));
    ValidatedScalarMultivariateFunctionPatch finishing_time = evolution_time + time_coordinate;

    LOGGING_PRINTLN_AT(1, "evolution_time=" << (evolution_time));
    LOGGING_PRINTLN_AT(1, "finishing_time=" << (finishing_time));
    result.spacetime_dependent_evolution_time = evolution_time;
    result.spacetime_dependent_finishing_time = finishing_time;
    result.parameter_dependent_evolution_time = unchecked_compose(evolution_time, join(starting_state_function, starting_time_function));
    result.parameter_dependent_finishing_time = unchecked_compose(finishing_time, join(starting_state_function, starting_time_function));

    apply_clock_unwinding(result, crossings, initial_set, starting_time_function,
                          starting_time_range, step_size, patch_factory, allow_unwind);

    LOGGING_PRINTLN_AT(1, "step_kind=" << result.step_kind << ", finishing_kind=" << result.finishing_kind);
    LOGGING_PRINTLN_AT(1, "parameter_dependent_evolution_time=" << result.parameter_dependent_evolution_time);
    LOGGING_PRINTLN_AT(1, "parameter_dependent_finishing_time=" << result.parameter_dependent_finishing_time << "\n");
    return result;
}

} // namespace EvolutionDetail
} // namespace Ariadne
