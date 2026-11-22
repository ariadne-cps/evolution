/***************************************************************************
 *            hybrid_crossing_analysis.cpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#include "hybrid_crossing_analysis.hpp"

#include "function/function.hpp"
#include "hybrid/hybrid_set.hpp"
#include "logging/logging.hpp"

#include "equation-solving/equation_solver.hpp"

namespace Ariadne::EvolutionDetail {
namespace {
const FloatDP analysis_zero(0, dp);
}

Map<DiscreteEvent, CrossingData> analyze_crossings(
    const Set<DiscreteEvent>& active_events,
    const EffectiveVectorMultivariateFunction& dynamic,
    const Map<DiscreteEvent, EffectiveScalarMultivariateFunction>& guards,
    const ValidatedVectorMultivariateFunctionPatch& flow,
    const HybridEnclosure& initial_set, const EquationSolverInterface& solver) {
    LOGGING_SCOPE_CREATE;

    Map<DiscreteEvent, CrossingData> crossings;

    ExactBoxType flow_spacial_domain = project(flow.domain(), range(0, flow.argument_size() - 1u));
    ExactIntervalType flow_time_domain = flow.domain()[flow.argument_size() - 1u];
    UpperBoxType flow_bounds = flow.range();
    for (Set<DiscreteEvent>::ConstIterator event_iter = active_events.begin();
         event_iter != active_events.end(); ++event_iter) {
        const DiscreteEvent event = *event_iter;
        LOGGING_PRINTLN_AT(1, "event=" << event);
        const EffectiveScalarMultivariateFunction& guard = guards[event];
        LOGGING_PRINTLN_AT(2, "guard=" << guard);

        // Compute the value of the guard function g over the flow of $\dot(x)=f(x)$
        UpperIntervalType guard_range = apply(guard, flow_bounds);
        LOGGING_PRINTLN_AT(2, "guard_range=" << guard_range);
        // Compute the derivative of the guard function g along flow lines of $\dot(x)=f(x)$
        // This is given by the Lie derivative at a point x, defined as $L_{f}g(x) = (\nabla g\cdot f)(x)$
        EffectiveScalarMultivariateFunction guard_derivative = lie_derivative(guard, dynamic);
        UpperIntervalType guard_derivative_range = apply(guard_derivative, flow_bounds);
        LOGGING_PRINTLN_AT(2, "guard_derivative_range=" << guard_derivative_range);
        if (definitely(guard_derivative_range.lower_bound() > analysis_zero)) {
            // If the derivative $L_{f}g$is strictly positive over the bounding box for the flow,
            // then the guard function is strictly increasing.
            // There is at most one crossing with the guard, and the time of this
            // crossing must be the time of the event along the trajectory.
            // The crossing time $\gamma(x_0)$ given the initial state can usually be computed
            // by solving the equation $g(\phi(x_0,\gamma(x_0))) = 0$
            ValidatedScalarMultivariateFunctionPatch crossing_time;
            try {
                crossing_time = solver.implicit(compose(guard, flow), flow_spacial_domain, flow_time_domain);
                if (possibly(crossing_time.error() > 1e-8_pr)) {
                    LOGGING_PRINTLN_AT(2, event << ": crossing_time: error=" << crossing_time.error() << ", range=" << crossing_time.range());
                }
                crossings[event] = CrossingData(CrossingGeometry::INCREASING, CrossingTimeCertificate(crossing_time));
                LOGGING_PRINTLN_AT(2, "crossing_time=" << crossing_time);
            } catch (const EquationSolverException& e) {
                // If the crossing time cannot be computed, then we can still
                // use the fact that the crossing occurs as soon as $g(x(t))=0$.
                LOGGING_PRINTLN_AT(2, "Error in computing crossing time for event " << *event_iter << ": " << e.what());
                crossings[event] = CrossingData(CrossingGeometry::INCREASING);
            } catch (const std::runtime_error& e) {
                LOGGING_PRINTLN_AT(2, "Unexpected error in computing crossing time for event " << *event_iter << ": " << e.what());
                throw;
            }
        } else if (definitely(guard_derivative_range.upper_bound() < analysis_zero)) {
            // If the guard derivative is strictly negative over the bounding box for the flow,
            // then the guard function is strictly decreasing.
            // This means that the event is either initially active, or does not occur.
            // There is no need to compute a crossing time.
            crossings[event] = CrossingData(CrossingGeometry::DECREASING);
        } else {
            // If the derivative of the guard function along flow lines cannot be shown
            // to have a definite sign over the entire flow box, then try to compute
            // the sign of the second derivative $L_{f}^{2}g(x)=L_{f}L_{f}g(x)$.
            ValidatedScalarMultivariateFunction guard_second_derivative = lie_derivative(guard_derivative, dynamic);
            UpperIntervalType guard_second_derivative_bounds_range = apply(guard_second_derivative, flow_bounds);
            UpperIntervalType guard_second_derivative_flow_range = compose(guard_second_derivative, flow).range();
            UpperIntervalType guard_second_derivative_range = intersection(guard_second_derivative_bounds_range, guard_second_derivative_flow_range);
            LOGGING_PRINTLN_AT(2, "guard_second_derivative_range=" << guard_second_derivative_range);
            if (definitely(guard_second_derivative_range.lower_bound() > analysis_zero)) {
                // If the second derivative is positive, then either
                //    (i) the event is immediately active
                //   (ii) the event is never active, or
                //  (iii) the event is initially inactive, but becomes active
                //        due to a transverse crossing.
                //   (iv) the initial state is on the boundary of the guard
                //        set, possibly with the flow tangent to this set
                // It is hard compute the crossing time, even in case (iii),
                // due to the singularity due to the tangency in (iv). However,
                // we do know that in (iii), the event occurs when $t>0$ and
                // $g(\phi(x_0,t))=0$. The crossing time is not computed.
                crossings[event] = CrossingData(CrossingGeometry::CONVEX);
            } else if (definitely(guard_second_derivative_range.upper_bound() < analysis_zero)) {
                // If the second derivative is negative, then the guard
                // values $g(x(t))$ are concave along flow lines. There are
                // five main cases:
                //   (i) The event is initially active.
                //  (ii) The event is not initially active, but later becomes active.
                // (iii) The event is never active, but would become active if
                //       flowing backwards in time.
                //  (iv) The event is never active, and the maximum value of
                //       the guard along the flow lines is zero.
                // Additionally, there is the degenerate case
                //  (vi) At some point in the (forward) flow, the state touches
                //       the guard set at a point of tangency.
                // Due to the presence of the tangency, the event time is
                // not a smooth function of the initial state. Further, since
                // some trajectories cross the boundary of the guard set twice,
                // the condition $g(x(t))=0$ is not sufficient for determining
                // the jump time. A necessary and sufficient condition,
                // assuming the event is not initially active, is that
                // $g(x)=0$ and $L_{f}g(x)\geq 0$. Alternatively, we can compute
                // the <em>critical time</em> $|mu(x_0) at which the guard value
                // reaches a maximum. A necessary and sufficient condition
                // is then $g(\phi(x_0,t))=0$ and $t<= \mu(x_0)$.
                //   Note that while $g(\phi(x_0,t))=0$ and $L_{f}f(\phi(x_0,t))\geq0$ is a
                // necessary and sufficient condition for a crossing, there
                // is no necessary and sufficient condition for no crossing
                // which does not involve the critical time. A necessary and
                // sufficient condition for no crossing involving the critical
                // time is $(g(\phi(x_0,t))<=0 /\ t<=\mu(x_0)) \/ g(\phi(x_0,\mu(x_0)))<=0$
                try {
                    ValidatedScalarMultivariateFunctionPatch critical_time = solver.implicit(compose(guard_derivative, flow), flow_spacial_domain, flow_time_domain);
                    UpperIntervalType critical_time_range = critical_time.range();
                    LOGGING_PRINTLN_AT(2, "critical_time_range=" << critical_time_range);
                    if (possibly(critical_time.error() > 1e-8_pr)) {
                        LOGGING_PRINTLN_AT(2, event << ": critical_time: error=" << critical_time.error() << ", range=" << critical_time.range());
                    }

                    HybridEnclosure evolve_set_at_critical_time = initial_set;
                    evolve_set_at_critical_time.apply_space_evolve_step(flow, critical_time);
                    UpperBoxType evolve_bounds_at_critical_time = evolve_set_at_critical_time.bounding_box().euclidean_set();
                    UpperIntervalType guard_range_at_critical_time = intersection(apply(guard, evolve_bounds_at_critical_time), evolve_set_at_critical_time.range_of(guard));
                    LOGGING_PRINTLN_AT(2, "evolve_bounds_at_critical_time=" << evolve_bounds_at_critical_time);
                    LOGGING_PRINTLN_AT(2, "guard_range_at_critical_time=" << guard_range_at_critical_time);
                    if (definitely(guard_range_at_critical_time.upper_bound() < 0)) {
                        // No crossing
                        // Omit this event from the returned crossing map.
                    } else if (definitely(guard_range_at_critical_time.lower_bound() > 0)) {
                        LOGGING_PRINTLN_AT(2, "guard range eventually positive");
                        // Transverse crossing
                        // FIXME: Find a more reliable way of solving the implicit equation for the crossing time
                        //   which takes into account the fact that the derivative over the domain goes negative
                        static const Rational INTERVAL_REDUCTION_FACTOR(15, 16);
                        ValidatedScalarMultivariateFunctionPatch reduced_critical_time = INTERVAL_REDUCTION_FACTOR * critical_time;
                        HybridEnclosure evolve_set_at_reduced_critical_time = initial_set;
                        evolve_set_at_reduced_critical_time.apply_space_evolve_step(flow, reduced_critical_time);
                        HybridEnclosure evolve_set_at_upper_reduced_critical_time = initial_set;
                        evolve_set_at_upper_reduced_critical_time.apply_fixed_evolve_step(flow, static_cast<StepSizeType>(cast_exact(reduced_critical_time.range().upper_bound())));
                        LOGGING_PRINTLN_AT(2, "guard_range_at_initial_time=" << initial_set.range_of(guard));
                        LOGGING_PRINTLN_AT(2, "guard_range_at_reduced_critical_time=" << evolve_set_at_reduced_critical_time.range_of(guard));
                        LOGGING_PRINTLN_AT(2, "guard_derivative_range_at_initial_time=" << initial_set.range_of(guard_derivative));
                        LOGGING_PRINTLN_AT(2, "guard_derivative_range_at_upper_reduced_critical_time=" << evolve_set_at_upper_reduced_critical_time.range_of(guard_derivative));
                        LOGGING_PRINTLN_AT(2, "guard_derivative_range_at_reduced_critical_time=" << evolve_set_at_reduced_critical_time.range_of(guard_derivative));
                        UpperIntervalType crossing_flow_time_range(0, critical_time_range.upper_bound());
                        LOGGING_PRINTLN_AT(2, "crossing_flow_time_range=" << crossing_flow_time_range);
                        ExactIntervalType crossing_flow_time_domain = cast_exact_interval(INTERVAL_REDUCTION_FACTOR * UpperIntervalType(0, critical_time_range.upper_bound()));
                        LOGGING_PRINTLN_AT(2, "crossing_flow_time_domain=" << crossing_flow_time_domain);
                        try {
                            // Use the configured solver without changing its settings.
                            ValidatedScalarMultivariateFunctionPatch crossing_time = solver.implicit(compose(guard, flow), flow_spacial_domain, crossing_flow_time_domain);
                            UpperIntervalType crossing_time_range = crossing_time.range();
                            LOGGING_PRINTLN_AT(2, "crossing_time_range=" << crossing_time_range);
                            crossings[event] = CrossingData(CrossingGeometry::CONCAVE, CrossingTimeCertificate(crossing_time));
                        } catch (const EquationSolverException& e) {
                            LOGGING_PRINTLN_AT(2, "Error in computing crossing time over reduced time interval for event " << *event_iter << ": " << e.what());
                            // Retain the maximum certificate to block continuous evolution at critical time.
                            crossings[event] = CrossingData(GuardMaximumCertificate(
                                critical_time, guard_range_at_critical_time, evolve_bounds_at_critical_time));
                        }
                    } else {
                        crossings[event] = CrossingData(GuardMaximumCertificate(
                            critical_time, guard_range_at_critical_time, evolve_bounds_at_critical_time));
                    }
                } catch (const EquationSolverException& e) {
                    LOGGING_PRINTLN_AT(2, "Error in computing critical time for event " << *event_iter << ":\n  " << e.what());
                    crossings[event] = CrossingData(CrossingGeometry::CONCAVE);
                } catch (const std::runtime_error& e) {
                    LOGGING_PRINTLN_AT(2, "Unexpected error in computing critical time for event " << *event_iter << ":\n  " << e.what());
                    crossings[event] = CrossingData(CrossingGeometry::CONCAVE);
                }
            } else {
                // The crossing cannot be shown to be one of the kinds mentioned
                // above. A theoretically exact expression for the crossing
                // set is generally not available.
                crossings[event] = CrossingData(CrossingGeometry::UNRESOLVED);
            }
        }
    }
    if (!crossings.empty()) {
        LOGGING_PRINTLN("crossings: " << crossings);
    }
    return crossings;
}

} // namespace Ariadne::EvolutionDetail
