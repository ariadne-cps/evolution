/***************************************************************************
 *            hybrid_continuous_time_approximate_evolver.cpp
 *
 *  Copyright  2008-20  Luca Geretti, Mirko Albanese
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

/** @file src/hybrid_continuous_time_approximate_evolver.cpp
 *  @brief Implement hybrid simulator operations.
 */

#include "evolution/hybrid_continuous_time_approximate_evolver.hpp"

#include <stdexcept>

#include "algebra/algebra.hpp"
#include "evolution/hybrid_orbit.hpp"
#include "evolution/hybrid_time.hpp"
#include "function/formula.hpp"
#include "function/function.hpp"
#include "function/functional.hpp"
#include "function/taylor_model.hpp"
#include "hybrid/hybrid_set.hpp"
#include "integration/runge_kutta_integrator.hpp"
#include "logging/logging.hpp"
#include "simulation_time.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/space.hpp"
#include "symbolic/valuation.hpp"
#include "system/hybrid_automaton_interface.hpp"
#include "system/odae.hpp"
#include "system/odai.hpp"
#include "utility/array.hpp"
#include "utility/container.hpp"
#include "utility/stlio.hpp"
#include "utility/tuple.hpp"

namespace Ariadne {

template<class T>
class Orbit;

namespace {
void validate_simulation_mode(const HybridAutomatonInterface& system, DiscreteLocation location) {
    if (system.inclusion(location))
        throw std::invalid_argument("Hybrid ODAI approximate evolution requires an input sampling policy; not yet supported");
    if (const auto* equation = system.equation(location))
        if (equation->parameter_space().size() != 0u)
            throw std::invalid_argument("Hybrid approximate evolution does not yet support constant model parameters");
}
} // namespace

HybridContinuousTimeEvolver<ApproximateTag>::HybridContinuousTimeEvolver(const SystemType& system)
    : _sys_ptr(system.clone()), _configuration(new ConfigurationType()) {
}

Map<DiscreteEvent, EffectiveScalarMultivariateFunction>
    HybridContinuousTimeEvolver<ApproximateTag>::_guard_functions(const DiscreteLocation& location) const {
    Set<DiscreteEvent> events = _sys_ptr->events(location);
    Map<DiscreteEvent, EffectiveScalarMultivariateFunction> guards;
    for (Set<DiscreteEvent>::ConstIterator iter = events.begin(); iter != events.end(); ++iter) {
        EventKind kind = _sys_ptr->event_kind(location, *iter);
        if (kind != EventKind::INVARIANT && kind != EventKind::PROGRESS)
            guards.insert(*iter, _sys_ptr->guard_function(location, *iter));
    }
    return guards;
}

Bool HybridContinuousTimeEvolver<ApproximateTag>::_satisfies_invariants(const DiscreteLocation& location, const Point<FloatDPApproximation>& point, const FloatDPApproximation& time) const {
    Set<DiscreteEvent> events = _sys_ptr->events(location);
    for (const auto& event: events) {
        EventKind kind = _sys_ptr->event_kind(location, event);
        if (kind == EventKind::INVARIANT || kind == EventKind::PROGRESS) {
            // Both predicates are supplied by invariant_function with p(x)<=0.
            auto predicate = _sys_ptr->invariant_function(location, event);
            if (probably(evaluate(predicate, EvolutionDetail::simulation_arguments(Vector<FloatDPApproximation>(point), _sys_ptr->time_variable(location).has_value(), time)) > 0)) {
                return false;
            }
        }
    }
    return true;
}

template<class X>
Point<X> make_point(const HybridPoint<X>& hpt, const RealSpace& sspc) {
    if (hpt.space() == sspc) {
        return hpt.point();
    }
    const Map<RealVariable, X> values = hpt.values();
    Point<X> pt(sspc.dimension(), [&](SizeType i) { return values[sspc.variable(i)]; });
    return pt;
}

template<class X>
Void make_point_list(const HybridPoint<X>& list, Vector<Point<X>>& return_point_list, const RealSpace& sspc, SizeType index) {
    if (list.space() == sspc) {
        return_point_list[index] = list.point();
    } else {
        const Map<RealVariable, X> values = list.values();
        Point<X> pt(sspc.dimension(), [&](SizeType j) { return values[sspc.variable(j)]; });
        return_point_list[index] = pt;
    }
}

Vector<HybridPoint<FloatDPApproximation>> create_hybrid_point_list(BinaryTreePaving& grid, DiscreteLocation loc, RealSpace spc, DiscretizationHybridType _dtype, Nat fineness, const UpperBoxType& sampling_box) {
    if (_dtype == DiscretizationHybridType::HMince) {
        grid.mince(fineness);
    } else {
        grid.recombine();
    }

    BinaryTreePaving::ConstIterator iter = grid.begin();
    Vector<HybridPoint<FloatDPApproximation>> approxPointList(grid.size(), HybridPoint<FloatDPApproximation>());
    SizeType k(0);
    for (; iter != grid.end(); ++iter) {
        // Outer approximation may include cells crossing the requested boundary.
        // Sample their intersection so every initial point belongs to the box.
        auto sampled = intersection(UpperBoxType(iter->box()), sampling_box);
        auto midpoint = sampled.midpoint();
        HybridPoint<FloatDPApproximation> pt(loc, spc, midpoint);
        approxPointList.at(k) = pt;
        k++;
    }

    return approxPointList;
}

BinaryTreePaving create_hybrid_grid(UpperBoxType box) {
    Point<FloatDP> tmpPointCenter = Point(box.centre());
    Vector<FloatDP> origin(tmpPointCenter.dimension(), FloatDP(0, dp));
    Vector<FloatDP> lengths(box.widths().size(), FloatDP(0, dp));
    for (SizeType i = 0; i < tmpPointCenter.dimension(); i++) {
        origin[i] = tmpPointCenter[i];
        lengths[i] = box.widths().at(i).raw();
    }
    Grid grid(origin, lengths);
    BinaryTreePaving gridPaving(grid);

    return gridPaving;
}

template<class X>
HybridPoint<X> make_hybrid_state_auxiliary_point(const DiscreteLocation& location, const Point<X>& spt,
                                                 const RealSpace& sspc, const RealSpace& aspc, const RealSpace& saspc, const EffectiveVectorMultivariateFunction& auxiliary_function, bool timed = false, X time = X(0, dp)) {
    Point<X> sapt(saspc.dimension(), spt.zero_element());
    Point<X> apt = evaluate(auxiliary_function, EvolutionDetail::simulation_arguments(Vector<X>(spt), timed, time));
    for (SizeType i = 0; i != sapt.size(); ++i) {
        RealVariable var = saspc.variable(i);
        sapt[i] = sspc.contains(var) ? spt[sspc[var]] : apt[aspc[var]];
    }
    return HybridPoint<X>(location, saspc, sapt);
}

template<class X>
Vector<HybridPoint<X>> make_hybrid_state_auxiliary_point_list(const Vector<DiscreteLocation>& location, const Vector<Point<X>>& mpt,
                                                              const Vector<RealSpace>& sspc, const Vector<RealSpace>& aspc, const Vector<RealSpace>& saspc, const Vector<EffectiveVectorMultivariateFunction>& auxiliary_function, const HybridAutomatonInterface& system) {
    Vector<HybridPoint<X>> point_list(mpt.size(), HybridPoint<X>());
    for (SizeType index = 0; index < point_list.size(); index++) {
        Point<X> sapt(saspc[index].dimension(), mpt[index].zero_element());
        Point<X> apt = evaluate(auxiliary_function[index], EvolutionDetail::simulation_arguments(Vector<X>(mpt[index]), system.time_variable(location[index]).has_value(), X(0, dp)));
        for (SizeType i = 0; i != sapt.size(); ++i) {
            RealVariable var = saspc[index].variable(i);
            sapt[i] = sspc[index].contains(var) ? mpt.at(index)[sspc.at(index)[var]] : apt[aspc.at(index)[var]];
        }

        point_list[index] = HybridPoint<X>(location[index], saspc[index], sapt);
    }

    return point_list;
}

inline FloatDPApproximation evaluate(const EffectiveScalarMultivariateFunction& f, const Vector<FloatDPApproximation>& x) { return f(x); }
inline Vector<FloatDPApproximation> evaluate(const EffectiveVectorMultivariateFunction& f, const Vector<FloatDPApproximation>& x) { return f(x); }

/*
auto HybridContinuousTimeEvolver<ApproximateTag>::orbit(const HybridBoundedConstraintSet& initial_set, const TerminationType& termination) const
    -> Orbit<HybridApproximatePointType>
{
    DiscreteLocation first_location = *initial_set.locations().begin();
    RealSpace spc = this->_sys_ptr->continuous_state_space(first_location);
    Point<FloatDPApproximation> pt = initial_set.euclidean_set(first_location,spc).bounding_box().midpoint();
    return orbit(HybridApproximatePointType(first_location,spc,pt),termination);
}
*/

auto HybridContinuousTimeEvolver<ApproximateTag>::orbit(const HybridRealPoint& init_pt, const TerminationType& termination) const
    -> Orbit<HybridApproximatePointType> {
    return orbit(HybridApproximatePointType(init_pt, dp), termination);
}

auto HybridContinuousTimeEvolver<ApproximateTag>::orbit(const HybridApproximatePointType& init_pt,
                                                        const TerminationType& termination) const
    -> Orbit<HybridApproximatePointType> {
    LOGGING_SCOPE_CREATE;

    const HybridAutomatonInterface& system = *_sys_ptr;

    HybridTime t(0.0_x, 0);
    Dyadic h(cast_exact(configuration().step_size()));
    HybridTime tmax(termination.maximum_time(), termination.maximum_steps());

    DiscreteLocation location = init_pt.location();
    validate_simulation_mode(system, location);
    RealSpace continuous_state_space = system.continuous_state_space(location);
    RealSpace continuous_auxiliary_space = system.continuous_auxiliary_space(location);
    RealSpace continuous_state_auxiliary_space = system.state_auxiliary_space()[location];
    EffectiveVectorMultivariateFunction auxiliary_function = system.auxiliary_function(location);
    ApproximatePointType point = make_point(init_pt, continuous_state_space);
    ApproximatePointType next_point;
    List<DiscreteEvent> event_trace;

    Orbit<HybridApproximatePoint> orbit(make_hybrid_state_auxiliary_point(location, point, continuous_state_space,
                                                                          continuous_auxiliary_space, continuous_state_auxiliary_space, system.auxiliary_function(location), system.time_variable(location).has_value()));

    EffectiveVectorMultivariateFunction dynamic = system.dynamic_function(location);
    Map<DiscreteEvent, EffectiveScalarMultivariateFunction> guards = _guard_functions(location);

    RungeKutta4Integrator integrator(configuration().step_size().get_d());

    bool reached_time_horizon = false;
    while (!reached_time_horizon && possibly(t < tmax) && (event_trace.empty() || !termination.terminating_events().contains(event_trace.back()))) {
        std::streamsize old_precision = std::clog.precision();
        LOGGING_PRINTLN_AT(1,
                           "t=" << std::setw(4) << std::left << t.continuous_time().compute_get(Effort(0u), double_precision)
                                << " #e=" << std::left << t.discrete_time()
                                << " p=" << point
                                << " l=" << std::left << location
                                << " e=" << std::left << event_trace
                                << std::setprecision(old_precision));

        if (not _satisfies_invariants(location, point, FloatDPApproximation(t.continuous_time(), dp))) {
            LOGGING_PRINTLN("invariant/progress condition not satisfied, stopping evolution.");
            break;
        }

        Bool enabled = false;
        DiscreteEvent event;
        for (Map<DiscreteEvent, EffectiveScalarMultivariateFunction>::ConstIterator guard_iter = guards.begin(); guard_iter != guards.end(); ++guard_iter) {
            if (probably(evaluate(guard_iter->second, EvolutionDetail::simulation_arguments(Vector<FloatDPApproximation>(point), system.time_variable(location).has_value(), FloatDPApproximation(t.continuous_time(), dp))) > 0)) {
                enabled = true;
                event = guard_iter->first;
                break;
            }
        }

        if (enabled) {
            DiscreteLocation target = system.target(location, event);
            EffectiveVectorMultivariateFunction reset = system.reset_function(location, event);
            const auto reset_arguments = EvolutionDetail::simulation_arguments(Vector<FloatDPApproximation>(point), system.time_variable(location).has_value(), FloatDPApproximation(t.continuous_time(), dp));
            validate_simulation_mode(system, target);
            location = target;
            continuous_state_space = system.continuous_state_space(location);
            continuous_auxiliary_space = system.continuous_auxiliary_space(location);
            continuous_state_auxiliary_space = system.state_auxiliary_space()[location];
            auxiliary_function = system.auxiliary_function(location);
            next_point = reset(reset_arguments);
            event_trace.push_back(event);

            LOGGING_PRINTLN_AT(1, "event " << event << " enabled: next point " << next_point << ", on location " << target);

            dynamic = system.dynamic_function(location);
            guards = _guard_functions(location);

            t._discrete_time += 1;
        } else {
            reached_time_horizon = !probably(FloatDPApproximation(t.continuous_time(), dp) + FloatDPApproximation(h, dp) < FloatDPApproximation(tmax.continuous_time(), dp));
            const Real duration = reached_time_horizon ? tmax.continuous_time() - t.continuous_time() : Real(h);
            const FloatDPApproximation dt(duration, dp);
            next_point = system.time_variable(location)
                           ? integrator.step(dynamic, point, FloatDPApproximation(t.continuous_time(), dp), dt)
                           : integrator.step(dynamic, point, dt);
            t._continuous_time = reached_time_horizon ? tmax.continuous_time() : t._continuous_time + duration;
        }
        point = next_point;
        orbit.insert(t, make_hybrid_state_auxiliary_point(location, point,
                                                          continuous_state_space, continuous_auxiliary_space, continuous_state_auxiliary_space, auxiliary_function, system.time_variable(location).has_value(), FloatDPApproximation(t.continuous_time(), dp)));
    }

    return orbit;
}

auto HybridContinuousTimeEvolver<ApproximateTag>::orbit(const HybridBoundedConstraintSet& initial_set, const TerminationType& termination) const
    -> OrbitListType {
    DiscreteLocation first_location = *initial_set.locations().begin();
    RealSpace spc = this->_sys_ptr->continuous_state_space(first_location);
    Box<UpperIntervalType> box = initial_set.bounding_box().euclidean_set(first_location);
    HybridUpperBox Hbox = HybridUpperBox(first_location, spc, box);
    return orbit(Hbox, first_location, this->_sys_ptr->state_space(), termination);
}
/*
auto HybridContinuousTimeEvolver<ApproximateTag>::orbit(const Vector<HybridRealPoint>& initial_set, const TerminationType& termination) const
    -> Orbit<HybridApproximateListPointType>
{
    HybridApproximateListPointType point_list = HybridApproximateListPointType(initial_set.size(),HybridApproximatePointType(initial_set.at(0),dp));
    for(SizeType i=0; i<initial_set.size(); i++)
    {
        point_list.at(i) = HybridApproximatePointType(initial_set.at(i),dp);
    }
    return orbit(point_list,termination);
}
*/

auto HybridContinuousTimeEvolver<ApproximateTag>::orbit(const HybridUpperBox& initial_box, DiscreteLocation loc, HybridSpace space, const TerminationType& termination) const
    -> Orbit<HybridApproximateListPointType> {
    auto box = initial_box.euclidean_set(space[loc]);
    auto lengths = box.widths();
    Nat box_width_null = 0;
    for (SizeType i = 0; i < lengths.size(); ++i) {
        if (lengths[i].get_d() > 0) {
            continue;
        }
        ++box_width_null;
    }
    if (box_width_null == lengths.size()) {
        auto midpoint = box.midpoint();
        HybridApproximateListPointType pointList(1, HybridApproximatePointType(loc, space[loc], midpoint));
        return orbit(pointList, termination);
    }
    if (box_width_null > 0) {
        // A partially flat box needs positive grid lengths in every coordinate.
        box = widen(box, FloatDPUpperBound(0.0001_q, dp));
    }
    BinaryTreePaving gridPaving = create_hybrid_grid(box);
    gridPaving.adjoin_outer_approximation(box, _configuration->fineness());
    HybridApproximateListPointType pointList = create_hybrid_point_list(gridPaving, loc, this->_sys_ptr->continuous_state_space(loc), _configuration->d_type(), _configuration->fineness(), initial_box.euclidean_set(space[loc]));
    return orbit(pointList, termination);
}

auto HybridContinuousTimeEvolver<ApproximateTag>::orbit(const HybridApproximateListPointType& init_list,
                                                        const TerminationType& termination) const
    -> Orbit<HybridApproximateListPointType> {
    LOGGING_SCOPE_CREATE;

    const HybridAutomatonInterface& system = *_sys_ptr;

    Vector<HybridTime> t(init_list.size(), HybridTime(0.0_x, 0));
    Dyadic h(cast_exact(configuration().step_size()));
    HybridTime tmax(termination.maximum_time(), termination.maximum_steps());

    Vector<ApproximatePointType> point_list(init_list.size(), ApproximatePointType());

    Vector<ApproximatePointType> next_point(init_list.size(), ApproximatePointType());
    Vector<List<DiscreteEvent>> event_trace(init_list.size(), List<DiscreteEvent>());

    Vector<RealSpace> continuous_state_space(init_list.size(), RealSpace());
    Vector<RealSpace> continuous_auxiliary_space(init_list.size(), RealSpace());
    Vector<RealSpace> continuous_state_auxiliary_space(init_list.size(), RealSpace());
    Vector<EffectiveVectorMultivariateFunction> auxiliary_function(init_list.size(), EffectiveVectorMultivariateFunction());
    Vector<DiscreteLocation> location(init_list.size(), DiscreteLocation());

    for (SizeType i = 0; i < init_list.size(); i++) {
        location[i] = init_list[i].location();
        validate_simulation_mode(system, location[i]);
        continuous_state_space[i] = system.continuous_state_space(location[i]);
        continuous_auxiliary_space[i] = system.continuous_auxiliary_space(location[i]);
        continuous_state_auxiliary_space[i] = system.state_auxiliary_space()[location[i]];
        auxiliary_function[i] = system.auxiliary_function(location[i]);
        make_point_list(init_list[i], point_list, continuous_state_space[i], i);
    }

    OrbitListType orbitList(make_hybrid_state_auxiliary_point_list(location, point_list, continuous_state_space,
                                                                   continuous_auxiliary_space, continuous_state_auxiliary_space, auxiliary_function, system));

    for (SizeType i = 0; i < init_list.size(); i++) {
        bool reached_time_horizon = false;
        while (!reached_time_horizon && possibly(t[i] < tmax) && (event_trace[i].empty() || !termination.terminating_events().contains(event_trace[i].back()))) {
            EffectiveVectorMultivariateFunction dynamic = system.dynamic_function(location[i]);
            Map<DiscreteEvent, EffectiveScalarMultivariateFunction> guards = _guard_functions(location[i]);

            RungeKutta4Integrator integrator(configuration().step_size().get_d());

            std::streamsize old_precision = std::clog.precision();
            LOGGING_PRINTLN_AT(1,
                               "t=" << std::setw(4) << std::left << t[i].continuous_time().compute_get(Effort(0u), double_precision)
                                    << " #e=" << std::left << t[i].discrete_time()
                                    << " p=" << point_list[i]
                                    << " l=" << std::left << location[i]
                                    << " e=" << std::left << event_trace[i]
                                    << std::setprecision(old_precision));

            if (not _satisfies_invariants(location[i], point_list[i], FloatDPApproximation(t[i].continuous_time(), dp))) {
                LOGGING_PRINTLN("invariant/progress condition not satisfied, stopping evolution.");
                break;
            }

            Bool enabled = false;
            DiscreteEvent event;

            for (Map<DiscreteEvent, EffectiveScalarMultivariateFunction>::ConstIterator guard_iter = guards.begin(); guard_iter != guards.end(); ++guard_iter) {
                if (probably(evaluate(guard_iter->second, EvolutionDetail::simulation_arguments(Vector<FloatDPApproximation>(point_list[i]), system.time_variable(location[i]).has_value(), FloatDPApproximation(t[i].continuous_time(), dp))) > 0)) {
                    enabled = true;
                    event = guard_iter->first;
                    break;
                }
            }

            if (enabled) {
                DiscreteLocation target = system.target(location[i], event);
                EffectiveVectorMultivariateFunction reset = system.reset_function(location[i], event);
                const auto reset_arguments = EvolutionDetail::simulation_arguments(Vector<FloatDPApproximation>(point_list[i]), system.time_variable(location[i]).has_value(), FloatDPApproximation(t[i].continuous_time(), dp));
                validate_simulation_mode(system, target);
                location[i] = target;
                continuous_state_space[i] = system.continuous_state_space(location[i]);
                continuous_auxiliary_space[i] = system.continuous_auxiliary_space(location[i]);
                continuous_state_auxiliary_space[i] = system.state_auxiliary_space()[location[i]];
                auxiliary_function[i] = system.auxiliary_function(location[i]);
                next_point[i] = reset(reset_arguments);
                event_trace[i].push_back(event);

                LOGGING_PRINTLN_AT(1, "event " << event << " enabled: next point " << next_point[i] << ", on location " << target);

                dynamic = system.dynamic_function(location[i]);
                guards = _guard_functions(location[i]);

                t[i]._discrete_time += 1;
            } else {
                reached_time_horizon = !probably(FloatDPApproximation(t[i].continuous_time(), dp) + FloatDPApproximation(h, dp) < FloatDPApproximation(tmax.continuous_time(), dp));
                const Real duration = reached_time_horizon ? tmax.continuous_time() - t[i].continuous_time() : Real(h);
                const FloatDPApproximation dt(duration, dp);
                next_point[i] = system.time_variable(location[i])
                                  ? integrator.step(dynamic, point_list[i], FloatDPApproximation(t[i].continuous_time(), dp), dt)
                                  : integrator.step(dynamic, point_list[i], dt);
                t[i]._continuous_time = reached_time_horizon ? tmax.continuous_time() : t[i]._continuous_time + duration;
            }

            point_list[i] = next_point[i];

            orbitList.insert(t[i], make_hybrid_state_auxiliary_point(location[i], point_list[i], continuous_state_space[i], continuous_auxiliary_space[i], continuous_state_auxiliary_space[i], auxiliary_function[i], system.time_variable(location[i]).has_value(), FloatDPApproximation(t[i].continuous_time(), dp)), i);
        }
    }

    return orbitList;
}

HybridContinuousTimeEvolverConfiguration<ApproximateTag>::HybridContinuousTimeEvolverConfiguration() : _step_size(0.125, double_precision), _fineness(0), _discretization_type(DiscretizationHybridType::HMince) {}

OutputStream& HybridContinuousTimeEvolverConfiguration<ApproximateTag>::_write(OutputStream& os) const {
    os << "HybridContinuousTimeEvolverConfiguration<ApproximateTag>("
       << "\n  step_size=" << step_size()
       << "\n  fineness=" << fineness()
       << "\n  discretization type=" << d_type()
       << "\n)";
    return os;
}

} // namespace Ariadne

