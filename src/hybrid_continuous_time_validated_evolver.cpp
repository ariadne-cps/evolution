/***************************************************************************
 *            hybrid_continuous_time_validated_evolver.cpp
 *
 *  Copyright  2009-20  Pieter Collins
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

/** @file src/hybrid_continuous_time_validated_evolver.cpp
 *  @brief Implement hybrid evolver evolution operations.
 */

#include "evolution/hybrid_continuous_time_validated_evolver.hpp"

#include <algorithm>
#include <typeinfo>

#include "algebra/polynomial.hpp"
#include "algebra/vector.hpp"
#include "continuous_integration.hpp"
#include "evolution/hybrid_continuous_time_validated_evolver.hpp"
#include "evolution/hybrid_time.hpp"
#include "evolution/orbit.hpp"
#include "evolution_time.hpp"
#include "function/function.hpp"
#include "function/function_model.hpp"
#include "function/function_patch.hpp"
#include "hybrid_crossing_analysis.hpp"
#include "hybrid_guard_predicate.hpp"
#include "hybrid_outer_flow.hpp"
#include "hybrid_parameters.hpp"
#include "hybrid_step_results.hpp"
#include "hybrid_timing.hpp"
#include "hybrid_transition_handlers.hpp"
#include "integration/integrator.hpp"
#include "integration/selected_inclusion_integrator.hpp"
#include "integration/validated_step_integrator.hpp"
#include "numeric/numeric.hpp"
#include "paving/binary_tree_paving.hpp"
#include "system/continuous_time_hybrid_automaton_interface.hpp"
#include "system/odae.hpp"
#include "system/odai.hpp"

#include "equation-solving/equation_solver.hpp"
#include "equation-solving/interval_newton_method.hpp"

namespace {

} // namespace

namespace Ariadne {

// TODO: Move this functionality to Numeric module
inline PositiveValidatedUpperNumber abs(PositiveValidatedUpperNumber y) { return y; }
inline PositiveValidatedUpperNumber mag(PositiveValidatedUpperNumber y) { return y; }

static const FloatDP hybrid_zero = {0, dp};

inline auto operator+(Int n, FloatDPBounds x) -> decltype(FloatDP(n, dp) + x) { return FloatDP(n, dp) + x; }

template<class PR>
inline auto operator*(FloatUpperBound<PR> x1, PositiveFloat<PR> x2) -> FloatUpperBound<PR> {
    return FloatUpperBound<PR>(mul(up, x1._u, x2._v));
}
template<class PR>
inline auto operator*(PositiveFloatUpperBound<PR> x1, PositiveFloat<PR> x2) -> PositiveFloatUpperBound<PR> {
    return PositiveFloatUpperBound<PR>(mul(up, x1._u, x2._v));
}

inline auto operator>(FloatDPLowerBound x, Real r) -> decltype(x > declval<FloatDPBounds>()) { return x > FloatDPBounds(r, dp); }
inline auto operator>=(FloatDPLowerBound x, Real r) -> decltype(x >= declval<FloatDPBounds>()) { return x >= FloatDPBounds(r, dp); }
inline auto operator>(FloatDPUpperBound x, Real r) -> decltype(x > declval<FloatDPBounds>()) { return x > FloatDPBounds(r, dp); }
inline auto operator<=(FloatDPUpperBound x, Real r) -> decltype(x <= declval<FloatDPBounds>()) { return x <= FloatDPBounds(r, dp); }

static const DiscreteEvent final_event("_tmax_");

typedef Vector<FloatDP> ExactFloatVector;
typedef Vector<ExactIntervalType> ExactIntervalVectorType;

Set<DiscreteEvent> blocking_events(const Map<DiscreteEvent, TransitionData>& transitions);
Set<DiscreteEvent> activating_events(const Map<DiscreteEvent, TransitionData>& transitions);

OutputStream& operator<<(OutputStream& os, const StepKind& stpk) {
    switch (stpk) {
        case StepKind::CONSTANT_EVOLUTION_TIME:
            os << "CONSTANT_EVOLUTION_TIME";
            break;
        case StepKind::SPACE_DEPENDENT_EVOLUTION_TIME:
            os << "SPACE_DEPENDENT_EVOLUTION_TIME";
            break;
        case StepKind::TIME_DEPENDENT_EVOLUTION_TIME:
            os << "TIME_DEPENDENT_EVOLUTION_TIME";
            break;
        case StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME:
            os << "SPACETIME_DEPENDENT_EVOLUTION_TIME";
            break;
        case StepKind::PARAMETER_DEPENDENT_EVOLUTION_TIME:
            os << "PARAMETER_DEPENDENT_EVOLUTION_TIME";
            break;
        case StepKind::PARAMETER_DEPENDENT_FINISHING_TIME:
            os << "PARAMETER_DEPENDENT_FINISHING_TIME";
            break;
        case StepKind::SPACETIME_DEPENDENT_FINISHING_TIME:
            os << "SPACETIME_DEPENDENT_FINISHING_TIME";
            break;
        case StepKind::CONSTANT_FINISHING_TIME:
            os << "CONSTANT_FINISHING_TIME";
            break;
        default:
            os << "UNKNOWN";
            break;
    }
    return os;
}

OutputStream& operator<<(OutputStream& os, const FinishingKind& fnshk) {
    switch (fnshk) {
        case FinishingKind::BEFORE_FINAL_TIME:
            os << "BEFORE_FINAL_TIME";
            break;
        case FinishingKind::AT_FINAL_TIME:
            os << "AT_FINAL_TIME";
            break;
        case FinishingKind::AFTER_FINAL_TIME:
            os << "AFTER_FINAL_TIME";
            break;
        case FinishingKind::STRADDLE_FINAL_TIME:
            os << "STRADDLE_FINAL_TIME";
            break;
        default:
            os << "UNKNOWN";
            break;
    }
    return os;
}

OutputStream& operator<<(OutputStream& os, const TransitionData& transition) {
    os << "kind=" << transition.event_kind << ", guard=" << transition.guard_function << ", "
                                                                                         "target="
       << transition.target;
    if (is_activating(transition.event_kind)) {
        os << ", reset=" << transition.reset_function;
    }
    return os;
}

OutputStream& operator<<(OutputStream& os, const TimingData& timing) {
    os << "step_kind=" << timing.step_kind << ", finishing_kind=" << timing.finishing_kind << ", step_size=" << timing.step_size << ", "
       << "final_time=" << timing.final_time;
    if (timing.step_kind == StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME) {
        os << ", spacial_evolution_time=" << timing.spacetime_dependent_evolution_time;
    } else if (timing.step_kind == StepKind::PARAMETER_DEPENDENT_FINISHING_TIME) {
        os << ", parameter_dependent_finishing_time=" << timing.parameter_dependent_finishing_time;
    }
    os << ", parameter_dependent_evolution_time=" << timing.parameter_dependent_evolution_time;
    return os;
}

Bool is_blocking(EventKind evk) { return EvolutionDetail::transition_policy(evk).blocks_flow; }

Bool is_activating(EventKind evk) { return EvolutionDetail::transition_policy(evk).activates_reset; }

// Extract the blocking events.
Set<DiscreteEvent> blocking_events(const Map<DiscreteEvent, TransitionData>& transitions) {
    Set<DiscreteEvent> events;
    for (Map<DiscreteEvent, TransitionData>::ConstIterator transition_iter = transitions.begin();
         transition_iter != transitions.end(); ++transition_iter) {
        if (is_blocking(transition_iter->second.event_kind)) {
            events.insert(transition_iter->first);
        }
    }
    return events;
}

// Extract the activating events.
Set<DiscreteEvent> activating_events(const Map<DiscreteEvent, TransitionData>& transitions) {
    Set<DiscreteEvent> events;
    for (Map<DiscreteEvent, TransitionData>::ConstIterator transition_iter = transitions.begin();
         transition_iter != transitions.end(); ++transition_iter) {
        if (is_activating(transition_iter->second.event_kind)) {
            events.insert(transition_iter->first);
        }
    }
    return events;
}

Orbit<HybridEnclosure>
    HybridContinuousTimeEvolverBase<ValidatedTag>::
        orbit(const HybridExactBoxType& initial_box,
              const HybridTerminationCriterion& termination,
              Semantics semantics) const {
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN_AT(1, "initial_box=" << initial_box);
    HybridEnclosure initial_enclosure(initial_box, EnclosureConfiguration(this->function_factory()));
    EvolutionDetail::prepare_hybrid_mode(this->system(), initial_enclosure);
    LOGGING_PRINTLN_AT(1, "initial_enclosure=" << initial_enclosure);
    return this->orbit(initial_enclosure, termination, semantics);
}

Orbit<HybridEnclosure>
    HybridContinuousTimeEvolverBase<ValidatedTag>::
        orbit(const HybridBoxSet& initial_box,
              const HybridTerminationCriterion& termination,
              Semantics semantics) const {
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN_AT(1, "initial_box=" << initial_box);
    HybridEnclosure initial_enclosure(initial_box, this->system().continuous_state_space(initial_box.location()), EnclosureConfiguration(this->function_factory()));
    EvolutionDetail::prepare_hybrid_mode(this->system(), initial_enclosure);
    LOGGING_PRINTLN_AT(1, "initial_enclosure=" << initial_enclosure);
    return this->orbit(initial_enclosure, termination, semantics);
}

Orbit<HybridEnclosure>
    HybridContinuousTimeEvolverBase<ValidatedTag>::
        orbit(const HybridBoundedConstraintSet& initial_set,
              const HybridTerminationCriterion& termination,
              Semantics semantics) const {
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN_AT(1, "initial_set=" << initial_set);
    HybridEnclosure initial_enclosure(initial_set, this->system().continuous_state_space(initial_set.location()), EnclosureConfiguration(this->function_factory()));
    EvolutionDetail::prepare_hybrid_mode(this->system(), initial_enclosure);
    LOGGING_PRINTLN_AT(1, "initial_enclosure=" << initial_enclosure);
    return this->orbit(initial_enclosure, termination, semantics);
}

Orbit<HybridEnclosure>
    HybridContinuousTimeEvolverBase<ValidatedTag>::
        orbit(const HybridEnclosure& initial,
              const HybridTerminationCriterion& termination,
              Semantics semantics) const {
    this->_validate_mode(initial.location(), semantics);
    ARIADNE_PRECONDITION(this->system().state_auxiliary_space().has_location(initial.location()))
    ARIADNE_PRECONDITION(this->system().continuous_state_space(initial.location()) == initial.state_space())
    ARIADNE_PRECONDITION(initial.auxiliary_space().size() == 0u || this->system().continuous_auxiliary_space(initial.location()) == initial.auxiliary_space())
    HybridEnclosure prepared = initial;
    EvolutionDetail::prepare_hybrid_mode(this->system(), prepared);
    auto backend = _backend_factory_ptr->create(prepared, termination.maximum_time());
    backend->set_semantics(semantics);
    backend->run([this, &termination, semantics](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics& statistics) {
        return this->_process_working_set(item, termination, semantics, statistics);
    });
    return backend->release_orbit();
}

HybridContinuousTimeEvolverBase<ValidatedTag>::HybridContinuousTimeEvolverBase(const SystemType& system) {
    this->_create(system, FunctionFactoryType(make_taylor_function_patch_factory()));
}

HybridContinuousTimeEvolverBase<ValidatedTag>::HybridContinuousTimeEvolverBase(const SystemType& system,
                                                                               const FunctionFactoryType& factory) {
    this->_create(system, factory);
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::_create(
    const SystemType& system,
    const FunctionFactoryType& factory) {
    this->_sys_ptr = std::shared_ptr<SystemType>(system.clone());
    this->_function_factory_ptr = factory.managed_pointer();
    this->_solver_ptr = std::shared_ptr<EquationSolverInterface>(new IntervalNewtonSolver(1e-8, 12));
    this->_backend_factory_ptr = std::make_shared<FlatHybridEvolutionBackendFactory>();
    this->ALLOW_CREEP = true;
    this->ALLOW_UNWIND = false;
}

const HybridContinuousTimeEvolverBase<ValidatedTag>::SystemType&
    HybridContinuousTimeEvolverBase<ValidatedTag>::system() const {
    return *this->_sys_ptr;
}

HybridContinuousTimeEvolverBase<ValidatedTag>::ConfigurationType&
    HybridContinuousTimeEvolverBase<ValidatedTag>::configuration() {
    return *this->_configuration_ptr;
}

const HybridContinuousTimeEvolverBase<ValidatedTag>::ConfigurationType&
    HybridContinuousTimeEvolverBase<ValidatedTag>::configuration() const {
    return *this->_configuration_ptr;
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::set_function_factory(const FunctionFactoryType& factory) {
    this->_function_factory_ptr = factory.managed_pointer();
}

const HybridContinuousTimeEvolverBase<ValidatedTag>::FunctionFactoryType
    HybridContinuousTimeEvolverBase<ValidatedTag>::function_factory() const {
    return FunctionFactoryType(this->_function_factory_ptr);
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::set_integrator(const IntegratorInterface& integrator) {
    this->_integrator_ptr = std::shared_ptr<IntegratorInterface>(integrator.clone());
    this->_inclusion_integrator_ptr.reset();
    this->_box_integrator_ptr.reset();
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::set_integrator(const IntegratorBase& base, const InputApproximation& reach,
                                                                   const InputApproximation& evolve) {
    if (dynamic_cast<const PreconditionedGradedTaylorSeriesIntegrator*>(&base) != nullptr)
        throw std::invalid_argument("Hybrid inclusion input approximations do not support preconditioning");
    auto selected = std::make_shared<SelectedParametricMultifunctionIntegrator>(base, reach, evolve);
    this->set_integrator(static_cast<const IntegratorInterface&>(base));
    this->_inclusion_integrator_ptr = std::move(selected);
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::set_integrator(const FirstOrderMultifunctionIntegrator& boxes) {
    this->_box_integrator_ptr = std::make_shared<FirstOrderMultifunctionIntegrator>(boxes);
    this->_inclusion_integrator_ptr.reset();
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::_validate_mode(DiscreteLocation location, Semantics semantics) const {
    if (this->system().time_variable(location)) {
        if (semantics != Semantics::UPPER)
            throw std::invalid_argument("Timed hybrid dynamics/events currently require UPPER semantics");
        if (dynamic_cast<const PreconditionedGradedTaylorSeriesIntegrator*>(this->_integrator_ptr.get()) != nullptr)
            throw std::invalid_argument("Timed hybrid evolution does not support preconditioning");
    }
    if (EvolutionDetail::hybrid_parameter_space(this->system(), location).dimension() != 0u) {
        if (semantics != Semantics::UPPER)
            throw std::invalid_argument("Parameterized hybrid modes require UPPER semantics");
        if (dynamic_cast<const PreconditionedGradedTaylorSeriesIntegrator*>(this->_integrator_ptr.get()) != nullptr)
            throw std::invalid_argument("Parameterized hybrid evolution does not support preconditioning");
        EvolutionDetail::validate_outer_event_kinds(this->system(), location);
    }
    if (const auto* inclusion = this->system().inclusion(location)) {
        if (semantics != Semantics::UPPER)
            throw std::invalid_argument("Hybrid inclusion modes require UPPER semantics");
        // UPPER may restrict reach by PROGRESS boundaries and construct
        // conservative URGENT/PERMISSIVE/IMPACT candidates. Boundary-only
        // impact filtering retains tangencies; LOWER remains unsupported.
        EvolutionDetail::validate_outer_event_kinds(this->system(), location);
        if (!this->_inclusion_integrator_ptr && !this->_box_integrator_ptr)
            throw std::invalid_argument("Select an explicit integration method for hybrid inclusion modes");
        if (this->_inclusion_integrator_ptr && inclusion->has_input_function() && inclusion->number_of_inputs() == 0u)
            throw std::invalid_argument("The selected hybrid inclusion method requires free numerical inputs");
    }
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::set_solver(const EquationSolverInterface& solver) {
    this->_solver_ptr = std::shared_ptr<EquationSolverInterface>(solver.clone());
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::set_evolution_backend_factory(const HybridEvolutionBackendFactoryInterface& factory) {
    _backend_factory_ptr = std::shared_ptr<const HybridEvolutionBackendFactoryInterface>(factory.clone());
}

HybridContinuousTimeEvolverBase<ValidatedTag>::EnclosureType
    HybridContinuousTimeEvolverBase<ValidatedTag>::enclosure(const HybridExactBox& initial_box) const {
    HybridEnclosure result(initial_box, EnclosureConfiguration(this->function_factory()));
    EvolutionDetail::prepare_hybrid_mode(this->system(), result);
    return result;
}

HybridContinuousTimeEvolverBase<ValidatedTag>::EnclosureType
    HybridContinuousTimeEvolverBase<ValidatedTag>::enclosure(const HybridBoundedConstraintSet& initial_set) const {
    HybridEnclosure result(initial_set, this->system().continuous_state_space(initial_set.location()), EnclosureConfiguration(this->function_factory()));
    EvolutionDetail::prepare_hybrid_mode(this->system(), result);
    return result;
}

struct EvolutionStepData {
    EvolutionStepData() : progress(false), finishing(false), events() {}
    Bool progress;
    Bool finishing;
    Set<DiscreteEvent> events;
};

Void HybridContinuousTimeEvolverBase<ValidatedTag>::
    _log_summary(const HybridEnclosure& starting_set, const HybridEvolutionStatistics& statistics) const {
    UpperBoxType starting_bounding_box = starting_set.state_bounding_box();
    UpperIntervalType starting_time_range = starting_set.time_range();
    UpperIntervalType starting_dwell_time_range = starting_set.dwell_time_range();
    std::streamsize old_precision = std::clog.precision();
    LOGGING_PRINTLN(
        "#r=" << std::setw(5) << std::left << statistics.reach
              << " #i=" << std::setw(5) << std::left << statistics.intermediate
              << " #f=" << std::setw(4) << std::left << statistics.final
              << " #e=" << std::setw(3) << std::left << starting_set.previous_events().size()
              << " #p=" << std::setw(2) << std::left << starting_set.number_of_parameters()
              << " #c=" << std::setw(1) << std::left << starting_set.number_of_constraints()
              << " t=[" << std::setw(6) << std::setprecision(3) << std::left << std::fixed << starting_time_range.lower_bound()
              << "," << std::setw(6) << std::left << std::fixed << starting_time_range.upper_bound() << "]" << std::flush
              << " dwt=[" << std::setw(6) << std::setprecision(3) << std::left << std::fixed << starting_dwell_time_range.lower_bound()
              << "," << std::setw(6) << std::left << std::fixed << starting_dwell_time_range.upper_bound() << "]" << std::flush
              << " c=" << starting_bounding_box.centre()
              << " r=" << std::setw(4) << std::fixed << std::setprecision(3) << starting_bounding_box.radius()
              << " te=" << std::setw(7) << std::scientific << std::setprecision(1) << starting_set.time_function().error() << std::flush
              << " se=" << std::setw(7) << std::scientific << std::setprecision(1) << sup_norm(starting_set.state_function().errors()) << std::fixed << std::flush
              << " l=" << std::left << starting_set.location()
              << " e=" << starting_set.previous_events() << std::setprecision(old_precision));
}

Map<DiscreteEvent, TransitionData>
    HybridContinuousTimeEvolverBase<ValidatedTag>::
        _extract_transitions(const DiscreteLocation& location) const {
    LOGGING_SCOPE_CREATE;
    Map<DiscreteEvent, TransitionData> transitions;
    Set<DiscreteEvent> events = this->system().events(location);
    if (events.empty())
        return transitions;
    // An inclusion has no single-valued field for a Lie derivative.
    // Extract supported guards and non-reset blockers without one.
    if (this->system().inclusion(location) != nullptr || this->system().time_variable(location) || EvolutionDetail::hybrid_parameter_space(this->system(), location).dimension() != 0u)
        return EvolutionDetail::outer_flow_transitions(this->system(), location);

    const EffectiveVectorMultivariateFunction& dynamic = this->system().dynamic_function(location);
    for (Set<DiscreteEvent>::ConstIterator event_iter = events.begin();
         event_iter != events.end(); ++event_iter) {
        DiscreteEvent event = *event_iter;
        EventKind event_kind = this->system().event_kind(location, event);
        LOGGING_PRINTLN_AT(1, "event=" << event << ", kind=" << event_kind);
        EffectiveScalarMultivariateFunction constraint_function;
        if (is_activating(event_kind)) {
            constraint_function = this->system().guard_function(location, event);
        } else {
            constraint_function = this->system().invariant_function(location, event);
        }
        LOGGING_PRINTLN_AT(1, "constraint_function=" << constraint_function);
        EffectiveScalarMultivariateFunction constraint_flow_derivative_function = lie_derivative(constraint_function, dynamic);
        EffectiveVectorMultivariateFunction reset_function;
        DiscreteLocation target;
        RealSpace target_space;
        EffectiveVectorMultivariateFunction target_auxiliary_function;
        RealSpace target_auxiliary_space;
        if (is_activating(event_kind)) {
            reset_function = this->system().reset_function(location, event);
            target = this->system().target(location, event);
            target_space = this->system().continuous_state_space(target);
            target_auxiliary_function = this->system().auxiliary_function(target);
            target_auxiliary_space = this->system().continuous_auxiliary_space(target);
        }
        TransitionData transition_data = {event, event_kind, constraint_function, constraint_flow_derivative_function, target, reset_function, target_space, target_auxiliary_function, target_auxiliary_space};
        if (is_activating(event_kind))
            EvolutionDetail::configure_transition_target(this->system(), transition_data);
        transitions.insert(event, transition_data);
    }
    LOGGING_PRINTLN("transitions=" << transitions);
    return transitions;
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::
    _apply_invariants(HybridEnclosure& initial_set,
                      const Map<DiscreteEvent, TransitionData>& transitions) const {
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN_AT(1, "initial_set=" << initial_set);
    LOGGING_PRINTLN("transitions=" << transitions);
    EvolutionDetail::apply_state_invariants(initial_set, transitions);
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::
    _process_starting_events(const HybridEnclosure& initial_set,
                             const Map<DiscreteEvent, TransitionData>& transitions,
                             Semantics semantics,
                             HybridEvolutionStepResult& result) const {
    LOGGING_SCOPE_CREATE;
    HybridEnclosure invariant_set = initial_set;

    // Apply restrictions due to invariants
    this->_apply_invariants(invariant_set, transitions);

    // Set the flowable set, storing the invariant set as a base for jumps
    HybridEnclosure flowable_set = invariant_set;

    // Compute possibly initially active events
    Set<DiscreteEvent> events = transitions.keys();
    for (Map<DiscreteEvent, TransitionData>::ConstIterator transition_iter = transitions.begin();
         transition_iter != transitions.end(); ++transition_iter) {
        DiscreteEvent event = transition_iter->first;
        const TransitionData& transition = transition_iter->second;
        // Assume that an impact cannot occur immediately after any other event.
        // This is true after the impact, since $L_{f}g$ changes sign.
        // After a different event, this is not so clear, though it should be
        // a modelling error for the set to intersect the interior of the guard
        // FIXME: Need to consider impact which really can occur immediately due to mapping to boundary.
        // TODO: The condition for an impact to occur is $L_{f}g>0$.
        //       Since this is now available, we should implement this!
        const auto policy = EvolutionDetail::transition_policy(transition.event_kind);
        if (policy.inspect_initial_activation) {
            const auto activation = EvolutionDetail::classify_initial_guard(invariant_set, transition);
            if (activation != EvolutionDetail::HybridInitialGuardActivation::EXCLUDED) {
                if (policy.resets_initially && EvolutionDetail::permits_initial_reset(activation, semantics)) {
                    HybridEnclosure immediate_jump_set = invariant_set;
                    EvolutionDetail::restrict_source_guard(
                        immediate_jump_set, event, transition.guard_function, transition.source_arguments);
                    if (!definitely(immediate_jump_set.is_empty())) {
                        // UPPER retains possible (including spurious) jumps.
                        // LOWER never selects only the guard-active part of a
                        // source whose whole activation has not been certified.
                        LOGGING_PRINTLN_AT(1, event << ": " << transition.event_kind << ", immediate");
                        result.record_reach(immediate_jump_set);
                        EvolutionDetail::apply_transition_reset(immediate_jump_set, transition);
                        LOGGING_PRINTLN_AT(1, "immediate_jump_set=" << immediate_jump_set);
                        result.record_intermediate(immediate_jump_set);
                        const auto resolution = policy.blocks_flow &&
                                                        activation == EvolutionDetail::HybridInitialGuardActivation::FULL_SOURCE
                                                  ? HybridEvolutionTransitionResolution::FULL_SOURCE
                                                  : HybridEvolutionTransitionResolution::UNRESOLVED;
                        result.transition_to(immediate_jump_set, event, transition.event_kind, std::nullopt, resolution);
                    }
                }
                // Preserve the existing blocking restriction; refinement of
                // a partition is not a substitute for universal activation.
                if (policy.blocks_flow) {
                    flowable_set.new_invariant(event, transition.guard_function, transition.source_arguments);
                }
            }
        }
    }

    // Put the flowable set in the starting sets for ordinary evolution
    if (!definitely(flowable_set.is_empty())) {
        result.continue_with(flowable_set);
    }
}

ValidatedIntegrationStep
    HybridContinuousTimeEvolverBase<ValidatedTag>::
        _compute_integration_step(EffectiveVectorMultivariateFunction dynamic,
                                  const ExactBoxType& initial_box,
                                  const StepSizeType& maximum_step_size) const {
    LOGGING_SCOPE_CREATE;

    IntegratorInterface& integrator = *this->_integrator_ptr;

    // The shared adapter preserves adaptive error control and caps the
    // certified interval at the requested maximum, even for overshooting solvers.
    const FlowStepIntegratorAdapter step_integrator(integrator);
    return step_integrator.step(ValidatedIntegrationRequest(
        dynamic, initial_box, suggest(maximum_step_size)));
}

Set<DiscreteEvent>
    HybridContinuousTimeEvolverBase<ValidatedTag>::
        _compute_active_events(const EffectiveVectorMultivariateFunction& dynamic,
                               const Map<DiscreteEvent, EffectiveScalarMultivariateFunction>& guards,
                               const ValidatedVectorMultivariateFunctionPatch& flow,
                               const HybridEnclosure& starting_set) const {
    LOGGING_SCOPE_CREATE;
    Set<DiscreteEvent> events = guards.keys();
    Set<DiscreteEvent> active_events;
    HybridEnclosure reach_set = starting_set;
    UpperBoxType flow_bounds = flow.range();
    LOGGING_PRINTLN("flow_bounds=" << flow_bounds);
    reach_set.apply_full_reach_step(flow);
    LOGGING_PRINTLN_AT(1, "reach_set=" << reach_set);
    for (Set<DiscreteEvent>::Iterator event_iter = events.begin(); event_iter != events.end(); ++event_iter) {
        const DiscreteEvent event = *event_iter;
        const EffectiveScalarMultivariateFunction& guard_function = guards[event];
        EffectiveScalarMultivariateFunction flow_derivative = lie_derivative(guard_function, dynamic);
        LOGGING_PRINTLN_AT(1, "event=" << event << ", guard=" << guard_function << ", flow_derivative=" << flow_derivative);
        // First try a simple test based on the bounding box
        UpperIntervalType guard_range = apply(guard_function, reach_set.state_bounding_box());
        // For IMPACT events, also look at flow direction
        UpperIntervalType flow_derivative_range = apply(flow_derivative, reach_set.state_bounding_box());
        LOGGING_PRINTLN_AT(1, "guard_range=" << guard_range << ", flow_derivative_range=" << flow_derivative_range);
        if (possibly(guard_range.upper_bound() >= hybrid_zero)) {
            // Now make a set containing the complement of the constraint,
            // and test for emptiness. If the set is empty, then the guard is
            // not satisfied anywhere.
            HybridEnclosure test_set = reach_set;
            test_set.new_activation(event, guard_function);
            if (!definitely(test_set.is_empty())) {
                active_events.insert(*event_iter);
            }
        }
    }
    return active_events;
}

Map<DiscreteEvent, CrossingData>
    HybridContinuousTimeEvolverBase<ValidatedTag>::_compute_crossings(
        const Set<DiscreteEvent>& active_events, const EffectiveVectorMultivariateFunction& dynamic,
        const Map<DiscreteEvent, EffectiveScalarMultivariateFunction>& guards,
        const FlowFunctionModel& flow, const HybridEnclosure& initial_set) const {
    return EvolutionDetail::analyze_crossings(active_events, dynamic, guards, flow, initial_set, *this->_solver_ptr);
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::
    _recondition(HybridEnclosure& set) const {
    set.recondition();
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::
    _apply_reach_step(HybridEnclosure& set,
                      const ValidatedVectorMultivariateFunctionPatch& flow,
                      const TimingData& timing_data) const {
    set.apply_parameter_reach_step(flow, timing_data.parameter_dependent_evolution_time);
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::
    _apply_evolve_step(HybridEnclosure& set,
                       const ValidatedVectorMultivariateFunctionPatch& flow,
                       const TimingData& timing_data) const {
    switch (timing_data.step_kind) {
        case StepKind::CONSTANT_EVOLUTION_TIME:
        case StepKind::SPACE_DEPENDENT_EVOLUTION_TIME:
        case StepKind::TIME_DEPENDENT_EVOLUTION_TIME:
        case StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME:
        case StepKind::PARAMETER_DEPENDENT_EVOLUTION_TIME:
            set.apply_parameter_evolve_step(flow, timing_data.parameter_dependent_evolution_time);
            break;
        case StepKind::PARAMETER_DEPENDENT_FINISHING_TIME:
        case StepKind::SPACETIME_DEPENDENT_FINISHING_TIME:
        case StepKind::CONSTANT_FINISHING_TIME:
            set.apply_finishing_parameter_evolve_step(flow, timing_data.parameter_dependent_finishing_time);
            break;
        default:
            ARIADNE_FAIL_MSG("Unhandled step kind " << timing_data.step_kind);
    }
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::
    _process_jumped_set(const HybridEnclosure& jumped_set,
                        const TerminationType& termination,
                        Semantics semantics,
                        HybridEvolutionStepResult& result) const {
    const Real final_time = termination.maximum_time();
    const Integer maximum_steps = termination.maximum_steps();
    const Set<DiscreteEvent>& terminating_events = termination.terminating_events();

    HybridEnclosure initial_set = jumped_set;
    LOGGING_PRINTLN_AT(1, "initial_set=" << initial_set);

    LOGGING_RUN_MUTED(auto time_range = initial_set.time_range());
    if (definitely(time_range.lower_bound() > final_time)) {
        ARIADNE_WARN("initial_set.time_range()=" << time_range << " which exceeds final time=" << final_time);
        result.stop_reason = HybridEvolutionStopReason::PAST_HORIZON;
        return;
    }

    // Extract starting location
    const DiscreteLocation location = initial_set.location();

    // Cache dynamic and constraint functions
    LOGGING_RUN_AT(1, auto transitions = this->_extract_transitions(location));
    Set<DiscreteEvent> events = transitions.keys();

    LOGGING_PRINTLN_AT(1, "initial_set " << initial_set);
    LOGGING_PRINTLN_AT(2, "transitions=" << transitions);

    if (definitely(time_range.lower_bound() >= final_time)) {
        // The closed horizon includes points at t=T, including a reset that
        // lands exactly there. Retain its invariant-valid slice before any
        // further continuous or discrete evolution.
        this->_apply_invariants(initial_set, transitions);
        initial_set.set_time(final_time);
        result.stop_reason = HybridEvolutionStopReason::TIME_LIMIT;
        if (!definitely(initial_set.is_empty())) {
            result.record_final(initial_set);
        }
        return;
    }

    // Test if maximum number of steps has been exceeded; if so, the set should be discarded.
    // NOTE: We could also place a test for the maximum number of steps being reaches which computing jump sets
    // This is not done since the maximum_steps information is not passed to the _apply_evolution_step(...) method.
    if (initial_set.previous_events().size() >= maximum_steps) {
        LOGGING_PRINTLN_AT(2, "initial_set has undergone more than the maximum number of events (" << maximum_steps << ")");
        LOGGING_RUN_AT(1, this->_apply_invariants(initial_set, transitions));
        result.record_final(initial_set);
        result.stop_reason = HybridEvolutionStopReason::DISCRETE_LIMIT;
        return;
    }

    // Test if a terminating event has been reached.
    if (initial_set.previous_events().size() >= 1 && terminating_events.contains(initial_set.previous_events().back())) {
        LOGGING_PRINTLN_AT(2, "initial_set has undergone event " << initial_set.previous_events().back());
        LOGGING_RUN_AT(1, this->_apply_invariants(initial_set, transitions));
        result.record_final(initial_set);
        result.stop_reason = HybridEvolutionStopReason::TERMINATING_EVENT;
        return;
    }

    // Process the initially active events; cut out active points to leave initial flowable set.
    LOGGING_PRINTLN_AT(2, "initial_set has not reached the discrete termination condition");

    this->_process_starting_events(initial_set, transitions, semantics, result);
}

HybridEvolutionStepResult HybridContinuousTimeEvolverBase<ValidatedTag>::
    _process_working_set(const HybridEvolutionWorkItem& working_set,
                         const HybridTerminationCriterion& termination,
                         const Semantics& semantics,
                         const HybridEvolutionStatistics& statistics) const {
    LOGGING_SCOPE_CREATE;
    this->_validate_mode(working_set.enclosure.location(), semantics);
    HybridEvolutionStepResult result;
    switch (working_set.phase) {
        case HybridEvolutionPhase::PROCESS_STARTING_EVENTS:
            this->_process_jumped_set(working_set.enclosure, termination, semantics, result);
            break;
        case HybridEvolutionPhase::EVOLVE_CONTINUOUSLY:
            this->_log_summary(working_set.enclosure, statistics);
            this->_evolution_step(working_set.enclosure, termination.maximum_time(), semantics, result);
            break;
        default:
            ARIADNE_FAIL_MSG("Unknown hybrid evolution processing phase");
    }
    if (result.successors.empty() && result.stop_reason == HybridEvolutionStopReason::NONE) {
        result.stop_reason = HybridEvolutionStopReason::NO_CONTINUATION;
    }
    return result;
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::
    _evolution_step(const HybridEnclosure& set,
                    const Real& final_time,
                    const Semantics& semantics,
                    HybridEvolutionStepResult& result) const {
    HybridEnclosure starting_set = set;

    const DiscreteLocation& location = starting_set.location();
    const auto* inclusion = this->system().inclusion(location);
    this->_validate_mode(location, semantics);
    LOGGING_RUN_AT(1, auto transitions = this->_extract_transitions(location));

    LOGGING_PRINTLN_AT(1, "starting_set=" << starting_set);
    LOGGING_PRINTLN_AT(2, "starting_time=" << starting_set.time_function());
    if (definitely(starting_set.is_empty())) {
        LOGGING_PRINTLN_AT(1, "Empty starting_set " << starting_set);
        result.stop_reason = HybridEvolutionStopReason::EMPTY;
        return;
    }

    if (definitely(starting_set.time_range().lower_bound() >= final_time)) {
        ARIADNE_WARN("starting_set.time_range()=" << starting_set.time_range() << " which exceeds final time=" << final_time);
        result.stop_reason = HybridEvolutionStopReason::PAST_HORIZON;
        return;
    }

    // Compute the bounding box of the enclosure
    const ExactBoxType starting_bounding_box = cast_exact_box(starting_set.state_bounding_box());
    LOGGING_PRINTLN_AT(2, "starting_bounding_box=" << starting_bounding_box);

    const SizeType protected_parameters = starting_set.continuous_set().number_of_protected_parameters();
    // Allow the error axes introduced by reconditioning state, absolute time
    // and dwell time, so a structural step can actually relieve the pressure.
    const SizeType parameter_limit =
        (starting_set.continuous_set().configuration().reconditioning_num_blocks() + 1u) *
            starting_set.state_dimension() +
        protected_parameters + 2u;
    const bool parameter_pressure = inclusion != nullptr && starting_set.number_of_parameters() > parameter_limit;

    // Test to see if set requires reconditioning
    if (this->_configuration_ptr->enable_reconditioning() &&
        (parameter_pressure || possibly(norm(starting_set.state_function().errors()) > this->_configuration_ptr->maximum_spacial_error()) || possibly(starting_set.time_function().error() > this->_configuration_ptr->maximum_spacial_error()) || possibly(starting_set.continuous_set().dwell_time_function().error() > this->_configuration_ptr->maximum_spacial_error()))) {
        HybridEnclosure reconditioned_set = starting_set;
        this->_recondition(reconditioned_set);
        LOGGING_PRINTLN("RECONDITION location=" << location
                                                << " t=" << starting_set.time_range()
                                                << " parameters=" << starting_set.number_of_parameters()
                                                << "->" << reconditioned_set.number_of_parameters()
                                                << " radius=" << starting_set.continuous_set().radius()
                                                << "->" << reconditioned_set.continuous_set().radius()
                                                << " state_error=" << sup_norm(starting_set.state_function().errors())
                                                << "->" << sup_norm(reconditioned_set.state_function().errors())
                                                << " time_error=" << starting_set.time_function().error()
                                                << "->" << reconditioned_set.time_function().error());
        result.continue_with(reconditioned_set, HybridEvolutionSuccessorKind::RECONDITIONING);
        return;
    }

    // Handle a set that is too large, based on semantics
    if (possibly(starting_bounding_box.radius() > this->_configuration_ptr->maximum_enclosure_radius())) {
        if (semantics == Semantics::LOWER) {
            LOGGING_PRINTLN_AT(1, "Set too large, discarding");
            result.stop_reason = HybridEvolutionStopReason::LOWER_SIZE_LIMIT;
            return;
        } else if (this->_configuration_ptr->enable_subdivisions()) {
            List<HybridEnclosure> split_sets = starting_set.split();
            for (SizeType i = 0; i != split_sets.size(); ++i) {
                if (!definitely(split_sets[i].is_empty())) {
                    result.continue_with(split_sets[i], HybridEvolutionSuccessorKind::SUBDIVISION);
                }
            }
            LOGGING_PRINTLN("SPLIT location=" << location
                                              << " t=" << starting_set.time_range()
                                              << " pieces=" << result.successors.size()
                                              << " radius=" << starting_bounding_box.radius()
                                              << " limit=" << this->_configuration_ptr->maximum_enclosure_radius()
                                              << " parameters=" << starting_set.number_of_parameters());
            return;
        }
    }

    if (inclusion != nullptr) {
        this->_evolve_inclusion(starting_set, final_time, result);
        return;
    }
    if (this->system().time_variable(location) || EvolutionDetail::hybrid_parameter_space(this->system(), location).dimension() != 0u) {
        this->_evolve_outer_equation(starting_set, final_time, result);
        return;
    }
    const auto dynamic = this->system().dynamic_function(location);

    Map<DiscreteEvent, EffectiveScalarMultivariateFunction> guard_functions;
    for (Map<DiscreteEvent, TransitionData>::ConstIterator transition_iter = transitions.begin();
         transition_iter != transitions.end(); ++transition_iter) {
        guard_functions.insert(transition_iter->first, transition_iter->second.guard_function);
    }
    LOGGING_PRINTLN_AT(1, "guards=" << guard_functions);

    // Compute flow and actual time step size used
    // The shared bridge composes a patch over the entire enclosure parameter
    // domain. Certify that full state image, including points later excluded
    // by guards, rather than extrapolating a patch outside its certified box.
    const auto integration_domain = hull(starting_bounding_box,
                                         cast_exact_box(UpperBoxType(starting_set.state_function().range())));
    const auto step = this->_compute_integration_step(dynamic, integration_domain, this->configuration().maximum_step_size());
    // Crossing times and Lie derivatives currently refer to a single local
    // autonomous state flow. Do not discard patches or interpret input axes
    // as state coordinates when that capability is unavailable.
    const SizeType dimension = starting_set.state_function().result_size();
    if (step.time_step() <= StepSizeType(0u) || step.time_step() > this->configuration().maximum_step_size() ||
        step.reach_models().size() != 1u || !step.has_evolve_model() ||
        step.has_explicit_time() || step.constant_parameter_domain().dimension() != 0u ||
        step.initial_time() != StepSizeType(0u) || step.reach_time_offset() != StepSizeType(0u) ||
        step.reach_time_index() != std::optional<SizeType>(dimension)) {
        throw std::invalid_argument("Hybrid event analysis requires one structured autonomous local reach patch and endpoint");
    }
    const FlowFunctionModel flow_model(step.reach_models()[0]);
    if (flow_model.argument_size() != dimension + 1u || flow_model.result_size() != dimension ||
        flow_model.time_domain() != ExactIntervalType(0u, step.time_step()) ||
        step.evolve_model().argument_size() != dimension ||
        !subset(integration_domain, flow_model.space_domain()) ||
        !subset(integration_domain, step.evolve_model().domain())) {
        throw std::invalid_argument("Hybrid event analysis cannot interpret the integration step's coordinate layout");
    }
    LOGGING_PRINTLN_AT(1, "flow_model.domain()=" << flow_model.domain() << " flow_model.range()=" << flow_model.range());

    // Compute possibly active events
    Set<DiscreteEvent> active_events =
        this->_compute_active_events(dynamic, guard_functions, flow_model, starting_set);
    if (!active_events.empty()) {
        Map<DiscreteEvent, EventKind> active_events_log_data;
        for (Set<DiscreteEvent>::ConstIterator event_iter = active_events.begin(); event_iter != active_events.end(); ++event_iter) {
            active_events_log_data.insert(*event_iter, transitions[*event_iter].event_kind);
        }
        LOGGING_PRINTLN_AT(1, "active_events: " << active_events_log_data);
    }

    // Compute the kind of crossing (increasing, convex, etc);
    Map<DiscreteEvent, CrossingData> crossings =
        this->_compute_crossings(active_events, dynamic, guard_functions, flow_model, starting_set);
    // Crossing computation is read-only; prune conservative candidates here,
    // after iteration has finished, using the events that can actually cross.
    active_events = crossings.keys();
    LOGGING_PRINTLN_AT(1, "crossings=" << crossings);

    // Compute end conditions for flow
    TimingData timing_data = this->_estimate_timing(active_events, Real(final_time), flow_model, crossings, transitions, starting_set);
    LOGGING_PRINTLN_AT(1, "timing_data: " << timing_data);

    // Apply the time step
    HybridEnclosure reach_set, evolve_set;
    this->_apply_evolution_step(starting_set, step, flow_model, timing_data, crossings, dynamic, transitions, semantics, result);
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::_evolve_inclusion(
    const HybridEnclosure& initial, const Real& final_time, HybridEvolutionStepResult& result) const {
    const auto& inclusion = *this->system().inclusion(initial.location());
    // Inclusion candidates use only validated outer reach and state guards.
    // No root-time map, Lie derivative or full-source coverage is inferred.
    const auto invariants = this->_extract_transitions(initial.location());
    const auto horizon = FloatDPBounds(final_time, dp);
    const StepSizeType remaining = StepSizeType(horizon.upper_raw()) -
                                   StepSizeType(initial.time_range().lower_bound().raw());
    const auto maximum_step = EvolutionDetail::limit_step_to_remaining_time(
        this->configuration().maximum_step_size(), remaining);
    const auto& source = initial.continuous_set();
    const auto parameters = source.constant_parameter_function(inclusion.parameter_space().dimension());
    // The field is autonomous. The bridge advances each source time/dwell
    // map, including a parameter-dependent arrival time after a reset.
    const auto numerical = [&]() -> EvolutionDetail::AnchoredIntegrationStep {
        if (inclusion.has_explicit_time()) {
            return this->_inclusion_integrator_ptr
                     ? EvolutionDetail::compute_clocked_inclusion_step(inclusion, source, parameters, maximum_step, *this->_inclusion_integrator_ptr)
                     : EvolutionDetail::compute_clocked_inclusion_step(inclusion, source, parameters, maximum_step, *this->_box_integrator_ptr);
        }
        const auto step = this->_inclusion_integrator_ptr
                            ? EvolutionDetail::compute_inclusion_step(inclusion, source, parameters, StepSizeType(0u), maximum_step, *this->_inclusion_integrator_ptr)
                            : EvolutionDetail::compute_inclusion_step(inclusion, source, parameters, StepSizeType(0u), maximum_step, *this->_box_integrator_ptr);
        return {step, join(source.state_function(), parameters)};
    }();
    const auto& step = numerical.step;
    const auto applied = EvolutionDetail::apply_integration_step(source, step, numerical.arguments, Semantics::UPPER,
                                                                 this->_box_integrator_ptr ? EvolutionDetail::IntegrationStepRepresentation::ALLOW_BOX_ENCLOSURES
                                                                                           : EvolutionDetail::IntegrationStepRepresentation::REQUIRE_STRUCTURED);
    result = EvolutionDetail::publish_outer_step(initial, applied, invariants, final_time,
                                                 std::min(step.time_step().get_d(), remaining.get_d()));
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::_evolve_outer_equation(
    const HybridEnclosure& initial, const Real& final_time, HybridEvolutionStepResult& result) const {
    const auto remaining = StepSizeType(FloatDPBounds(final_time, dp).upper_raw()) -
                           StepSizeType(initial.time_range().lower_bound().raw());
    const auto maximum_step = EvolutionDetail::limit_step_to_remaining_time(configuration().maximum_step_size(), remaining);
    const auto numerical = EvolutionDetail::compute_clocked_equation_step(
        EvolutionDetail::prepare_equation_field(system().dynamic_function(initial.location()), initial.state_dimension(),
                                                system().time_variable(initial.location()).has_value(), EvolutionDetail::hybrid_parameter_space(system(), initial.location()).dimension()),
        initial.continuous_set(),
        initial.continuous_set().constant_parameter_function(EvolutionDetail::hybrid_parameter_space(system(), initial.location()).dimension()), maximum_step, *this->_integrator_ptr);
    const auto applied = EvolutionDetail::apply_integration_step(initial.continuous_set(), numerical.step,
                                                                 numerical.arguments, Semantics::UPPER);
    result = EvolutionDetail::publish_outer_step(initial, applied, _extract_transitions(initial.location()), final_time,
                                                 std::min(numerical.step.time_step().get_d(), remaining.get_d()));
}

Void HybridContinuousTimeEvolverBase<ValidatedTag>::
    _apply_evolution_step(const HybridEnclosure& starting_set,
                          const ValidatedIntegrationStep& step,
                          const ValidatedVectorMultivariateFunctionPatch& flow,
                          const TimingData& timing_data,
                          const Map<DiscreteEvent, CrossingData>& crossings,
                          const EffectiveVectorMultivariateFunction& dynamic,
                          const Map<DiscreteEvent, TransitionData>& transitions,
                          const Semantics& semantics,
                          HybridEvolutionStepResult& result) const {
    LOGGING_SCOPE_CREATE;

    // The sole caller, _evolution_step, has already rejected an empty set.
    // Flow/crossing/timing estimation does not mutate starting_set.
    EvolutionStepData _step_data;

    // Compute events enabling transitions and events blocking continuous evolution
    Set<DiscreteEvent> active_events = crossings.keys();
    Set<DiscreteEvent> activating_events = intersection(Ariadne::activating_events(transitions), active_events);
    Set<DiscreteEvent> blocking_events = intersection(Ariadne::blocking_events(transitions), active_events);

    LOGGING_PRINTLN("activating_events=" << activating_events);
    LOGGING_PRINTLN("blocking_events=" << blocking_events);

    const ValidatedScalarMultivariateFunctionPatch evolve_step_time = timing_data.parameter_dependent_evolution_time;
    const ValidatedScalarMultivariateFunctionPatch reach_step_time = embed(starting_set.parameter_domain(), timing_data.evolution_time_coordinate);

    LOGGING_PRINTLN_AT(1, "evolve_step_time=" << evolve_step_time << "\n")
    LOGGING_PRINTLN_AT(1, "reach_step_time=" << reach_step_time << "\n")

    // Compute the reach and evolve sets, without introducing bounds due to the final time.
    List<HybridEnclosure> reach_sets = {starting_set};
    List<HybridEnclosure> evolve_sets = {starting_set};

    const auto elapsed = timing_data.parameter_dependent_evolution_time.range();
    const FloatDP full_duration(step.time_step(), dp);
    if (semantics == Semantics::UPPER &&
        elapsed.lower_bound().raw() == full_duration && elapsed.upper_bound().raw() == full_duration) {
        // At the certified endpoint use its independent map, rather than a
        // slice of a possibly wider tube. Preserve hybrid labels and history.
        const auto applied = EvolutionDetail::apply_integration_step(
            starting_set.continuous_set(), step, EvolutionDetail::IntegrationStepCoordinates::STATE, semantics);
        reach_sets.front().replace_continuous_set(applied.reach.front());
        evolve_sets.front().replace_continuous_set(applied.evolved);
        if (timing_data.finishing_kind == FinishingKind::AT_FINAL_TIME) {
            evolve_sets.front().set_time(timing_data.final_time);
        }
    } else {
        // Event-limited or parameter-dependent stopping times lie within the
        // tube; the endpoint certificate only covers the complete duration.
        _apply_reach_step(reach_sets.front(), flow, timing_data);
        _apply_evolve_step(evolve_sets.front(), flow, timing_data);
    }

    // This one continuous computation can publish several outputs or successors.
    // Credit its elapsed flow span once; structural steps never reach this path.
    const double flow_duration = timing_data.parameter_dependent_evolution_time.range().upper_bound().raw().get_d();
    const double remaining_duration = timing_data.final_time.get_d() - starting_set.time_range().lower_bound().get_d();
    result.continuous_evolution_time = std::max(0.0, std::min(flow_duration, remaining_duration));

    LOGGING_PRINTLN_AT(1, "flow_reach_set=" << reach_sets.front() << "\n")
    LOGGING_PRINTLN_AT(1, "flow_evolve_set=" << evolve_sets.front() << "\n")

    // Apply constraints on reach and evolve sets due to invariants and urgent guards
    for (Set<DiscreteEvent>::ConstIterator event_iter = blocking_events.begin();
         event_iter != blocking_events.end(); ++event_iter) {
        const TransitionData& transition_data = transitions[*event_iter];
        const CrossingData& crossing_data = crossings[*event_iter];

        const EvolutionDetail::HybridGuardContext context{
            starting_set, flow, transition_data, crossing_data, semantics, this->function_factory()};
        EvolutionDetail::apply_flow_guard(reach_sets, context, reach_step_time);
        if (!EvolutionDetail::defer_endpoint_guard(crossing_data)) {
            EvolutionDetail::apply_flow_guard(evolve_sets, context, evolve_step_time);
        }
    }

    // Every endpoint is coupled to its pre-deferred continuation checkpoint.
    // Deferred constraints modify only the reported component.
    auto endpoint_branches = EvolutionDetail::pair_endpoint_branches(evolve_sets);
    for (Set<DiscreteEvent>::ConstIterator event_iter = blocking_events.begin();
         event_iter != blocking_events.end(); ++event_iter) {
        const TransitionData& transition_data = transitions[*event_iter];
        const CrossingData& crossing_data = crossings[*event_iter];
        if (EvolutionDetail::defer_endpoint_guard(crossing_data)) {
            const EvolutionDetail::HybridGuardContext context{
                starting_set, flow, transition_data, crossing_data, semantics, this->function_factory()};
            EvolutionDetail::restrict_reported_endpoints(endpoint_branches, context, evolve_step_time);
        }
    }

    // Publish continuous contributions before transitions and their resets.
    const auto publication = EvolutionDetail::publish_step_contributions(
        reach_sets, endpoint_branches, timing_data, result);
    _step_data.progress = publication.progress;
    _step_data.finishing = publication.finishing;

    LOGGING_PRINTLN("Computing jump sets");

    // Compute jump sets
    for (Set<DiscreteEvent>::ConstIterator event_iter = activating_events.begin();
         event_iter != activating_events.end(); ++event_iter) {
        DiscreteEvent event = *event_iter;

        // Compute active set
        const EvolutionDetail::HybridGuardContext context{
            starting_set, flow, transitions[event], crossings[event], semantics, this->function_factory()};
        auto candidate = EvolutionDetail::make_jump_candidate(context, dynamic, timing_data);
        List<HybridEnclosure> jump_sets{std::move(candidate.enclosure)};
        const auto jump_step_time = std::move(candidate.elapsed_time);

        // Apply maximum time bound, as this will be applied after the next flow step
        for (List<HybridEnclosure>::Iterator jump_set_iter = jump_sets.begin(); jump_set_iter != jump_sets.end(); ++jump_set_iter) {
            HybridEnclosure& _jump_set = *jump_set_iter;
            if (possibly(_jump_set.time_range().upper_bound() > timing_data.final_time)) {
                _jump_set.bound_time(timing_data.final_time);
                if (definitely(not _jump_set.is_empty())) {
                    ARIADNE_WARN("Explicitly bounding time in jump set");
                }
            }
        }

        LOGGING_PRINTLN_AT(1, "Applying guard to jump sets");
        // Apply blocking conditions for other active events
        for (Set<DiscreteEvent>::ConstIterator other_event_iter = blocking_events.begin();
             other_event_iter != blocking_events.end(); ++other_event_iter) {
            LOGGING_PRINTLN_AT(2, "jump_sets.size()=" << jump_sets.size());
            DiscreteEvent other_event = *other_event_iter;
            if (other_event != event) {
                const TransitionData& other_transition_data = transitions[other_event];
                const CrossingData& other_crossing_data = crossings[other_event];
                const EvolutionDetail::HybridGuardContext other_context{
                    starting_set, flow, other_transition_data, other_crossing_data, semantics, this->function_factory()};
                EvolutionDetail::apply_flow_guard(jump_sets, other_context, jump_step_time);
            }
        }

        LOGGING_PRINTLN_AT(1, event << ": " << transitions[event].event_kind << ", " << crossings[event].summary());
        // Apply reset
        for (List<HybridEnclosure>::Iterator jump_set_iter = jump_sets.begin(); jump_set_iter != jump_sets.end(); ++jump_set_iter) {
            HybridEnclosure& _jump_set = *jump_set_iter;
            if (!definitely(_jump_set.is_empty())) {
                const DiscreteLocation& target = transitions[event].target;
                LOGGING_PRINTLN_AT(2, "target=" << target << ", auxiliary_space=" << this->system().continuous_auxiliary_space(target) << ", auxiliary_function=" << this->system().auxiliary_function(target));
                const TransitionData& transition = transitions[event];
                EvolutionDetail::apply_transition_reset(_jump_set, transition);
                result.transition_to(_jump_set, event, transition.event_kind, crossings[event].summary());
                _step_data.events.insert(event);
                LOGGING_PRINTLN_AT(1, "jump_set=" << _jump_set);
            }
        }
    }

    if (result.successors.empty() && _step_data.finishing) {
        result.stop_reason = HybridEvolutionStopReason::TIME_LIMIT;
    }
    if (_step_data.finishing || !_step_data.progress || !_step_data.events.empty()) {
        LOGGING_PRINTLN_AT(2, timing_data.finishing_kind << " " << _step_data.events << (_step_data.progress ? " progress" : " no progress"));
    }
}

TimingData
    HybridContinuousTimeEvolverBase<ValidatedTag>::
        _estimate_timing(Set<DiscreteEvent>&,
                         Real final_time,
                         const FlowFunctionModel& flow,
                         Map<DiscreteEvent, CrossingData>&,
                         const Map<DiscreteEvent, TransitionData>&,
                         const HybridEnclosure& initial_set) const {
    // Compute the evolution time for the given step.
    LOGGING_SCOPE_CREATE;
    const StepSizeType step_size = static_cast<StepSizeType>(flow.domain()[flow.domain().size() - 1].upper_bound());
    TimingData result;
    result.step_kind = StepKind::CONSTANT_EVOLUTION_TIME;
    result.finishing_kind = FinishingKind::STRADDLE_FINAL_TIME;
    result.step_size = step_size;
    result.final_time = final_time;
    result.evolution_time_domain = ExactIntervalType(0, step_size);
    result.evolution_time_coordinate = this->function_factory().create_coordinate({result.evolution_time_domain}, 0u);
    result.parameter_dependent_evolution_time = this->function_factory().create_constant(initial_set.parameter_domain(), result.step_size);
    LOGGING_PRINTLN("timing_data=" << result);
    return result;
}

HybridContinuousTimeEvolverBaseConfiguration<ValidatedTag>::HybridContinuousTimeEvolverBaseConfiguration(HybridContinuousTimeEvolverBase<ValidatedTag>& evolver)
    : _evolver(evolver) {
    set_flow_accuracy(1e-5);
    set_maximum_step_size(1.0);
    set_maximum_enclosure_radius(100.0);
    set_maximum_spacial_error(1e-2);
    set_enable_reconditioning(true);
    set_enable_subdivisions(true);
}

Void HybridContinuousTimeEvolverBaseConfiguration<ValidatedTag>::set_flow_accuracy(const ApproximateRealType value) {
    _evolver._integrator_ptr = std::shared_ptr<GradedTaylorSeriesIntegrator>(new GradedTaylorSeriesIntegrator(StepMaximumError(cast_exact(value))));
    _flow_accuracy = cast_exact(value);
}

OutputStream&
    HybridContinuousTimeEvolverBaseConfiguration<ValidatedTag>::_write(OutputStream& os) const {
    os << "HybridContinuousTimeEvolverBaseConfiguration<ValidatedTag>"
       << "(\n  flow_accuracy=" << flow_accuracy()
       << ",\n  maximum_step_size=" << maximum_step_size()
       << ",\n  maximum_enclosure_radius=" << maximum_enclosure_radius()
       << ",\n  maximum_spacial_error=" << maximum_spacial_error()
       << ",\n  enable_reconditioning=" << enable_reconditioning()
       << ",\n  enable_subdivisions=" << enable_subdivisions()
       << "\n)";
    return os;
}

HybridContinuousTimeEvolver<ValidatedTag>::HybridContinuousTimeEvolver(const SystemType& system)
    : HybridContinuousTimeEvolverBase<ValidatedTag>(system) {
    this->_configuration_ptr.reset(new HybridContinuousTimeEvolverConfiguration<ValidatedTag>(*this));
}

HybridContinuousTimeEvolver<ValidatedTag>::HybridContinuousTimeEvolver(
    const SystemType& system,
    const FunctionFactoryType& factory)
    : HybridContinuousTimeEvolverBase<ValidatedTag>(system, factory) {
    this->_configuration_ptr.reset(new HybridContinuousTimeEvolverConfiguration<ValidatedTag>(*this));
}

HybridContinuousTimeEvolver<ValidatedTag>::HybridContinuousTimeEvolver(const HybridContinuousTimeEvolver<ValidatedTag>& other)
    : HybridContinuousTimeEvolverBase<ValidatedTag>(other.system(), other.function_factory()) {
    this->_copy_configuration_from(other);
}

HybridContinuousTimeEvolver<ValidatedTag>& HybridContinuousTimeEvolver<ValidatedTag>::operator=(const HybridContinuousTimeEvolver<ValidatedTag>& other) {
    if (this != &other) {
        this->_sys_ptr.reset(other.system().clone());
        this->set_function_factory(other.function_factory());
        this->_copy_configuration_from(other);
    }
    return *this;
}

Void HybridContinuousTimeEvolver<ValidatedTag>::_copy_configuration_from(const HybridContinuousTimeEvolver<ValidatedTag>& other) {
    // The configuration stores its owning evolver by reference. Rebuild that
    // relationship so copies remain independent and outlive the original.
    auto integrator = std::shared_ptr<IntegratorInterface>(other._integrator_ptr->clone());
    auto solver = std::shared_ptr<EquationSolverInterface>(other._solver_ptr->clone());
    this->_configuration_ptr.reset(new HybridContinuousTimeEvolverConfiguration<ValidatedTag>(*this));
    const auto& source = other.configuration();
    auto& target = this->configuration();
    target.set_flow_accuracy(source.flow_accuracy());
    target.set_maximum_step_size(source.maximum_step_size());
    target.set_maximum_enclosure_radius(source.maximum_enclosure_radius());
    target.set_maximum_spacial_error(source.maximum_spacial_error());
    target.set_enable_reconditioning(source.enable_reconditioning());
    target.set_enable_subdivisions(source.enable_subdivisions());
    this->_integrator_ptr = std::move(integrator);
    this->_inclusion_integrator_ptr = other._inclusion_integrator_ptr;
    this->_box_integrator_ptr = other._box_integrator_ptr;
    this->_solver_ptr = std::move(solver);
    this->ALLOW_CREEP = other.ALLOW_CREEP;
    this->ALLOW_UNWIND = other.ALLOW_UNWIND;
}

TimingData
    HybridContinuousTimeEvolver<ValidatedTag>::
        _estimate_timing(Set<DiscreteEvent>&,
                         Real final_time,
                         const FlowFunctionModel& flow,
                         Map<DiscreteEvent, CrossingData>& crossings,
                         const Map<DiscreteEvent, TransitionData>& transitions,
                         const HybridEnclosure& initial_set) const {
    return EvolutionDetail::estimate_hybrid_timing(
        final_time, flow, crossings, transitions, initial_set, this->system(),
        this->function_factory(), *this->_solver_ptr, ALLOW_CREEP, ALLOW_UNWIND);
}

HybridContinuousTimeEvolverConfiguration<ValidatedTag>::HybridContinuousTimeEvolverConfiguration(HybridContinuousTimeEvolver<ValidatedTag>& evolver)
    : HybridContinuousTimeEvolverBaseConfiguration<ValidatedTag>(evolver) {
}

HybridContinuousTimeEvolverFactory<ValidatedTag>::HybridContinuousTimeEvolverFactory()
    : HybridContinuousTimeEvolverFactory<ValidatedTag>(ValidatedFunctionPatchFactory(make_taylor_function_patch_factory())) {
}

HybridContinuousTimeEvolverFactory<ValidatedTag>::HybridContinuousTimeEvolverFactory(const ValidatedFunctionPatchFactory& factory)
    : _function_factory_ptr(factory.managed_pointer()) {
}

HybridContinuousTimeEvolver<ValidatedTag>*
    HybridContinuousTimeEvolverFactory<ValidatedTag>::create(const ContinuousTimeHybridAutomatonInterface& system) const {
    return new HybridContinuousTimeEvolver<ValidatedTag>(system, ValidatedFunctionPatchFactory(_function_factory_ptr));
}

} // namespace Ariadne
