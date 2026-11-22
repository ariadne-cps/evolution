/***************************************************************************
 *            continuous_time_evolver.cpp
 *
 *  Copyright  2008-20  Alberto Casagrande, Pieter Collins
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

/** @file src/continuous_time_evolver.cpp
 *  @brief Implement vector field evolver operations.
 */

#include "evolution/continuous_time_evolver.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>

#include "algebra/vector.hpp"
#include "continuous_integration.hpp"
#include "evolution/enclosure.hpp"
#include "evolution/orbit.hpp"
#include "evolution_time.hpp"
#include "function/constraint.hpp"
#include "function/function.hpp"
#include "function/functional.hpp"
#include "function/taylor_function.hpp"
#include "integration/integrator.hpp"
#include "integration/validated_step_integrator.hpp"
#include "logging/logging.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/space.hpp"
#include "system/odae.hpp"
#include "utility/array.hpp"
#include "utility/container.hpp"
#include "utility/macros.hpp"
#include "utility/stlio.hpp"
#include "utility/tuple.hpp"

#include "expression-set/expression_set.hpp"

namespace Ariadne {

// TODO: Move to Numeric module
inline PositiveValidatedUpperNumber abs(PositiveValidatedUpperNumber y) { return y; }

namespace {

EffectiveVectorMultivariateFunction equation_field(const OrdinaryDifferentialAlgebraicEquation& system) {
    if (system.has_explicit_time() || system.parameter_space().dimension() == 0u) {
        return system.function();
    }
    // The numerical timed-parameter overload expects [state,time,parameters].
    // An autonomous parameterised field ignores that argument; no clock or
    // parameter is added to the physical state or differentiated.
    const SizeType n = system.dimension();
    const SizeType np = system.parameter_space().dimension();
    const auto coordinates = EffectiveScalarMultivariateFunction::coordinates(n + 1u + np);
    List<EffectiveScalarMultivariateFunction> arguments;
    for (SizeType i = 0u; i != n; ++i) {
        arguments.append(coordinates[i]);
    }
    for (SizeType i = 0u; i != np; ++i) {
        arguments.append(coordinates[n + 1u + i]);
    }
    return compose(system.function(), EffectiveVectorMultivariateFunction(arguments));
}

ValidatedVectorMultivariateFunctionPatch constant_parameter_coordinates(const Enclosure& set, SizeType count) {
    Array<SizeType> indices(count);
    SizeType remaining = count;
    for (SizeType i = set.number_of_parameters(); i != 0u && remaining != 0u; --i) {
        if (set.variable_kinds()[i - 1u] == EnclosureVariableKind::PARAMETER) {
            indices[--remaining] = i - 1u;
        }
    }
    if (remaining != 0u) {
        throw std::logic_error("Continuous evolution lost its trajectory-constant parameter coordinates");
    }
    const auto& factory = set.configuration().function_factory();
    auto result = factory.create_zeros(count, set.parameter_domain());
    for (SizeType i = 0u; i != count; ++i) {
        result[i] = factory.create_coordinate(set.parameter_domain(), indices[i]);
    }
    return result;
}

template<class ES>
List<ES> subdivide(const ES& enclosure) {
    List<ES> result;
    Pair<ES, ES> split = enclosure.split();
    result.append(split.first);
    result.append(split.second);
    return result;
}

struct CarriedExpansionSnapshot {
    static constexpr std::size_t degree_slots = 32u;
    std::array<unsigned long long, degree_slots> below_tight_count{};
    std::array<unsigned long long, degree_slots> bridge_count{};
    std::array<unsigned long long, degree_slots> above_loose_count{};
    std::array<double, degree_slots> below_tight_abs_mass{};
    std::array<double, degree_slots> bridge_abs_mass{};
    std::array<double, degree_slots> above_loose_abs_mass{};
};

CarriedExpansionSnapshot
    carried_expansion_snapshot(const ValidatedVectorMultivariateTaylorFunctionModelDP& mapping) {
    CarriedExpansionSnapshot profile;
    for (SizeType component = 0u; component != mapping.result_size(); ++component) {
        const auto component_mapping = mapping[component];
        const auto& expansion = component_mapping.model().expansion();
        for (const auto& term: expansion) {
            auto bucket_degree = static_cast<std::size_t>(term.index().degree());
            if (bucket_degree >= CarriedExpansionSnapshot::degree_slots) {
                bucket_degree = CarriedExpansionSnapshot::degree_slots - 1u;
            }
            const double magnitude = std::abs(term.coefficient().get_d());
            if (magnitude < 3e-14) {
                ++profile.below_tight_count[bucket_degree];
                profile.below_tight_abs_mass[bucket_degree] += magnitude;
            } else if (magnitude < 1e-12) {
                ++profile.bridge_count[bucket_degree];
                profile.bridge_abs_mass[bucket_degree] += magnitude;
            } else {
                ++profile.above_loose_count[bucket_degree];
                profile.above_loose_abs_mass[bucket_degree] += magnitude;
            }
        }
    }
    return profile;
}

Void print_carried_expansion_snapshot(
    OutputStream& output, const String& stage,
    SizeType step,
    const ValidatedVectorMultivariateTaylorFunctionModelDP& mapping,
    const Sweeper<FloatDP>& sweeper) {
    const auto profile = carried_expansion_snapshot(mapping);
    for (std::size_t bucket_degree = 0u;
         bucket_degree != CarriedExpansionSnapshot::degree_slots; ++bucket_degree) {
        const auto total_count =
            profile.below_tight_count[bucket_degree] + profile.bridge_count[bucket_degree] + profile.above_loose_count[bucket_degree];
        if (total_count == 0u) {
            continue;
        }
        output << "[CarriedExpansionSnapshot]"
               << " sweeper=" << sweeper
               << " stage=" << stage
               << " step=" << step
               << " degree=" << bucket_degree
               << " below_3e-14_count=" << profile.below_tight_count[bucket_degree]
               << " bridge_3e-14_to_1e-12_count=" << profile.bridge_count[bucket_degree]
               << " above_1e-12_count=" << profile.above_loose_count[bucket_degree]
               << " below_3e-14_abs_mass=" << profile.below_tight_abs_mass[bucket_degree]
               << " bridge_3e-14_to_1e-12_abs_mass=" << profile.bridge_abs_mass[bucket_degree]
               << " above_1e-12_abs_mass=" << profile.above_loose_abs_mass[bucket_degree]
               << std::endl;
    }
}

} // namespace

namespace EvolutionDetail {

Dyadic limit_step_to_remaining_time(const Dyadic& configured, const Dyadic& remaining) {
    const FloatDP representable(remaining, FloatDP::ROUND_DOWNWARD, DoublePrecision());
    const Dyadic binary64_step(representable);
    if (binary64_step == Dyadic(0u)) {
        throw std::runtime_error(
            "ContinuousTimeEvolver: remaining time is too small for a positive binary64 step");
    }
    return min(configured, binary64_step);
}

} // namespace EvolutionDetail

ContinuousDynamics::ContinuousDynamics(const OrdinaryDifferentialAlgebraicEquation& system)
    : _declaration(SharedPointer<const OrdinaryDifferentialAlgebraicEquation>(system.clone())) {}
ContinuousDynamics::ContinuousDynamics(const OrdinaryDifferentialAlgebraicInclusion& system)
    : _declaration(SharedPointer<const OrdinaryDifferentialAlgebraicInclusion>(system.clone())) {}
ContinuousDynamics* ContinuousDynamics::clone() const { return new ContinuousDynamics(*this); }
SizeType ContinuousDynamics::dimension() const {
    return std::visit([](const auto& system) { return system->dimension(); }, _declaration);
}
RealSpace ContinuousDynamics::state_space() const {
    return std::visit([](const auto& system) -> RealSpace { return system->state_space(); }, _declaration);
}
RealSpace ContinuousDynamics::auxiliary_space() const {
    const auto* value = equation();
    return value != nullptr ? value->auxiliary_space() : RealSpace{};
}
RealSpace ContinuousDynamics::state_auxiliary_space() const {
    const auto* value = equation();
    return value != nullptr ? value->state_auxiliary_space() : state_space();
}
const RealSpace& ContinuousDynamics::parameter_space() const {
    return std::visit([](const auto& system) -> const RealSpace& { return system->parameter_space(); }, _declaration);
}
const BoxDomainType& ContinuousDynamics::constant_parameter_domain() const {
    return std::visit([](const auto& system) -> const BoxDomainType& { return system->constant_parameter_domain(); }, _declaration);
}
bool ContinuousDynamics::has_explicit_time() const {
    return std::visit([](const auto& system) { return system->has_explicit_time(); }, _declaration);
}
const OrdinaryDifferentialAlgebraicEquation* ContinuousDynamics::equation() const {
    const auto* value = std::get_if<SharedPointer<const OrdinaryDifferentialAlgebraicEquation>>(&_declaration);
    return value != nullptr ? value->get() : nullptr;
}
const OrdinaryDifferentialAlgebraicInclusion* ContinuousDynamics::inclusion() const {
    const auto* value = std::get_if<SharedPointer<const OrdinaryDifferentialAlgebraicInclusion>>(&_declaration);
    return value != nullptr ? value->get() : nullptr;
}

ContinuousTimeEvolver::ContinuousTimeEvolver(const OrdinaryDifferentialAlgebraicEquation& system, const IntegratorInterface& i)
    : _system(new SystemType(system)), _integrator(i.clone()),
      _function_factory(dynamic_cast<const IntegratorBase&>(i).function_factory()), _configuration(new ConfigurationType()) {
    if ((system.has_explicit_time() || system.parameter_space().dimension() != 0u) &&
        system.auxiliary_space().size() != 0u) {
        // Enclosure currently represents auxiliary outputs as state-only
        // functions. Time- or parameter-dependent outputs must not be evaluated
        // using a map with the wrong argument space.
        throw std::invalid_argument(
            "ContinuousTimeEvolver does not yet support auxiliary outputs with explicit time or constant parameters");
    }
    if ((system.has_explicit_time() || system.parameter_space().dimension() != 0u) &&
        dynamic_cast<const PreconditionedGradedTaylorSeriesIntegrator*>(&i) != nullptr) {
        throw std::invalid_argument("Preconditioned continuous evolution does not yet support explicit time or constant parameters");
    }
}

ContinuousTimeEvolver::ContinuousTimeEvolver(const OrdinaryDifferentialAlgebraicInclusion& system,
                                             const IntegratorBase& integrator,
                                             const InputApproximation& reach, const InputApproximation& evolve)
    : _system(new SystemType(system)),
      _inclusion_integrator(new SelectedParametricMultifunctionIntegrator(integrator, reach, evolve)),
      _function_factory(integrator.function_factory()), _configuration(new ConfigurationType()) {
    if (dynamic_cast<const PreconditionedGradedTaylorSeriesIntegrator*>(&integrator) != nullptr) {
        throw std::invalid_argument("Preconditioned continuous integration does not support inclusion input approximations");
    }
    if (system.has_input_function() && system.number_of_inputs() == 0u) {
        throw std::invalid_argument("The selected high-order inclusion method requires free inputs; select a capable solver explicitly");
    }
}

ContinuousTimeEvolver::ContinuousTimeEvolver(const OrdinaryDifferentialAlgebraicInclusion& system,
                                             const FirstOrderMultifunctionIntegrator& integrator)
    : _system(new SystemType(system)), _box_integrator(new FirstOrderMultifunctionIntegrator(integrator)),
      _function_factory(new TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-10))),
      _configuration(new ConfigurationType()) {}

ValidatedIntegrationStep ContinuousTimeEvolver::_inclusion_step(
    const EnclosureType& current, const TimeStepType& time, const StepSizeType& maximum_step) const {
    const auto& inclusion = *this->system().inclusion();
    const auto parameters = constant_parameter_coordinates(current, inclusion.parameter_space().dimension());
    if (this->_inclusion_integrator) {
        return EvolutionDetail::compute_inclusion_step(inclusion, current.euclidean_set(), parameters, time, maximum_step, *this->_inclusion_integrator);
    }
    return EvolutionDetail::compute_inclusion_step(inclusion, current.euclidean_set(), parameters, time, maximum_step, *this->_box_integrator);
}

auto ContinuousTimeEvolver::enclosure(const ExactBoxType& box) const -> EnclosureType {
    return EnclosureType(box, this->system().state_space(), EnclosureConfiguration(this->function_factory()));
}

auto ContinuousTimeEvolver::function_factory() const -> const FunctionFactoryType& {
    return this->_function_factory;
}

auto ContinuousTimeEvolver::orbit(const RealVariablesBox& initial_set, const TimeType& time, Semantics semantics) const -> Orbit<EnclosureType> {
    auto enclosure = EnclosureType(initial_set, this->system().state_space(), EnclosureConfiguration(this->function_factory()));
    // A time-dependent field without auxiliary outputs still carries a
    // zero-output auxiliary function with arguments [state,time]. Enclosure
    // already has a state-arity zero-output mapping, so leave it unchanged.
    if (this->system().auxiliary_space().size() != 0u) {
        enclosure.set_auxiliary(this->system().auxiliary_space(), this->system().equation()->auxiliary_function());
    }
    return orbit(enclosure, time, semantics);
}

auto ContinuousTimeEvolver::orbit(const RealExpressionBoundedConstraintSet& initial_set, const TimeType& time, Semantics semantics) const -> Orbit<EnclosureType> {
    auto enclosure = EnclosureType(initial_set.euclidean_set(this->system().state_space()), this->system().state_space(), EnclosureConfiguration(this->function_factory()));
    // A time-dependent field without auxiliary outputs still carries a
    // zero-output auxiliary function with arguments [state,time]. Enclosure
    // already has a state-arity zero-output mapping, so leave it unchanged.
    if (this->system().auxiliary_space().size() != 0u) {
        enclosure.set_auxiliary(this->system().auxiliary_space(), this->system().equation()->auxiliary_function());
    }
    return orbit(enclosure, time, semantics);
}

auto ContinuousTimeEvolver::orbit(const EnclosureType& initial_set, const TimeType& time, Semantics semantics) const -> Orbit<EnclosureType> {
    LOGGING_SCOPE_CREATE
    ARIADNE_PRECONDITION(this->system().state_auxiliary_space() == initial_set.state_auxiliary_space())

    if (this->system().inclusion() != nullptr && semantics != Semantics::UPPER) {
        throw std::invalid_argument("Continuous inclusion evolution requires UPPER semantics");
    }
    EnclosureType initial = initial_set;
    if (this->system().parameter_space().dimension() != 0u) {
        if (semantics == Semantics::LOWER) {
            throw std::invalid_argument("Parameterised continuous evolution currently requires UPPER semantics");
        }
        initial.adjoin_parameters(this->system().constant_parameter_domain(), EnclosureVariableKind::PARAMETER);
    }
    auto result = std::make_shared<SynchronisedOrbit>(initial);
    WorkloadType workload([time](const TimedEnclosureType& timed_enclosure, SharedPointer<ProgressIndicator> indicator) {
        indicator->update_current(timed_enclosure.first.get_d());
        indicator->update_final(time.get_d());
    },
                          std::bind_front(&ContinuousTimeEvolver::_process_timed_enclosure, this), time, semantics, result);
    _append_initial_set(workload, TimeStepType(0u), initial);
    workload.process();

    return std::move(*result);
}

Void ContinuousTimeEvolver::
    _append_initial_set(WorkloadType& workload, const TimeStepType& initial_time, const EnclosureType& current_set) const {
    if (possibly(current_set.euclidean_set().bounding_box().radius() > this->_configuration->maximum_enclosure_radius())) {
        LOGGING_PRINTLN_AT(1, "set is too large, splitting")
        Pair<EnclosureType, EnclosureType> split_sets = current_set.split();
        if (!definitely(split_sets.first.is_empty())) {
            _append_initial_set(workload, initial_time, split_sets.first);
        }
        if (!definitely(split_sets.second.is_empty())) {
            _append_initial_set(workload, initial_time, split_sets.second);
        }
    } else {
        workload.append({initial_time, current_set});
    }
}

Void ContinuousTimeEvolver::
    _process_timed_enclosure(WorkloadType::Access& workload,
                             const TimedEnclosureType& current_timed_set,
                             const TimeType& maximum_time,
                             Semantics semantics,
                             SharedPointer<SynchronisedOrbit> result) const {
    LOGGING_SCOPE_CREATE
    TimeStepType current_time = current_timed_set.first;
    EnclosureType current_set = current_timed_set.second;
    FloatDPUpperBound current_set_radius = current_set.euclidean_set().bounding_box().radius();

    LOGGING_PRINTLN("#r=" << std::setw(5) << std::left << result->reach_size()
                          << " t=" << std::setw(7) << std::fixed << current_time.get_d()
                          << " p=" << std::setw(4) << std::left << current_set.number_of_parameters()
                          << " r=" << std::setw(7) << current_set.radius()
                          << " c=" << current_set.centre())

    if (definitely(current_time >= maximum_time)) {
        result->adjoin_final(current_set);
    } else if (semantics == Semantics::UPPER and this->_configuration->enable_subdivisions() and decide(current_set_radius > this->_configuration->maximum_enclosure_radius())) {
        // Subdivide
        List<EnclosureType> subdivisions = subdivide(current_set);
        for (SizeType i = 0; i != subdivisions.size(); ++i) {
            const EnclosureType& subdivided_set_model = subdivisions[i];
            workload.append({current_time, subdivided_set_model});
        }
    } else if (semantics == Semantics::LOWER and decide(current_set_radius > this->_configuration->maximum_enclosure_radius())) {
        LOGGING_PRINTLN("Terminating lower evolution at time " << current_time << " and set " << current_set << " due to maximum radius being exceeded.")
    } else {
        this->_process_timed_enclosure_step(workload, current_timed_set, maximum_time, semantics, result);
    }
}

Void ContinuousTimeEvolver::
    _process_timed_enclosure_step(WorkloadType::Access& workload,
                                  const TimedEnclosureType& working_timed_set_model,
                                  const TimeType& maximum_time,
                                  Semantics semantics,
                                  SharedPointer<SynchronisedOrbit> result) const {
    LOGGING_SCOPE_CREATE
    typedef EffectiveVectorMultivariateFunction FunctionType;

    EnclosureType current_set = working_timed_set_model.second;
    TimeStepType current_time = working_timed_set_model.first;
    SharedPointer<PreconditionedTaylorSeriesState> carried_preconditioned_state =
        working_timed_set_model.preconditioned_state;
    LOGGING_PRINTLN_AT(1, "working_timed_set_model time = " << current_time)

    LOGGING_PRINTLN("current_time = " << current_time)
    LOGGING_PRINTLN("current_set = " << current_set)

    LOGGING_PRINTLN("box = " << current_set.bounding_box())
    LOGGING_PRINTLN("radius = " << current_set.euclidean_set().bounding_box().radius())

    const IntegratorInterface* integrator = this->_integrator.get();
    const auto* preconditioned_integrator =
        dynamic_cast<const PreconditionedGradedTaylorSeriesIntegrator*>(integrator);

    // The preconditioned integrator carries the local-initial Taylor mapping
    // across steps. Standard reconditioning would deliberately discard part of
    // that dependence, so it is bypassed on this specialised path.
    const SizeType protected_parameters = static_cast<SizeType>(std::count(
        current_set.variable_kinds().begin(), current_set.variable_kinds().end(), EnclosureVariableKind::PARAMETER));
    const SizeType parameter_limit =
        (current_set.configuration().reconditioning_num_blocks() + 1u) * current_set.state_dimension() +
        protected_parameters;
    const bool too_many_parameters = this->system().inclusion() != nullptr &&
                                     current_set.number_of_parameters() > parameter_limit;
    if (preconditioned_integrator == nullptr && this->_configuration->enable_reconditioning() &&
        (too_many_parameters || possibly(norm(current_set.state_function().errors()) > this->_configuration->maximum_spacial_error()))) {
        current_set.recondition();
        workload.append({current_time, current_set});
        return;
    }

    /////////////// Main Evolution ////////////////////////////////

    // Limit the next step to the remaining evolution horizon. Evolution
    // horizons are Real, while integration steps are exact Dyadic values:
    // use validated bounds rather than an unchecked floating-point cast.
    // UPPER reaches at least the horizon; LOWER stays at or before it.
    const auto maximum_time_bounds = maximum_time.get(DoublePrecision());
    const TimeStepType dyadic_end_time =
        (semantics == Semantics::LOWER)
            ? TimeStepType(maximum_time_bounds.lower_raw())
            : TimeStepType(maximum_time_bounds.upper_raw());
    if (current_time >= dyadic_end_time) {
        result->adjoin_final(current_set);
        return;
    }

    // Set evolution parameters
    const StepSizeType configured_step_size = this->_configuration->maximum_step_size();
    const StepSizeType remaining_time = dyadic_end_time - current_time;
    // The integration machinery constructs FloatDP from its suggested
    // Dyadic step and requires exact binary64 representability. An exact
    // difference of accumulated dyadic times need not satisfy that.
    // Round the positive remainder downward before passing it to the
    // integrator, so neither LOWER nor UPPER can overshoot the chosen bound.
    const StepSizeType maximum_step_size = EvolutionDetail::limit_step_to_remaining_time(
        configured_step_size, remaining_time);

    // Get bounding boxes for time and space bounding_box
    auto current_set_bounds = cast_exact_box(current_set.euclidean_set().bounding_box());
    LOGGING_PRINTLN("current_set_bounds = " << current_set_bounds)

    StepSizeType step_size;
    EvolutionDetail::AppliedIntegrationStep applied{{}, current_set, false};
    SharedPointer<PreconditionedTaylorSeriesState> next_preconditioned_state;

    if (this->system().inclusion() != nullptr) {
        const auto common_step = this->_inclusion_step(current_set, current_time, maximum_step_size);
        step_size = common_step.time_step();
        const auto parameters = constant_parameter_coordinates(current_set, this->system().parameter_space().dimension());
        applied = EvolutionDetail::apply_integration_step(
            current_set, common_step, join(current_set.state_function(), parameters), semantics,
            this->_box_integrator != nullptr
                ? EvolutionDetail::IntegrationStepRepresentation::ALLOW_BOX_ENCLOSURES
                : EvolutionDetail::IntegrationStepRepresentation::REQUIRE_STRUCTURED);
    } else {
        const FunctionType dynamic = equation_field(*this->system().equation());
        const SizeType parameter_count = _system->parameter_space().dimension();
        if (_system->has_explicit_time() || parameter_count != 0u) {
            // Bound the equation over absolute time and its fixed-parameter domain.
            // In particular, the bound is a certified enclosure, not an
            // unverified estimate constructed from the initial-state box.
            EulerBounder bounder;
            StepSizeType bounded_step;
            UpperBoxType bounding_box;
            const auto parameters = constant_parameter_coordinates(current_set, parameter_count);
            const auto parameter_domain = cast_exact_box(UpperBoxType(parameters.range()));
            make_lpair(bounded_step, bounding_box) = bounder.compute(
                dynamic, current_set_bounds, current_time, parameter_domain,
                suggest(maximum_step_size));
            const FlowStepIntegratorAdapter step_integrator(*integrator);
            const auto common_step = step_integrator.step(ValidatedIntegrationRequest(
                dynamic, current_set_bounds, current_time, bounded_step,
                parameter_domain, bounding_box));
            step_size = common_step.time_step();
            applied = EvolutionDetail::apply_integration_step(
                current_set, common_step, join(current_set.state_function(), parameters), semantics);
        } else if (preconditioned_integrator != nullptr) {
            if (preconditioned_integrator->diagnostics() && result->reach_size() == 1u && preconditioned_integrator->preconditioning() == TaylorSeriesPreconditioning::QR) {
                PreconditionedGradedTaylorSeriesIntegrator identity_probe =
                    *preconditioned_integrator;
                identity_probe.set_preconditioning(TaylorSeriesPreconditioning::IDENTITY);
                identity_probe.set_diagnostics(true);

                PreconditionedGradedTaylorSeriesIntegrator qr_probe =
                    *preconditioned_integrator;
                qr_probe.set_preconditioning(TaylorSeriesPreconditioning::QR);
                qr_probe.set_diagnostics(true);

                auto identity_state =
                    identity_probe.precondition(current_set.state_function());
                auto qr_state =
                    qr_probe.precondition(current_set.state_function());

                std::cerr << "[PreconditionedSecondStepState]"
                          << " physical_error=" << current_set.state_function().error()
                          << " physical_range=" << current_set.state_function().range()
                          << " identity_domy=" << identity_state.local_domain()
                          << " identity_A=" << identity_state.linear_map()
                          << " identity_normalized_error=" << identity_state.normalised_mapping().error()
                          << " qr_domy=" << qr_state.local_domain()
                          << " qr_A=" << qr_state.linear_map()
                          << " qr_normalized_error=" << qr_state.normalised_mapping().error()
                          << std::endl;

                auto identity_step =
                    identity_probe.step(dynamic, identity_state, suggest(maximum_step_size));
                auto qr_step =
                    qr_probe.step(dynamic, qr_state, suggest(maximum_step_size));

                std::cerr << "[PreconditionedSecondStepResult]"
                          << " identity_h=" << identity_step.time_step()
                          << " identity_flow_error=" << identity_step.flowpipe_mapping().error()
                          << " qr_h=" << qr_step.time_step()
                          << " qr_flow_error=" << qr_step.flowpipe_mapping().error()
                          << std::endl;
            }

            PreconditionedTaylorSeriesState local_state =
                carried_preconditioned_state
                    ? *carried_preconditioned_state
                    : preconditioned_integrator->precondition(current_set.state_function());

            const SizeType diagnostic_step = result->reach_size() + 1u;
            const bool carried_snapshot_step =
                diagnostic_step == 500u || diagnostic_step == 1000u || diagnostic_step == 1500u || diagnostic_step == 2000u;
            if (preconditioned_integrator->carried_expansion_diagnostics() && carried_snapshot_step) {
                print_carried_expansion_snapshot(
                    std::cerr, "input", diagnostic_step,
                    dynamic_cast<const ValidatedVectorMultivariateTaylorFunctionModelDP&>(
                        local_state.normalised_mapping().reference()),
                    preconditioned_integrator->sweeper());
            }

            PreconditionedTaylorSeriesStep local_step =
                preconditioned_integrator->step(
                    dynamic, local_state, suggest(maximum_step_size));

            if (preconditioned_integrator->carried_expansion_diagnostics() && carried_snapshot_step) {
                print_carried_expansion_snapshot(
                    std::cerr, "output", diagnostic_step,
                    dynamic_cast<const ValidatedVectorMultivariateTaylorFunctionModelDP&>(
                        local_step.final_state().normalised_mapping().reference()),
                    preconditioned_integrator->sweeper());
            }

            step_size = local_step.time_step();

            const ValidatedVectorMultivariateFunctionPatch& final_mapping = local_step.evolved_mapping();
            const ValidatedIntegrationStep common_step(
                step_size, local_step.flowpipe_mapping().range(), final_mapping.range(),
                {}, {local_step.flowpipe_mapping()}, final_mapping, current_time,
                local_step.flowpipe_mapping().argument_size() - 1u, false, ExactBoxType(0u), current_time);
            applied = EvolutionDetail::apply_integration_step(current_set, common_step,
                                                              EvolutionDetail::IntegrationStepCoordinates::PARAMETERS, semantics);
            next_preconditioned_state =
                std::make_shared<PreconditionedTaylorSeriesState>(
                    local_step.final_state());

            // Temporary diagnostic: compare endpoint-first composition with the
            // old endpoint-after-full-flowpipe ordering.
            if (preconditioned_integrator->diagnostics() && result->reach_size() < 20u) {
                const auto& final_state = local_step.final_state();
                ValidatedVectorMultivariateFunctionPatch old_order_endpoint =
                    partial_evaluate(
                        local_step.flowpipe_mapping(),
                        local_step.flowpipe_mapping().argument_size() - 1u,
                        step_size);
                std::cerr << "[PreconditionedStepDiagnostic]"
                          << " step=" << result->reach_size()
                          << " carried_state=" << (carried_preconditioned_state ? 1 : 0)
                          << " t=" << current_time
                          << " h=" << step_size
                          << " state_error=" << current_set.state_function().error()
                          << " normalized_error=" << local_state.normalised_mapping().error()
                          << " normalized_range=" << local_state.normalised_mapping().range()
                          << " A=" << local_state.linear_map()
                          << " flowpipe_error=" << local_step.flowpipe_mapping().error()
                          << " old_order_endpoint_error=" << old_order_endpoint.error()
                          << " final_mapping_error=" << final_mapping.error()
                          << " final_normalized_error=" << final_state.normalised_mapping().error()
                          << " final_normalized_range=" << final_state.normalised_mapping().range()
                          << " final_A=" << final_state.linear_map()
                          << std::endl;
            }
        } else {
            const FlowStepIntegratorAdapter step_integrator(*integrator);
            const auto common_step = step_integrator.step(ValidatedIntegrationRequest(
                dynamic, current_set_bounds, suggest(maximum_step_size), current_time));
            step_size = common_step.time_step();
            applied = EvolutionDetail::apply_integration_step(current_set, common_step,
                                                              EvolutionDetail::IntegrationStepCoordinates::STATE, semantics);
        }
    }

    TimeStepType next_time = current_time + TimeStepType(step_size);

    for (const auto& reach: applied.reach) {
        result->adjoin_reach(EnclosureType(reach, current_set.state_space(), current_set.auxiliary_space()));
    }
    EnclosureType next_set(applied.evolved, current_set.state_space(), current_set.auxiliary_space());
    result->adjoin_intermediate(next_set);
    if (next_preconditioned_state) {
        workload.append({next_time, next_set, std::move(next_preconditioned_state)});
    } else {
        workload.append({next_time, next_set});
    }
}

ContinuousTimeEvolverConfiguration::ContinuousTimeEvolverConfiguration() {
    set_maximum_step_size(1);
    set_maximum_enclosure_radius(100.0);
    set_maximum_spacial_error(1e-2);
    set_enable_reconditioning(true);
    set_enable_subdivisions(false);
}

OutputStream&
    ContinuousTimeEvolverConfiguration::_write(OutputStream& os) const {
    os << "ContinuousTimeEvolverConfiguration("
       << "\n maximum_step_size=" << maximum_step_size()
       << ",\n maximum_enclosure_radius=" << maximum_enclosure_radius()
       << ",\n maximum_spacial_error=" << maximum_spacial_error()
       << ",\n enable_reconditioning=" << enable_reconditioning()
       << ",\n enable_subdivisions=" << enable_subdivisions()
       << "\n)";
    return os;
}

} // namespace Ariadne

