/***************************************************************************
 *            continuous_integration.cpp
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

/** @file src/continuous_integration.cpp
 *  @brief Shared validated integration step preparation and application.
 */

#include "continuous_integration.hpp"

#include <optional>
#include <stdexcept>

#include "evolution_time.hpp"
#include "function/function.hpp"
#include "integration/integrator.hpp"
#include "integration/selected_inclusion_integrator.hpp"
#include "integration/validated_step_integrator.hpp"
#include "system/odai.hpp"

namespace Ariadne::EvolutionDetail {
EffectiveVectorMultivariateFunction prepare_equation_field(const EffectiveVectorMultivariateFunction& field,
                                                           SizeType state_dimension, bool timed, SizeType parameter_count) {
    ARIADNE_PRECONDITION(field.argument_size() == state_dimension + (timed ? 1u : 0u) + parameter_count);
    if (timed || parameter_count == 0u)
        return field;
    const auto coordinates = EffectiveScalarMultivariateFunction::coordinates(state_dimension + 1u + parameter_count);
    List<EffectiveScalarMultivariateFunction> arguments;
    for (SizeType i = 0u; i != state_dimension; ++i)
        arguments.append(coordinates[i]);
    for (SizeType i = 0u; i != parameter_count; ++i)
        arguments.append(coordinates[state_dimension + 1u + i]);
    return compose(field, EffectiveVectorMultivariateFunction(arguments));
}

namespace {

ValidatedIntegrationStep at_autonomous_time(const ValidatedIntegrationStep& step, const StepSizeType& time) {
    // Autonomous fields are time-translation invariant. Keep the numerical
    // patch arguments unchanged and record their offset from physical time.
    return ValidatedIntegrationStep(step.time_step(), step.reach(), step.evolve(), {}, step.reach_models(),
                                    step.has_evolve_model() ? std::optional(step.evolve_model()) : std::nullopt,
                                    time, step.reach_time_index(), false, step.constant_parameter_domain(),
                                    step.reach_time_offset() + time - step.initial_time());
}

Pair<StepSizeType, UpperBoxType> interval_field_bound(
    const ValidatedIntervalTaylorBoxMultifunctionDP& field, const ExactBoxType& initial, StepSizeType h,
    bool structured) {
    const auto domain = field.domain();
    if (!definitely(subset(initial, domain))) {
        throw std::invalid_argument("Inclusion state lies outside the interval model's domain");
    }
    const UpperBoxType bound(domain);
    const auto values = [&]() -> UpperBoxType {
        if (structured) {
            // Polynomial expansion can lose coefficient correlations. Bound
            // the actual numerical input field used by the selected solver,
            // which can be wider than the original pointwise box extension.
            const auto data = field.inclusion_integrator_field();
            return UpperBoxType(apply(data.field, product(bound, UpperBoxType(data.inputs))));
        }
        return field.evaluate(bound);
    }();
    for (SizeType attempt = 0u; attempt != 64u; ++attempt) {
        UpperBoxType tube(initial.dimension());
        for (SizeType i = 0u; i != initial.dimension(); ++i) {
            tube[i] = UpperIntervalType(initial[i]) + UpperIntervalType(0u, h) * values[i];
        }
        if (definitely(inside(tube, domain))) {
            return {h, bound};
        }
        h = EvolutionDetail::limit_step_to_remaining_time(h, h / 2u);
    }
    throw BoundingNotFoundException("Cannot certify an inclusion step inside the interval model's domain");
}
ValidatedIntegrationStep inclusion_step_impl(
    const OrdinaryDifferentialAlgebraicInclusion& inclusion, const Enclosure& current,
    const ValidatedVectorMultivariateFunctionPatch& parameters,
    const StepSizeType& time, const StepSizeType& maximum_step,
    const SelectedParametricMultifunctionIntegrator* selected,
    const FirstOrderMultifunctionIntegrator* boxes) {
    const auto initial = hull(cast_exact_box(current.bounding_box()),
                              cast_exact_box(UpperBoxType(current.state_function().range())));
    if (parameters.domain() != current.parameter_domain() ||
        parameters.result_size() != inclusion.parameter_space().dimension()) {
        throw std::invalid_argument("Inclusion parameters must retain their original enclosure coordinates");
    }
    const auto run = [&](const ValidatedIntegrationRequest& request) {
        if (selected) {
            return StructuredInclusionStepIntegratorAdapter(*selected).step(request);
        }
        return MultiflowStepIntegratorAdapter(*boxes).step(request);
    };
    if (!inclusion.has_input_function()) {
        const ValidatedIntervalTaylorBoxMultifunctionDP field(inclusion.interval_taylor_multifunction());
        const auto bounded = interval_field_bound(field, initial, maximum_step, selected != nullptr);
        return at_autonomous_time(run(ValidatedIntegrationRequest(field, initial, bounded.first, bounded.second)), time);
    }
    const SizeType np = inclusion.parameter_space().dimension();
    const auto parameter_domain = cast_exact_box(UpperBoxType(parameters.range()));
    const auto uncertain_arguments = product(parameter_domain, inclusion.inputs());
    EulerBounder bounder;
    auto bounded = inclusion.has_explicit_time()
                     ? bounder.compute(inclusion.function(), initial, time, uncertain_arguments, suggest(maximum_step))
                     : bounder.compute(inclusion.function(), initial, uncertain_arguments, suggest(maximum_step));
    // Refinement can make Euler's bound coincide with the tube. The chosen
    // inclusion solver still verifies strict containment on this padded box.
    for (SizeType i = 0u; i != bounded.second.dimension(); ++i) {
        const FloatDPUpperBound padding = cast_positive(mag(bounded.second[i]) + FloatDPUpperBound(1u, dp)) *
                                          cast_positive(FloatDPUpperBound(0x1p-40_x, dp));
        bounded.second[i] = widen(bounded.second[i], padding);
    }
    // The timed box capability also accepts a zero-dimensional input domain.
    if (inclusion.has_explicit_time() || np != 0u || inclusion.number_of_inputs() == 0u) {
        auto function = inclusion.function();
        if (!inclusion.has_explicit_time()) {
            const SizeType n = inclusion.dimension();
            const auto coordinates = EffectiveScalarMultivariateFunction::coordinates(function.argument_size() + 1u);
            List<EffectiveScalarMultivariateFunction> arguments;
            for (SizeType i = 0u; i != function.argument_size(); ++i) {
                arguments.append(coordinates[i < n ? i : i + 1u]);
            }
            function = compose(function, EffectiveVectorMultivariateFunction(arguments));
        }
        const ValidatedTimedInputFieldBoxAdapterDP field(
            function, cast_exact_box(bounded.second), IntervalDomainType(time, time + bounded.first),
            parameter_domain, inclusion.inputs());
        return run(ValidatedIntegrationRequest(field, initial, time, bounded.first, parameter_domain, bounded.second));
    }
    const ValidatedInputFieldBoxAdapterDP field(inclusion.function(), inclusion.inputs(), cast_exact_box(bounded.second));
    return at_autonomous_time(run(ValidatedIntegrationRequest(field, initial, bounded.first, bounded.second)), time);
}

struct ModelContext {
    Enclosure enclosure;
    std::optional<ValidatedVectorMultivariateFunctionPatch> arguments;
};

ModelContext model_context(const Enclosure& initial, const ExactBoxType& model_domain,
                           IntegrationStepCoordinates coordinates,
                           const ValidatedVectorMultivariateFunctionPatch* argument_map) {
    const bool mapped = argument_map != nullptr || coordinates == IntegrationStepCoordinates::STATE;
    const SizeType base = argument_map != nullptr ? argument_map->result_size()
                                                  : (mapped ? initial.state_dimension() : initial.number_of_parameters());
    if (model_domain.dimension() < base) {
        throw std::invalid_argument("Integration model omits existing enclosure coordinates");
    }
    const ExactBoxType original(project(model_domain, range(0u, base)));
    const ExactBoxType added(project(model_domain, range(base, model_domain.dimension())));
    ModelContext result{initial, std::nullopt};
    if (mapped) {
        const auto map = argument_map != nullptr ? *argument_map : initial.state_function();
        if (!subset(cast_exact_box(UpperBoxType(map.range())), original)) {
            throw std::invalid_argument("Integration argument map exceeds the certified numerical domain");
        }
        result.arguments = embed(map, added);
    } else if (original != initial.parameter_domain()) {
        throw std::invalid_argument("Parameterized integration maps must retain the enclosure domain prefix");
    }
    result.enclosure.adjoin_parameters(added, EnclosureVariableKind::INPUT);
    if (mapped && added.dimension() != 0u) {
        const auto extra = initial.configuration().function_factory().create_projection(
            result.enclosure.parameter_domain(),
            range(initial.number_of_parameters(), result.enclosure.number_of_parameters()));
        result.arguments = join(*result.arguments, extra);
    }
    return result;
}

ValidatedVectorMultivariateFunctionPatch local_reach(
    const Enclosure& initial, const ValidatedVectorMultivariateFunctionPatch& tube,
    SizeType time_index, const StepSizeType& origin) {
    const auto domain = tube.domain();
    const ExactBoxType spatial = product(
        ExactBoxType(project(domain, range(0u, time_index))),
        ExactBoxType(project(domain, range(time_index + 1u, domain.dimension()))));
    const IntervalDomainType time(StepSizeType(domain[time_index].lower_bound()) - origin,
                                  StepSizeType(domain[time_index].upper_bound()) - origin);
    if (time_index + 1u == domain.dimension() && origin == StepSizeType(0u)) {
        return tube;
    }
    const auto local_domain = product(spatial, time);
    const auto& factory = initial.configuration().function_factory();
    const auto clock = factory.create_coordinate(local_domain, spatial.dimension());
    auto arguments = factory.create_zeros(tube.argument_size(), local_domain);
    for (SizeType i = 0u; i != tube.argument_size(); ++i) {
        arguments[i] = i == time_index
                         ? clock + origin
                         : factory.create_coordinate(local_domain, i < time_index ? i : i - 1u);
    }
    return compose(tube, arguments);
}

Enclosure apply_reach_model(const Enclosure& initial,
                            const ValidatedVectorMultivariateFunctionPatch& tube,
                            const ValidatedIntegrationStep& step, IntegrationStepCoordinates coordinates,
                            const ValidatedVectorMultivariateFunctionPatch* argument_map) {
    const SizeType index = step.reach_time_index().value_or(tube.argument_size());
    if (index >= tube.argument_size() || tube.result_size() != initial.state_dimension()) {
        throw std::invalid_argument("Reach model has incompatible state or time coordinates");
    }
    const auto flow = local_reach(initial, tube, index, step.initial_time() - step.reach_time_offset());
    const auto domain = flow.domain();
    const SizeType time_index = domain.dimension() - 1u;
    auto context = model_context(initial, ExactBoxType(project(domain, range(0u, time_index))),
                                 coordinates, argument_map);
    auto mapped_flow = flow;
    if (context.arguments) {
        const auto time_domain = domain[time_index];
        const auto parameters = product(context.enclosure.parameter_domain(), time_domain);
        const auto time = initial.configuration().function_factory().create_coordinate(
            parameters, context.enclosure.number_of_parameters());
        mapped_flow = compose(cast_unrestricted(flow), join(embed(*context.arguments, time_domain), time));
    }
    context.enclosure.apply_parameterised_full_reach_step(mapped_flow);
    return std::move(context.enclosure);
}

Enclosure apply_endpoint_model(const Enclosure& initial, const ValidatedIntegrationStep& step,
                               IntegrationStepCoordinates coordinates,
                               const ValidatedVectorMultivariateFunctionPatch* argument_map) {
    const auto& endpoint = step.evolve_model();
    auto context = model_context(initial, endpoint.domain(), coordinates, argument_map);
    const auto mapped = context.arguments ? compose(cast_unrestricted(endpoint), *context.arguments) : endpoint;
    context.enclosure.apply_parameterised_fixed_evolve_step(mapped, step.time_step());
    return std::move(context.enclosure);
}

Enclosure apply_box(const Enclosure& initial, const UpperBoxType& box,
                    const StepSizeType& duration, bool reach) {
    // A Cartesian outer box retains no state/input correlation. Keep the
    // original feasible parameters and add independent ERROR coordinates;
    // do not present these axes as a numerical trajectory certificate.
    if (!definitely(box.is_bounded())) {
        throw std::invalid_argument("A compact enclosure cannot represent an unbounded box certificate");
    }
    Enclosure result = initial;
    const SizeType first = initial.number_of_parameters();
    result.adjoin_parameters(cast_exact_box(box), EnclosureVariableKind::ERROR);
    const auto state = initial.configuration().function_factory().create_projection(
        result.parameter_domain(), range(first, result.number_of_parameters()));
    if (reach) {
        result.apply_parameterised_full_reach_step(
            embed(state, IntervalDomainType(StepSizeType(0u), duration)));
    } else {
        result.apply_parameterised_fixed_evolve_step(state, duration);
    }
    return result;
}

AppliedIntegrationStep apply_step(const Enclosure& initial, const ValidatedIntegrationStep& step,
                                  IntegrationStepCoordinates coordinates,
                                  const ValidatedVectorMultivariateFunctionPatch* argument_map,
                                  Semantics semantics, IntegrationStepRepresentation representation) {
    if (step.reach().dimension() != initial.state_dimension()) {
        throw std::invalid_argument("Integration step has an incompatible state dimension");
    }
    // Nullable guards must be separate from validated logical expressions,
    // whose overloaded operators may evaluate both operands.
    if (argument_map != nullptr) {
        if (argument_map->domain() != initial.parameter_domain() ||
            argument_map->result_size() != initial.state_dimension() + step.constant_parameter_domain().dimension()) {
            throw std::invalid_argument("Integration argument map must retain the enclosure domain and fixed parameters");
        }
    } else if (coordinates == IntegrationStepCoordinates::STATE && step.constant_parameter_domain().dimension() != 0u) {
        throw std::invalid_argument("A state-coordinate step cannot silently introduce trajectory-constant parameters");
    }
    const bool boxes = !step.has_reach_models() || !step.has_evolve_model();
    if (boxes) {
        if (representation != IntegrationStepRepresentation::ALLOW_BOX_ENCLOSURES || semantics != Semantics::UPPER) {
            throw std::invalid_argument("Box certificates require explicit permission for outer UPPER enclosures");
        }
        if (step.constant_parameter_domain().dimension() != 0u) {
            if (argument_map == nullptr) {
                throw std::invalid_argument("Box certificates with fixed parameters require an explicit argument map");
            }
            const ExactBoxType parameters(project(cast_exact_box(UpperBoxType(argument_map->range())),
                                                  range(initial.state_dimension(), argument_map->result_size())));
            if (!subset(parameters, step.constant_parameter_domain())) {
                throw std::invalid_argument("Box certificate parameters do not cover the enclosure arguments");
            }
        }
    }
    AppliedIntegrationStep result{{}, initial, boxes};
    for (const auto& patch: step.reach_models()) {
        result.reach.append(apply_reach_model(initial, patch, step, coordinates, argument_map));
    }
    if (!step.has_reach_models()) {
        result.reach.append(apply_box(initial, step.reach(), step.time_step(), true));
    }
    result.evolved = step.has_evolve_model()
                       ? apply_endpoint_model(initial, step, coordinates, argument_map)
                       : apply_box(initial, step.evolve(), step.time_step(), false);
    return result;
}
} // namespace

ValidatedIntegrationStep compute_inclusion_step(
    const OrdinaryDifferentialAlgebraicInclusion& inclusion, const Enclosure& current,
    const ValidatedVectorMultivariateFunctionPatch& parameters,
    const StepSizeType& time, const StepSizeType& maximum_step,
    const SelectedParametricMultifunctionIntegrator& solver) {
    return inclusion_step_impl(inclusion, current, parameters, time, maximum_step, &solver, nullptr);
}

ValidatedIntegrationStep compute_inclusion_step(
    const OrdinaryDifferentialAlgebraicInclusion& inclusion, const Enclosure& current,
    const ValidatedVectorMultivariateFunctionPatch& parameters,
    const StepSizeType& time, const StepSizeType& maximum_step,
    const FirstOrderMultifunctionIntegrator& solver) {
    return inclusion_step_impl(inclusion, current, parameters, time, maximum_step, nullptr, &solver);
}

namespace {
// The numerical parameter is a step-local conditioning coordinate: the source
// absolute clock. No synthetic ODAE/ODAI declaration or dot(time)=1 is created.
EffectiveVectorMultivariateFunction clocked_field(const EffectiveVectorMultivariateFunction& field, SizeType n) {
    const auto coordinates = EffectiveScalarMultivariateFunction::coordinates(field.argument_size() + 1u);
    List<EffectiveScalarMultivariateFunction> arguments;
    for (SizeType i = 0u; i != n; ++i)
        arguments.append(coordinates[i]);
    arguments.append(coordinates[n] + coordinates[n + 1u]);
    for (SizeType i = n + 1u; i != field.argument_size(); ++i)
        arguments.append(coordinates[i + 1u]);
    return compose(field, EffectiveVectorMultivariateFunction(arguments));
}

AnchoredIntegrationStep clocked_step(
    const EffectiveVectorMultivariateFunction& field, const ExactBoxType& inputs,
    const Enclosure& current, const ValidatedVectorMultivariateFunctionPatch& model_parameters,
    const StepSizeType& maximum_step, const IntegratorInterface* equation, const SelectedParametricMultifunctionIntegrator* selected,
    const FirstOrderMultifunctionIntegrator* boxes) {
    const auto n = current.state_dimension();
    if (field.result_size() != n || field.argument_size() != n + 1u + model_parameters.result_size() + inputs.dimension() ||
        model_parameters.domain() != current.parameter_domain())
        throw std::invalid_argument("Clocked integration requires [state, global time, inputs]");
    const auto initial = hull(cast_exact_box(current.bounding_box()), cast_exact_box(UpperBoxType(current.state_function().range())));
    const auto origin = product(ExactBoxType({cast_exact_interval(current.time_function().range())}),
                                cast_exact_box(UpperBoxType(model_parameters.range())));
    if (origin.is_empty() || !origin.is_bounded())
        throw std::invalid_argument("Clocked integration requires a finite global time chart");
    const auto dynamic = clocked_field(field, n);
    EulerBounder bounder;
    auto bounded = bounder.compute(dynamic, initial, StepSizeType(0u), product(origin, inputs), suggest(maximum_step));
    if (equation == nullptr) {
        for (SizeType i = 0u; i != bounded.second.dimension(); ++i) {
            const auto padding = cast_positive(mag(bounded.second[i]) + FloatDPUpperBound(1u, dp)) * cast_positive(FloatDPUpperBound(0x1p-40_x, dp));
            bounded.second[i] = widen(bounded.second[i], padding);
        }
    }
    const auto run = [&]() -> ValidatedIntegrationStep {
        if (equation)
            return FlowStepIntegratorAdapter(*equation).step(ValidatedIntegrationRequest(
                dynamic, initial, StepSizeType(0u), bounded.first, origin, bounded.second));
        const ValidatedTimedInputFieldBoxAdapterDP adapter(dynamic, cast_exact_box(bounded.second),
                                                           IntervalDomainType(StepSizeType(0u), bounded.first), origin, inputs);
        const ValidatedIntegrationRequest request(adapter, initial, StepSizeType(0u), bounded.first, origin, bounded.second);
        if (selected)
            return StructuredInclusionStepIntegratorAdapter(*selected).step(request);
        return MultiflowStepIntegratorAdapter(*boxes).step(request);
    };
    return {run(), join(join(current.state_function(), current.time_function()), model_parameters)};
}
} // namespace

AnchoredIntegrationStep compute_clocked_equation_step(const EffectiveVectorMultivariateFunction& field,
                                                      const Enclosure& current, const ValidatedVectorMultivariateFunctionPatch& parameters,
                                                      const StepSizeType& maximum_step, const IntegratorInterface& solver) {
    return clocked_step(field, ExactBoxType(0u), current, parameters, maximum_step, &solver, nullptr, nullptr);
}
AnchoredIntegrationStep compute_clocked_inclusion_step(const OrdinaryDifferentialAlgebraicInclusion& model,
                                                       const Enclosure& current, const ValidatedVectorMultivariateFunctionPatch& parameters,
                                                       const StepSizeType& maximum_step, const SelectedParametricMultifunctionIntegrator& solver) {
    if (!model.has_explicit_time() || model.parameter_space().dimension() != parameters.result_size())
        throw std::invalid_argument("Clocked inclusion step requires explicit time and the original model parameter map");
    return clocked_step(model.function(), model.inputs(), current, parameters, maximum_step, nullptr, &solver, nullptr);
}
AnchoredIntegrationStep compute_clocked_inclusion_step(const OrdinaryDifferentialAlgebraicInclusion& model,
                                                       const Enclosure& current, const ValidatedVectorMultivariateFunctionPatch& parameters,
                                                       const StepSizeType& maximum_step, const FirstOrderMultifunctionIntegrator& solver) {
    if (!model.has_explicit_time() || model.parameter_space().dimension() != parameters.result_size())
        throw std::invalid_argument("Clocked inclusion step requires explicit time and the original model parameter map");
    return clocked_step(model.function(), model.inputs(), current, parameters, maximum_step, nullptr, nullptr, &solver);
}

AppliedIntegrationStep apply_integration_step(const Enclosure& initial, const ValidatedIntegrationStep& step,
                                              IntegrationStepCoordinates coordinates, Semantics semantics,
                                              IntegrationStepRepresentation representation) {
    return apply_step(initial, step, coordinates, nullptr, semantics, representation);
}

AppliedIntegrationStep apply_integration_step(const Enclosure& initial, const ValidatedIntegrationStep& step,
                                              const ValidatedVectorMultivariateFunctionPatch& arguments, Semantics semantics,
                                              IntegrationStepRepresentation representation) {
    return apply_step(initial, step, IntegrationStepCoordinates::PARAMETERS, &arguments, semantics, representation);
}
} // namespace Ariadne::EvolutionDetail
