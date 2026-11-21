/***************************************************************************
 *            continuous_integration.cpp
 *  Copyright  2026  Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#include "continuous_integration.hpp"

#include <optional>
#include <stdexcept>

#include "function/function.hpp"
#include "integration/validated_integration_step.hpp"

namespace Ariadne::EvolutionDetail {
namespace {

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
