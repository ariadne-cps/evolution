/***************************************************************************
 *            continuous_integration.cpp
 *  Copyright  2026  Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#include "continuous_integration.hpp"

#include <stdexcept>

#include "evolution/enclosure.hpp"
#include "function/function.hpp"
#include "function/function_patch.hpp"
#include "integration/validated_integration_step.hpp"

namespace Ariadne::EvolutionDetail {

namespace {
void apply_step(Enclosure& reach, Enclosure& evolved,
                const ValidatedIntegrationStep& step,
                IntegrationStepCoordinates coordinates,
                const ValidatedVectorMultivariateFunctionPatch* argument_map) {
    if (step.reach_models().size() != 1u || !step.has_evolve_model()) {
        throw std::invalid_argument("Continuous enclosure bridge requires one certified reach patch and an endpoint map");
    }
    const auto& tube = step.reach_models()[0];
    const auto& endpoint = step.evolve_model();
    const SizeType arguments = argument_map != nullptr
                                 ? argument_map->result_size()
                                 : (coordinates == IntegrationStepCoordinates::STATE
                                        ? reach.state_dimension()
                                        : reach.number_of_parameters());
    // Keep nullable guards separate from validated logical expressions:
    // their overloaded operators may evaluate both operands.
    if (argument_map != nullptr) {
        if (argument_map->domain() != reach.parameter_domain() ||
            argument_map->domain() != evolved.parameter_domain() ||
            arguments != reach.state_dimension() + step.constant_parameter_domain().dimension()) {
            throw std::invalid_argument("Integration argument map must retain the enclosure domain and fixed parameters");
        }
    }
    // The numerical solver can put time between state and fixed-parameter
    // coordinates (first step) or after all retained parameters (continuation).
    // Do not confuse the time index with the number of enclosure parameters.
    const SizeType time_index = step.reach_time_index().value_or(tube.argument_size());
    const StepSizeType patch_time_origin = step.initial_time() - step.reach_time_offset();
    const IntervalDomainType patch_time_domain(patch_time_origin, patch_time_origin + step.time_step());
    const IntervalDomainType local_time_domain(StepSizeType(0u), step.time_step());
    if (tube.argument_size() != arguments + 1u || time_index >= tube.argument_size() ||
        endpoint.argument_size() != arguments || tube.domain()[time_index] != patch_time_domain ||
        tube.result_size() != reach.state_dimension() || endpoint.result_size() != evolved.state_dimension()) {
        throw std::invalid_argument("Integration maps do not match the requested enclosure coordinates and time interval");
    }
    // Removing only the time argument must recover every original state or
    // uncertainty coordinate, including trajectory-constant parameters.
    const ExactBoxType model_coordinates = product(
        ExactBoxType(project(tube.domain(), range(0u, time_index))),
        ExactBoxType(project(tube.domain(), range(time_index + 1u, tube.argument_size()))));
    if (argument_map == nullptr && coordinates == IntegrationStepCoordinates::PARAMETERS &&
        (endpoint.domain() != evolved.parameter_domain() ||
         model_coordinates != reach.parameter_domain())) {
        throw std::invalid_argument("Parameterized integration maps must retain the enclosure domain");
    }
    if (argument_map != nullptr) {
        if (endpoint.domain() != model_coordinates ||
            !subset(cast_exact_box(UpperBoxType(argument_map->range())), model_coordinates)) {
            throw std::invalid_argument("Integration argument map exceeds the certified numerical domain");
        }
    }

    // Enclosure requires [original coordinates, local time]. Normalize both
    // the argument ordering and the time origin using a validated composition.
    // This preserves every fixed-parameter coordinate, without boxification,
    // and does not alter the independent endpoint certificate.
    auto local_tube = tube;
    if (time_index != arguments || patch_time_origin != StepSizeType(0u)) {
        const auto& factory = reach.configuration().function_factory();
        const auto local_domain = product(model_coordinates, local_time_domain);
        const auto local_time = factory.create_coordinate(local_domain, arguments);
        const auto first_argument = time_index == 0u
                                      ? local_time + patch_time_origin
                                      : factory.create_coordinate(local_domain, 0u);
        auto original_arguments = ValidatedVectorMultivariateFunctionPatch(1u, first_argument);
        for (SizeType i = 1u; i != tube.argument_size(); ++i) {
            const auto argument = i == time_index
                                    ? local_time + patch_time_origin
                                    : factory.create_coordinate(local_domain,
                                                                i < time_index ? i : i - 1u);
            original_arguments = join(original_arguments, argument);
        }
        local_tube = compose(tube, original_arguments);
    }

    // Work on independent copies; a failure leaves both inputs unchanged.
    Enclosure next_reach = reach;
    Enclosure next_evolved = evolved;
    if (argument_map != nullptr) {
        const auto domain = product(argument_map->domain(), local_time_domain);
        const auto time = reach.configuration().function_factory().create_coordinate(
            domain, argument_map->argument_size());
        const auto tube_arguments = join(embed(*argument_map, local_time_domain), time);
        next_reach.apply_parameterised_full_reach_step(
            compose(cast_unrestricted(local_tube), tube_arguments));
        next_evolved.apply_parameterised_fixed_evolve_step(
            compose(cast_unrestricted(endpoint), *argument_map), step.time_step());
    } else if (coordinates == IntegrationStepCoordinates::STATE) {
        next_reach.apply_full_reach_step(local_tube);
        next_evolved.apply_parameterised_fixed_evolve_step(
            compose(cast_unrestricted(endpoint), next_evolved.state_function()), step.time_step());
    } else {
        next_reach.apply_parameterised_full_reach_step(local_tube);
        next_evolved.apply_parameterised_fixed_evolve_step(endpoint, step.time_step());
    }
    reach = std::move(next_reach);
    evolved = std::move(next_evolved);
}
} // namespace

void apply_integration_step(Enclosure& reach, Enclosure& evolved,
                            const ValidatedIntegrationStep& step, IntegrationStepCoordinates coordinates) {
    if (coordinates == IntegrationStepCoordinates::STATE) {
        if (step.constant_parameter_domain().dimension() != 0u) {
            throw std::invalid_argument("A state-coordinate step cannot silently introduce trajectory-constant parameters");
        }
    }
    apply_step(reach, evolved, step, coordinates, nullptr);
}

void apply_integration_step(Enclosure& reach, Enclosure& evolved,
                            const ValidatedIntegrationStep& step,
                            const ValidatedVectorMultivariateFunctionPatch& arguments) {
    apply_step(reach, evolved, step, IntegrationStepCoordinates::PARAMETERS, &arguments);
}
} // namespace Ariadne::EvolutionDetail
