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

void apply_integration_step(Enclosure& reach, Enclosure& evolved,
                            const ValidatedIntegrationStep& step, IntegrationStepCoordinates coordinates) {
    if (step.reach_models().size() != 1u || !step.has_evolve_model()) {
        throw std::invalid_argument("Continuous enclosure bridge requires one certified reach patch and an endpoint map");
    }
    const auto& tube = step.reach_models()[0];
    const auto& endpoint = step.evolve_model();
    const SizeType arguments = coordinates == IntegrationStepCoordinates::STATE
                                 ? reach.state_dimension()
                                 : reach.number_of_parameters();
    const StepSizeType patch_time_origin = step.initial_time() - step.reach_time_offset();
    const IntervalDomainType patch_time_domain(patch_time_origin, patch_time_origin + step.time_step());
    const IntervalDomainType local_time_domain(StepSizeType(0u), step.time_step());
    if (step.reach_time_index() != arguments || tube.argument_size() != arguments + 1u ||
        endpoint.argument_size() != arguments || tube.domain()[arguments] != patch_time_domain ||
        tube.result_size() != reach.state_dimension() || endpoint.result_size() != evolved.state_dimension()) {
        throw std::invalid_argument("Integration maps do not match the requested enclosure coordinates and time interval");
    }
    if (coordinates == IntegrationStepCoordinates::STATE &&
        step.constant_parameter_domain().dimension() != 0u) {
        throw std::invalid_argument("A state-coordinate step cannot silently introduce trajectory-constant parameters");
    }
    if (coordinates == IntegrationStepCoordinates::PARAMETERS &&
        (endpoint.domain() != evolved.parameter_domain() ||
         tube.domain() != product(reach.parameter_domain(), patch_time_domain))) {
        throw std::invalid_argument("Parameterized integration maps must retain the enclosure domain");
    }

    // Enclosure uses a local [0,h] temporal parameter. A validated numerical
    // reach patch may instead use absolute time (e.g. [t0,t0+h]). Rebase the
    // time argument without modifying the state or fixed-parameter coordinates.
    // The separately certified endpoint map remains unchanged.
    auto local_tube = tube;
    if (patch_time_origin != StepSizeType(0u)) {
        const auto& factory = reach.configuration().function_factory();
        const auto local_domain = product(
            ExactBoxType(project(tube.domain(), range(0u, arguments))), local_time_domain);
        const auto coordinates_map = factory.create_projection(local_domain, range(0u, arguments));
        const auto local_time = factory.create_coordinate(local_domain, arguments);
        local_tube = compose(tube, join(coordinates_map, local_time + patch_time_origin));
    }

    // Work on independent copies; a failure leaves both inputs unchanged.
    Enclosure next_reach = reach;
    Enclosure next_evolved = evolved;
    if (coordinates == IntegrationStepCoordinates::STATE) {
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
} // namespace Ariadne::EvolutionDetail
