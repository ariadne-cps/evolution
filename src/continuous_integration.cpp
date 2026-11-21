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
    if (step.has_explicit_time() || step.constant_parameter_domain().dimension() != 0u ||
        step.reach_models().size() != 1u || !step.has_evolve_model() ||
        step.reach_time_offset() != step.initial_time()) {
        throw std::invalid_argument("Continuous enclosure bridge requires one autonomous local-time reach patch and an endpoint map");
    }
    const auto& tube = step.reach_models()[0];
    const auto& endpoint = step.evolve_model();
    const SizeType arguments = coordinates == IntegrationStepCoordinates::STATE
                                 ? reach.state_dimension()
                                 : reach.number_of_parameters();
    const IntervalDomainType time_domain(StepSizeType(0u), step.time_step());
    if (step.reach_time_index() != arguments || tube.argument_size() != arguments + 1u ||
        endpoint.argument_size() != arguments || tube.domain()[arguments] != time_domain ||
        tube.result_size() != reach.state_dimension() || endpoint.result_size() != evolved.state_dimension()) {
        throw std::invalid_argument("Integration maps do not match the requested enclosure coordinates");
    }
    if (coordinates == IntegrationStepCoordinates::PARAMETERS &&
        (endpoint.domain() != evolved.parameter_domain() ||
         tube.domain() != product(reach.parameter_domain(), time_domain))) {
        throw std::invalid_argument("Parameterized integration maps must retain the enclosure domain");
    }
    // Apply both certificates independently. In particular, do not reconstruct
    // an endpoint from the tube: preconditioned solvers certify it separately.
    Enclosure next_reach = reach;
    Enclosure next_evolved = evolved;
    if (coordinates == IntegrationStepCoordinates::STATE) {
        next_reach.apply_full_reach_step(tube);
        next_evolved.apply_parameterised_fixed_evolve_step(
            compose(cast_unrestricted(endpoint), next_evolved.state_function()), step.time_step());
    } else {
        next_reach.apply_parameterised_full_reach_step(tube);
        next_evolved.apply_parameterised_fixed_evolve_step(endpoint, step.time_step());
    }
    reach = std::move(next_reach);
    evolved = std::move(next_evolved);
}
} // namespace Ariadne::EvolutionDetail
