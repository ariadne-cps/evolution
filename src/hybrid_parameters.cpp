/***************************************************************************
 *            hybrid_parameters.cpp
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

/** @file src/hybrid_parameters.cpp
 *  @brief Preserve physical parameter bindings across hybrid modes.
 */

#include "hybrid_parameters.hpp"

#include <stdexcept>

#include "system/odae.hpp"
#include "system/odai.hpp"

namespace Ariadne::EvolutionDetail {
const RealSpace& hybrid_parameter_space(const ContinuousTimeHybridAutomatonInterface& system, DiscreteLocation location) {
    if (const auto* equation = system.equation(location))
        return equation->parameter_space();
    if (const auto* inclusion = system.inclusion(location))
        return inclusion->parameter_space();
    static const RealSpace empty;
    return empty;
}
const ExactBoxType& hybrid_parameter_domain(const ContinuousTimeHybridAutomatonInterface& system, DiscreteLocation location) {
    if (const auto* equation = system.equation(location))
        return equation->constant_parameter_domain();
    if (const auto* inclusion = system.inclusion(location))
        return inclusion->constant_parameter_domain();
    static const ExactBoxType empty;
    return empty;
}
EnclosureFunctionArguments hybrid_function_arguments(const ContinuousTimeHybridAutomatonInterface& system, DiscreteLocation location) {
    const bool timed = system.time_variable(location).has_value();
    if (hybrid_parameter_space(system, location).dimension() != 0u)
        return timed ? EnclosureFunctionArguments::STATE_TIME_PARAMETERS : EnclosureFunctionArguments::STATE_PARAMETERS;
    return timed ? EnclosureFunctionArguments::STATE_TIME : EnclosureFunctionArguments::STATE;
}
namespace {
Void prepare_parameters(HybridEnclosure& set, const RealSpace& space, const ExactBoxType& domain) {
    if (space.dimension() == 0u)
        return;
    const auto& current = set.continuous_set();
    if (current.has_constant_parameter_binding()) {
        if (current.constant_parameter_space() != space)
            throw std::invalid_argument("Hybrid modes require the same ordered physical parameter tuple");
        if (!definitely(subset(current.constant_parameter_domain(), domain)))
            throw std::invalid_argument("Hybrid target parameter domain does not contain the bound physical domain");
    } else {
        auto bound = current;
        bound.adjoin_constant_parameters(space, domain);
        set.replace_continuous_set(bound);
    }
}
} // namespace

Void prepare_hybrid_mode(const ContinuousTimeHybridAutomatonInterface& system, HybridEnclosure& set) {
    const auto location = set.location();
    const auto& space = hybrid_parameter_space(system, location);
    prepare_parameters(set, space, hybrid_parameter_domain(system, location));
    set.set_auxiliary(system.continuous_auxiliary_space(location).variables(), system.auxiliary_function(location),
                      hybrid_function_arguments(system, location), space.dimension());
}
Void configure_transition_target(const ContinuousTimeHybridAutomatonInterface& system, TransitionData& transition) {
    const auto target = transition.target;
    transition.target_auxiliary_arguments = hybrid_function_arguments(system, target);
    transition.target_parameter_space = hybrid_parameter_space(system, target);
    transition.target_parameter_domain = hybrid_parameter_domain(system, target);
}
Void apply_transition_reset(HybridEnclosure& set, const TransitionData& transition) {
    prepare_parameters(set, transition.target_parameter_space, transition.target_parameter_domain);
    set.apply_reset(transition.event, transition.target, transition.target_space, transition.reset_function,
                    transition.target_auxiliary_space, transition.target_auxiliary_function,
                    transition.source_arguments, transition.target_auxiliary_arguments);
}
} // namespace Ariadne::EvolutionDetail
