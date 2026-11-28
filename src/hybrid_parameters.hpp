/***************************************************************************
 *            hybrid_parameters.hpp
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

/** @file src/hybrid_parameters.hpp
 *  @brief Preserve physical parameter bindings across hybrid modes.
 */

#ifndef ARIADNE_EVOLUTION_HYBRID_PARAMETERS_HPP
#define ARIADNE_EVOLUTION_HYBRID_PARAMETERS_HPP

#include "evolution/hybrid_continuous_time_evolver.hpp"

namespace Ariadne::EvolutionDetail {
const RealSpace& hybrid_parameter_space(const ContinuousTimeHybridAutomatonInterface&, DiscreteLocation);
const ExactBoxType& hybrid_parameter_domain(const ContinuousTimeHybridAutomatonInterface&, DiscreteLocation);
EnclosureFunctionArguments hybrid_function_arguments(const ContinuousTimeHybridAutomatonInterface&, DiscreteLocation);
Void prepare_hybrid_mode(const ContinuousTimeHybridAutomatonInterface&, HybridEnclosure&);
Void configure_transition_target(const ContinuousTimeHybridAutomatonInterface&, TransitionData&);
Void apply_transition_reset(HybridEnclosure&, const TransitionData&);
} // namespace Ariadne::EvolutionDetail

#endif
