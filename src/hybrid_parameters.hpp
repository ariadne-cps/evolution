/***************************************************************************
 *            hybrid_parameters.hpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#ifndef ARIADNE_EVOLUTION_HYBRID_PARAMETERS_HPP
#define ARIADNE_EVOLUTION_HYBRID_PARAMETERS_HPP

#include "evolution/hybrid_continuous_time_evolver.hpp"

namespace Ariadne::EvolutionDetail {
const RealSpace& hybrid_parameter_space(const HybridAutomatonInterface&, DiscreteLocation);
const ExactBoxType& hybrid_parameter_domain(const HybridAutomatonInterface&, DiscreteLocation);
EnclosureFunctionArguments hybrid_function_arguments(const HybridAutomatonInterface&, DiscreteLocation);
Void prepare_hybrid_mode(const HybridAutomatonInterface&, HybridEnclosure&);
Void configure_transition_target(const HybridAutomatonInterface&, TransitionData&);
Void apply_transition_reset(HybridEnclosure&, const TransitionData&);
} // namespace Ariadne::EvolutionDetail

#endif
