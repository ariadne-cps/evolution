/***************************************************************************
 *            ariadne-evolution.hpp
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

/** @file ariadne-evolution.hpp
 *  @brief Public umbrella header for evolution algorithms.
 */

#ifndef ARIADNE_ARIADNE_EVOLUTION_HPP
#define ARIADNE_ARIADNE_EVOLUTION_HPP

// Public Evolution workflows, including the installed C++ tutorial,
// require no secondary module or numerical implementation headers.
#include <iostream>

#include "evolution/configuration_interface.hpp"
#include "evolution/continuous_dynamics.hpp"
#include "evolution/continuous_time_approximate_evolver.hpp"
#include "evolution/continuous_time_evolver.hpp"
#include "evolution/discrete_time_approximate_evolver.hpp"
#include "evolution/discrete_time_evolver.hpp"
#include "evolution/discrete_time_inclusion_evolver.hpp"
#include "evolution/enclosure.hpp"
#include "evolution/enclosure_augmented_paving.hpp"
#include "evolution/enclosure_recombiner.hpp"
#include "evolution/enclosure_recombiner_interface.hpp"
#include "evolution/evolver_interface.hpp"
#include "evolution/first_order_pde.hpp"
#include "evolution/flowtube.hpp"
#include "evolution/hybrid_continuous_time_approximate_evolver.hpp"
#include "evolution/hybrid_continuous_time_evolver.hpp"
#include "evolution/hybrid_continuous_time_evolver_interface.hpp"
#include "evolution/hybrid_crossing.hpp"
#include "evolution/hybrid_discrete_time_approximate_evolver.hpp"
#include "evolution/hybrid_discrete_time_evolver.hpp"
#include "evolution/hybrid_discrete_time_predicate_evaluator.hpp"
#include "evolution/hybrid_discrete_time_step_executor.hpp"
#include "evolution/hybrid_discrete_time_step_image.hpp"
#include "evolution/hybrid_discrete_time_step_planner.hpp"
#include "evolution/hybrid_discrete_time_upper_scheduler.hpp"
#include "evolution/hybrid_enclosure.hpp"
#include "evolution/hybrid_evolution_backend.hpp"
#include "evolution/hybrid_evolution_tree.hpp"
#include "evolution/hybrid_orbit.hpp"
#include "evolution/hybrid_time.hpp"
#include "evolution/input_sampling.hpp"
#include "evolution/orbit.hpp"
#include "evolution/scalings.hpp"
#include "evolution/traits.hpp"
#include "function/function.hpp"
#include "hybrid/hybrid_set.hpp"
#include "integration/integrator.hpp"
#include "numeric/float_bounds.hpp"
#include "system/discrete_time_hybrid_automaton.hpp"

#include "ariadne-equation-solving.hpp"
#include "ariadne-integration.hpp"
#include "ariadne-system.hpp"

#endif // ARIADNE_ARIADNE_EVOLUTION_HPP
