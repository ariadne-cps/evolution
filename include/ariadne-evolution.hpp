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

#include "evolution/configuration_interface.hpp"
#include "evolution/differential_inclusion_evolver.hpp"
#include "evolution/enclosure.hpp"
#include "evolution/enclosure_augmented_paving.hpp"
#include "evolution/evolver_interface.hpp"
#include "evolution/first_order_pde.hpp"
#include "evolution/flowtube.hpp"
#include "evolution/hybrid_enclosure.hpp"
#include "evolution/hybrid_evolver.hpp"
#include "evolution/hybrid_evolver_interface.hpp"
#include "evolution/hybrid_orbit.hpp"
#include "evolution/hybrid_simulator.hpp"
#include "evolution/hybrid_time.hpp"
#include "evolution/iterated_map_evolver.hpp"
#include "evolution/orbit.hpp"
#include "evolution/scalings.hpp"
#include "evolution/traits.hpp"
#include "evolution/vector_field_evolver.hpp"
#include "evolution/vector_field_simulator.hpp"

#include "ariadne-equation-solving.hpp"
#include "ariadne-integration.hpp"
#include "ariadne-system.hpp"

#endif // ARIADNE_ARIADNE_EVOLUTION_HPP
