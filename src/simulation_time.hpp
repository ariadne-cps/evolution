/***************************************************************************
 *            simulation_time.hpp
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

/** @file src/simulation_time.hpp
 *  @brief Assemble external-time arguments for approximate evolution.
 */

#ifndef ARIADNE_EVOLUTION_SIMULATION_TIME_HPP
#define ARIADNE_EVOLUTION_SIMULATION_TIME_HPP
#include "algebra/vector.hpp"
namespace Ariadne::EvolutionDetail {
template<class X>
Vector<X> simulation_arguments(const Vector<X>& state, bool timed, const X& time) {
    return timed ? join(state, Vector<X>({time})) : state;
}
} // namespace Ariadne::EvolutionDetail
#endif
