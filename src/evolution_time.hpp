/***************************************************************************
 *            evolution_time.hpp
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

/** @file src/evolution_time.hpp
 *  @brief Time-step limiting for evolution algorithms.
 */

#pragma once

#include "numeric/dyadic.hpp"

namespace Ariadne::EvolutionDetail {

/** @brief Limit a configured step to the remaining positive time horizon,
 *  preserving exact binary64 representability required by the integrator.
 */
Dyadic limit_step_to_remaining_time(Dyadic const& configured,Dyadic const& remaining);

}
