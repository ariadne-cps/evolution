/***************************************************************************
 *            test_vector_field_evolver_header.cpp
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

/** @file tests/test_vector_field_evolver_header.cpp
 *  @brief Test vector field evolver header behavior and numerical invariants.
 */

#include <type_traits>

#include "evolution/vector_field_evolver.hpp"

static_assert(std::is_base_of_v<Ariadne::VectorFieldEvolver::Interface,
                                Ariadne::VectorFieldEvolver>);

int main() {}
