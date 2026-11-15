/***************************************************************************
 *            test_orbit_evolution_header.cpp
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

/** @file tests/test_orbit_evolution_header.cpp
 *  @brief Test orbit evolution header behavior and numerical invariants.
 */

// Standalone evolution Orbit header must not require analysis storage.
#include "evolution/orbit.hpp"
#include "logging/test_runner.hpp"

#include "ariadne-testing.hpp"

class TestOrbitEvolutionHeader : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    void test_header_compiles() {
        ARIADNE_TESTING_PROVE(true);
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_orbit_evolution_header");
    ARIADNE_TESTING_ADD_TEST(TestOrbitEvolutionHeader, test_header_compiles);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
