/***************************************************************************
 *            tests/test_located_projection.cpp
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

/** @file tests/test_located_projection.cpp
 *  @brief Test selection and continuous coordinate projection by location.
 */

#include "evolution/located_projection.hpp"
#include "logging/test_runner.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

class TestLocatedProjection : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_located_projection() {
        const DiscreteLocation q(1);
        const RealVariable x("x"), y("y");
        const RealSpace space({y, x});
        const Located<int> value(q, 37);
        Nat visited = 0u;
        with_projected_location(value, space, x, y,
            [&](const Projection2d& projection, const int& selected) {
                ARIADNE_TESTING_EXPECT_EQUAL(projection.n, 2u);
                ARIADNE_TESTING_EXPECT_EQUAL(projection.i, 1u);
                ARIADNE_TESTING_EXPECT_EQUAL(projection.j, 0u);
                ARIADNE_TESTING_EXPECT_EQUAL(selected, 37);
                ++visited;
            });
        ARIADNE_TESTING_EXPECT_EQUAL(visited, 1u);
        ARIADNE_TESTING_PROVE(value.location() == q);
    }

    Void test_hybrid_location_projection_without_materialisation() {
        const DiscreteLocation first(1), second(2);
        const RealVariable x("x"), y("y");
        MonolithicHybridSpace spaces;
        spaces.new_location(first, RealSpace({y, x}));
        spaces.new_location(second, RealSpace({x, y}));
        Nat constructed = 0u;
        Hybrid<int> family(HybridSpace(spaces),
            [&](const DiscreteLocation&, const RealSpace& space) {
                ++constructed;
                return static_cast<int>(space.dimension());
            });
        ARIADNE_TESTING_PROVE(not family.has_value(first));
        ARIADNE_TESTING_PROVE(not family.has_value(second));
        family[first] = 41;
        ARIADNE_TESTING_EXPECT_EQUAL(constructed, 1u);
        const Hybrid<int>& selected = family;

        Nat visited = 0u;
        with_projected_location(selected, first, x, y,
            [&](const Projection2d& projection, const int& value) {
                ARIADNE_TESTING_EXPECT_EQUAL(projection.n, 2u);
                ARIADNE_TESTING_EXPECT_EQUAL(projection.i, 1u);
                ARIADNE_TESTING_EXPECT_EQUAL(projection.j, 0u);
                ARIADNE_TESTING_EXPECT_EQUAL(value, 41);
                ++visited;
            });
        with_projected_location(selected, first, RealSpace({x, y}), x, y,
            [&](const Projection2d& projection, const int& value) {
                ARIADNE_TESTING_EXPECT_EQUAL(projection.i, 0u);
                ARIADNE_TESTING_EXPECT_EQUAL(projection.j, 1u);
                ARIADNE_TESTING_EXPECT_EQUAL(value, 41);
                ++visited;
            });
        ARIADNE_TESTING_EXPECT_EQUAL(visited, 2u);
        ARIADNE_TESTING_PROVE(not selected.has_value(second));
        ARIADNE_TESTING_EXPECT_EQUAL(constructed, 1u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(
            with_projected_location(selected, second, x, y,
                [&](const Projection2d&, const int&) {}));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_located_projection");
    ARIADNE_TESTING_ADD_TEST(TestLocatedProjection, test_located_projection);
    ARIADNE_TESTING_ADD_TEST(TestLocatedProjection, test_hybrid_location_projection_without_materialisation);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
