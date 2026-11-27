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

#include "drawing/box_draw.hpp"
#include "evolution/located_projection.hpp"
#include "logging/test_runner.hpp"
#include "rendering/null_rendering.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

// Record the actual geometric output rather than merely counting callbacks.
class ProjectedBoxCanvas final : public NullCanvas {
  public:
    Nat boundaries = 0u;
    List<Point2d> last_boundary;

    Void fill_boundary(const List<Point2d>& points) override {
        ++boundaries;
        last_boundary = points;
    }
};

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

    Void test_draw_located_box() {
        const DiscreteLocation location(1);
        const RealVariable x("x"), y("y");
        const RealSpace variable_order({y, x});
        const Located<ExactBoxType> located(
            location, ExactBoxType({
                          {-2, -1},
                          { 3,  5}
        }));

        ProjectedBoxCanvas canvas;
        with_projected_location(located, variable_order, x, y,
                                [&](const Projection2d& projection, const ExactBoxType& continuous_box) {
                                    Ariadne::draw(canvas, projection, continuous_box);
                                });
        ARIADNE_TESTING_EXPECT_EQUAL(canvas.boundaries, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(canvas.last_boundary.size(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(canvas.last_boundary[0u], Point2d(3.0, -2.0));
        ARIADNE_TESTING_EXPECT_EQUAL(canvas.last_boundary[2u], Point2d(5.0, -1.0));
        ARIADNE_TESTING_PROVE(located.location() == location);
    }

    Void test_draw_selected_hybrid_box() {
        const DiscreteLocation first(1), second(2);
        const RealVariable x("x"), y("y");
        MonolithicHybridSpace spaces;
        spaces.new_location(first, RealSpace({y, x}));
        spaces.new_location(second, RealSpace({x, y}));

        Hybrid<ExactBoxType> family(
            HybridSpace(spaces),
            [](const DiscreteLocation&, const RealSpace&) {
                return ExactBoxType({
                    {0, 0},
                    {0, 0}
                });
            });
        family[first] = ExactBoxType({
            {-2, -1},
            { 3,  5}
        });
        family[second] = ExactBoxType({
            {-9, -8},
            {-7, -6}
        });
        const Hybrid<ExactBoxType>& selected_family = family;

        ProjectedBoxCanvas first_canvas;
        with_projected_location(selected_family, first, x, y,
                                [&](const Projection2d& projection, const ExactBoxType& continuous_box) {
                                    Ariadne::draw(first_canvas, projection, continuous_box);
                                });
        ARIADNE_TESTING_EXPECT_EQUAL(first_canvas.boundaries, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(first_canvas.last_boundary[0u], Point2d(3.0, -2.0));
        ARIADNE_TESTING_EXPECT_EQUAL(first_canvas.last_boundary[2u], Point2d(5.0, -1.0));

        ProjectedBoxCanvas second_canvas;
        with_projected_location(selected_family, second, x, y,
                                [&](const Projection2d& projection, const ExactBoxType& continuous_box) {
                                    Ariadne::draw(second_canvas, projection, continuous_box);
                                });
        ARIADNE_TESTING_EXPECT_EQUAL(second_canvas.boundaries, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(second_canvas.last_boundary[0u], Point2d(-9.0, -7.0));
        ARIADNE_TESTING_EXPECT_EQUAL(second_canvas.last_boundary[2u], Point2d(-8.0, -6.0));
        ARIADNE_TESTING_EXPECT_EQUAL(selected_family.size(), 2u);
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
    ARIADNE_TESTING_ADD_TEST(TestLocatedProjection, test_draw_located_box);
    ARIADNE_TESTING_ADD_TEST(TestLocatedProjection, test_draw_selected_hybrid_box);
    ARIADNE_TESTING_ADD_TEST(TestLocatedProjection, test_hybrid_location_projection_without_materialisation);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
