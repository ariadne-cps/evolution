/***************************************************************************
 *            test_iterated_map_evolver.cpp
 *
 *  Copyright  2006-20  Pieter Collins
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

/** @file tests/test_iterated_map_evolver.cpp
 *  @brief Test iterated map evolver behavior and numerical invariants.
 */

#include <fstream>
#include "concise_logging.hpp"
#include <iostream>
#include <memory>

#include "algebra/algebra.hpp"
#include "algebra/differential.hpp"
#include "algebra/matrix.hpp"
#include "algebra/vector.hpp"
#include "evolution/enclosure.hpp"
#include "evolution/iterated_map_evolver.hpp"
#include "evolution/orbit.hpp"
#include "function/constraint.hpp"
#include "function/formula.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "function/taylor_model.hpp"
#include "geometry/box.hpp"
#include "labelled/figure.hpp"
#include "log_verbosity_scope.hpp"
#include "numeric/rational.hpp"
#include "set/list_set.hpp"
#include "system/iterated_map.hpp"
#include "utility/tuple.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;
using namespace std;

class TestIteratedMapEvolver : public Ariadne::Testing::ContextFixture {
  public:
    using Ariadne::Testing::ContextFixture::ContextFixture;

  public:
    Void test_stopping_policies() {
        RealVariable x("x"), other("other");
        IteratedMap map(List<PrimedRealAssignment>{next(x) = x + 1});
        IteratedMapEvolver evolver(map);
        auto initial = evolver.enclosure(ExactBoxType({
            {-1.0_x, 1.0_x}
        }));
        std::unique_ptr<IteratedMapEvolver> clone(evolver.clone());
        const auto& const_evolver = *clone;
        ARIADNE_TESTING_EXPECT_EQUAL(const_evolver.configuration().maximum_enclosure_radius(), 100.0);
        ARIADNE_TESTING_EXPECT_EQUAL(clone->system().state_space(), RealSpace({x}));
        auto zero = evolver.orbit(initial, Integer(0), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.reach().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.intermediate().size(), 0u);
        evolver.configuration().set_maximum_enclosure_radius(0.5);
        evolver.configuration().set_enable_premature_termination(true);
        auto stopped = evolver.orbit(initial, Integer(1), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(stopped.final().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(stopped.reach().size(), 0u);
        evolver.configuration().set_maximum_enclosure_radius(2.0);
        auto small = evolver.orbit(initial, Integer(1), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(small.final().size(), 1u);
        evolver.configuration().set_maximum_enclosure_radius(0.5);
        evolver.configuration().set_enable_premature_termination(false);
        {
            LogVerbosityScope diagnostics(100u);
            auto continued = evolver.orbit(initial, Integer(2), Semantics::LOWER);
            ARIADNE_TESTING_EXPECT_EQUAL(continued.final().size(), 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(continued.reach().size(), 2u);
            ARIADNE_TESTING_EXPECT_EQUAL(continued.intermediate().size(), 2u);
            ARIADNE_TESTING_PROVE(definitely(continued.final()[0].euclidean_set().subset(ExactBoxType({
                {ExactDouble(0.99), ExactDouble(3.01)}
            }))));
        }
        auto unsplit = evolver.orbit(initial, Integer(1), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(unsplit.final().size(), 1u);
        auto incompatible = LabelledEnclosure(initial.euclidean_set(), RealSpace({other}));
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.orbit(incompatible, Integer(1), Semantics::UPPER));
    }

    Void test_evolve() {
        typedef DoublePrecision PR;
        DoublePrecision pr;

        // Set up the map field
        // The Henon map \f$(x,y)\mapsto(a-x^2+by,x)
        Real a = 1.3_dec;
        Real b = 0.3_dec;
        RealVariable x("x"), y("y");
        IteratedMap henon({next(x) = a - x * x + b * y, next(y) = x});
        ARIADNE_TESTING_NOTE(henon);

        // Define the initial set
        RealExpressionBoundedConstraintSet initial_set({0.4_dec <= x <= 0.7_dec, 0.9_dec <= y <= 1.1_dec});

        ARIADNE_TESTING_NOTE(initial_set);

        // Function evaluation sanity check
        Vector<FloatApproximation<PR>> p({a, b}, pr);
        Vector<FloatApproximation<PR>> xa({0.5_x, 0.25_x}, pr);
        Vector<FloatApproximation<PR>> hxa = {p[0] - xa[0] * xa[0] + xa[1] * p[1], xa[0]};
        // Approximate equality is three-valued, even when the stored values
        // are identical. Check each computed representation explicitly.
        const auto evaluated = henon.update_function().evaluate(xa);
        ARIADNE_TESTING_REQUIRE(evaluated.size() == hxa.size());
        for (SizeType i = 0; i < hxa.size(); ++i) {
            ARIADNE_TESTING_EXPECT_SAME(evaluated[i], hxa[i]);
        }
        Matrix<FloatApproximation<PR>> dhxa = {
            {-2 * xa[0],      p[1]},
            { {1.0, pr}, {0.0, pr}}
        };
        const auto jacobian = henon.update_function().jacobian(xa);
        ARIADNE_TESTING_REQUIRE(jacobian.row_size() == dhxa.row_size());
        ARIADNE_TESTING_REQUIRE(jacobian.column_size() == dhxa.column_size());
        for (SizeType i = 0; i < dhxa.row_size(); ++i) {
            for (SizeType j = 0; j < dhxa.column_size(); ++j) {
                ARIADNE_TESTING_EXPECT_SAME(jacobian.at(i, j), dhxa.at(i, j));
            }
        }

        // Set up the evolution parameters and grid
        IteratedMap::TimeType time(6);

        // Set up the evaluators
        IteratedMapEvolver evolver(henon);
        evolver.configuration().set_maximum_enclosure_radius(0.25);

        ARIADNE_TESTING_NOTE(evolver.configuration());

        // Compute the orbit
        auto orbit = evolver.orbit(initial_set, time, Semantics::UPPER);

        auto initial_enclosure = LabelledEnclosure(cast_exact_box(initial_set.euclidean_set(henon.state_space()).bounding_box()),
                                                   henon.state_space(), EnclosureConfiguration(evolver.function_factory()));

        // Print the initial, evolve and reach sets
        LabelledFigure fig(Axes2d(-10 <= x <= 10, -10 <= y <= 10));
        fig << line_style(true) << fill_colour(cyan) << orbit.reach();
        fig << fill_colour(yellow) << orbit.final();
        fig << fill_colour(blue) << initial_enclosure;
        fig.write("test_iterated_map_evolver_xy");

        LabelledFigure fig2(Axes2d(0 <= TimeVariable() <= 10, -10 <= y <= 10));
        fig2 << orbit.reach();
        fig2.write("test_iterated_map_evolver_ty");
    }

    Void test_subdivide() {
        Real a = 1.3_dec;
        Real b = 0.3_dec;
        RealVariable x("x"), y("y");
        IteratedMap henon({next(x) = a - x * x + b * y, next(y) = x});

        RealExpressionBoundedConstraintSet initial_set({0.4_dec <= x <= 0.7_dec, 0.9_dec <= y <= 1.1_dec});

        IteratedMap::TimeType time(6);

        IteratedMapEvolver evolver(henon);
        evolver.configuration().set_maximum_enclosure_radius(0.1);
        evolver.configuration().set_enable_subdivisions(true);

        auto orbit = evolver.orbit(initial_set, time, Semantics::UPPER);

        ARIADNE_TESTING_PROVE(orbit.final().size() > 1);

        LabelledFigure fig(Axes2d(-10 <= x <= 10, -10 <= y <= 10));
        fig << line_style(true) << fill_colour(cyan) << orbit.reach();
        fig << fill_colour(yellow) << orbit.final();
        fig.write("test_iterated_map_evolver_xy_subdivisions");
    }
};

int main(int argc, const char* argv[]) {
    ARIADNE_TESTING_SUITE("test_iterated_map_evolver");
    ARIADNE_TESTING_ADD_TEST(TestIteratedMapEvolver, test_stopping_policies);
    ARIADNE_TESTING_ADD_TEST(TestIteratedMapEvolver, test_evolve);
    ARIADNE_TESTING_ADD_TEST(TestIteratedMapEvolver, test_subdivide);
    return run_evolution_test_suite(ariadne_test_suite, argc, argv);
}
