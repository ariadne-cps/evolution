/***************************************************************************
 *            test_discrete_time_evolver.cpp
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

/** @file tests/test_discrete_time_evolver.cpp
 *  @brief Test difference equation evolver behavior and numerical invariants.
 */

#include <fstream>
#include <iostream>
#include <memory>

#include "algebra/algebra.hpp"
#include "algebra/differential.hpp"
#include "algebra/matrix.hpp"
#include "algebra/vector.hpp"
#include "evolution/discrete_time_evolver.hpp"
#include "evolution/enclosure.hpp"
#include "evolution/orbit.hpp"
#include "function/constraint.hpp"
#include "function/formula.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "function/taylor_model.hpp"
#include "geometry/box.hpp"
#include "labelled/figure.hpp"
#include "log_verbosity_scope.hpp"
#include "logging/test_runner.hpp"
#include "numeric/rational.hpp"
#include "set/list_set.hpp"
#include "system/difference_algebraic_equation.hpp"
#include "utility/tuple.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;
using namespace std;

class TestDiscreteTimeEvolver : public Ariadne::Testing::ContextFixture {
  public:
    using Ariadne::Testing::ContextFixture::ContextFixture;

  public:
    Void test_stopping_policies() {
        RealVariable x("x"), other("other");
        DifferenceAlgebraicEquation map(List<PrimedRealAssignment>{next(x) = x + 1});
        DiscreteTimeEvolver<ValidatedTag> evolver(map);
        auto initial = evolver.enclosure(ExactBoxType({
            {-1.0_x, 1.0_x}
        }));
        std::unique_ptr<DiscreteTimeEvolver<ValidatedTag>> clone(evolver.clone());
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

    Void test_lower_witness_and_upper_subdivision() {
        const RealVariable x("x");
        const DifferenceAlgebraicEquation translation(List<PrimedRealAssignment>{next(x) = x + 1});
        DiscreteTimeEvolver<ValidatedTag> evolver(translation);
        const auto initial = evolver.enclosure(ExactBoxType({
            {-1.0_x, 1.0_x}
        }));
        evolver.configuration().set_enable_subdivisions(true);
        evolver.configuration().set_enable_premature_termination(false);
        evolver.configuration().set_maximum_enclosure_radius(0.75);

        // Ariadne LOWER is witness-preserving, not INNER: for x0=0,
        // x_n=n must remain inside every propagated unsplit enclosure.
        // A subdivision would not guarantee a witness in *each* child.
        const auto lower = evolver.orbit(initial, Integer(2), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.reach().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.intermediate().size(), 2u);
        const auto final_bounds = lower.final()[0u].euclidean_set().bounding_box()[0u];
        ARIADNE_TESTING_PROVE(final_bounds.lower_bound().raw().get_d() <= 2.0);
        ARIADNE_TESTING_PROVE(final_bounds.upper_bound().raw().get_d() >= 2.0);

        // UPPER may subdivide; the same configured radius must not cause
        // the LOWER algorithm to generate uncertified child enclosures.
        const auto upper = evolver.orbit(initial, Integer(1), Semantics::UPPER);
        ARIADNE_TESTING_PROVE(upper.final().size() > 1u);

        // Nonlinear image: a point witness survives even when the computed
        // enclosure is wider than the exact reachable set.
        const DifferenceAlgebraicEquation quadratic(List<PrimedRealAssignment>{next(x) = x * x});
        DiscreteTimeEvolver<ValidatedTag> nonlinear(quadratic);
        nonlinear.configuration().set_enable_subdivisions(true);
        nonlinear.configuration().set_enable_premature_termination(false);
        nonlinear.configuration().set_maximum_enclosure_radius(0.75);
        const auto nonlinear_lower = nonlinear.orbit(
            nonlinear.enclosure(ExactBoxType({
                {-1.0_x, 1.0_x}
        })),
            Integer(1), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(nonlinear_lower.final().size(), 1u);
        const auto quadratic_bounds =
            nonlinear_lower.final()[0u].euclidean_set().bounding_box()[0u];
        ARIADNE_TESTING_PROVE(quadratic_bounds.lower_bound().raw().get_d() <= 0.0);
        ARIADNE_TESTING_PROVE(quadratic_bounds.upper_bound().raw().get_d() >= 0.0);
    }

    Void test_lower_constrained_source_matches_continuous_admission() {
        const RealVariable x("x");
        const DifferenceAlgebraicEquation equation(
            List<PrimedRealAssignment>{next(x) = x + 1});
        DiscreteTimeEvolver<ValidatedTag> evolver(equation);
        // Retain the vector function owning the indexed component:
        // indexing a temporary may leave a borrowed function reference dangling.
        const auto identity = EffectiveVectorMultivariateFunction::identity(1u);
        ARIADNE_TESTING_EXPECT_EQUAL(identity.result_size(), 1u);
        const auto coordinate = identity[0u];
        auto initial = evolver.enclosure(ExactBoxType({
            {-1.0_x, 1.0_x}
        }));

        // This feasible constrained set has a witness x0=3/8, but none
        // of the centre or endpoint candidate values from the former
        // discrete-only witness search satisfies both constraints.
        // It must be accepted under the same initial-source contract
        // as the ordinary continuous validated LOWER evolver.
        initial.new_positive_state_constraint(coordinate - Real(1 / 4_q));
        initial.new_negative_state_constraint(coordinate - Real(1 / 2_q));

        const auto zero = evolver.orbit(initial, Integer(0), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.final().size(), 1u);

        const auto evolved = evolver.orbit(initial, Integer(1), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.reach().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.final()[0u].euclidean_set().number_of_constraints(), 2u);
        const auto bounds = evolved.final()[0u].euclidean_set().bounding_box()[0u];
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= 1.375);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 1.375);
    }

    Void test_upper_initial_auxiliary_outputs_and_resume() {
        const RealVariable x("x"), a("a"), other("other");
        const DifferenceAlgebraicEquation equation(
            List<PrimedRealAssignment>{next(x) = x / 2}, List<RealAssignment>{let(a) = 3 * x});
        DiscreteTimeEvolver<ValidatedTag> evolver(equation);
        const ExactBoxType box({
            {0, 1}
        });
        for (Nat entry: {0u, 1u, 2u}) {
            const auto source = entry == 0u ? evolver.enclosure(box)
                                            : LabelledEnclosure(box, RealSpace({x}), EnclosureConfiguration(evolver.function_factory()));
            const auto source_domain = source.parameter_domain();
            const auto zero = entry == 2u
                                ? evolver.orbit(RealExpressionBoundedConstraintSet({0 <= x <= 1}, {}), Integer(0), Semantics::UPPER)
                                : evolver.orbit(source, Integer(0), Semantics::UPPER);
            ARIADNE_TESTING_REQUIRE(zero.final().size() == 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(zero.initial().auxiliary_space(), RealSpace({a}));
            ARIADNE_TESTING_REQUIRE(zero.final()[0u].auxiliary_space() == RealSpace({a}));
            const auto first = evolver.orbit(zero.final()[0u], Integer(1), Semantics::UPPER);
            ARIADNE_TESTING_REQUIRE(first.final().size() == 1u);
            const auto resumed = evolver.orbit(first.final()[0u], Integer(1), Semantics::UPPER);
            ARIADNE_TESTING_REQUIRE(resumed.final().size() == 1u);
            for (const auto& final: {zero.final()[0u], first.final()[0u], resumed.final()[0u]}) {
                ARIADNE_TESTING_EXPECT_EQUAL(final.auxiliary_space(), RealSpace({a}));
                ARIADNE_TESTING_EXPECT_EQUAL(final.state_dimension(), 1u);
                const auto residual = final.auxiliary_function()[0u] - 3 * final.state_function()[0u];
                ARIADNE_TESTING_PROVE(subset(cast_exact_interval(residual.range()), ExactIntervalType(-0x1p-24_x, 0x1p-24_x)));
            }
            ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                             {0, 0.25_x}
            }),
                                         cast_exact_box(resumed.final()[0u].state_function().range())));
            ARIADNE_TESTING_EXPECT_EQUAL(source.parameter_domain(), source_domain);
        }
        auto incompatible = LabelledEnclosure(box, RealSpace({x}), EnclosureConfiguration(evolver.function_factory()));
        incompatible.set_auxiliary(RealSpace({other}), equation.auxiliary_function());
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.orbit(incompatible, Integer(0), Semantics::UPPER));
    }

    Void test_evolve() {
        typedef DoublePrecision PR;
        DoublePrecision pr;

        // Set up the map field
        // The Henon map \f$(x,y)\mapsto(a-x^2+by,x)
        Real a = 1.3_dec;
        Real b = 0.3_dec;
        RealVariable x("x"), y("y");
        DifferenceAlgebraicEquation henon({next(x) = a - x * x + b * y, next(y) = x});
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
        DifferenceAlgebraicEquation::TimeType time(6);

        // Set up the evaluators
        DiscreteTimeEvolver<ValidatedTag> evolver(henon);
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
        fig.write("test_discrete_time_evolver_xy");

        LabelledFigure fig2(Axes2d(0 <= TimeVariable() <= 10, -10 <= y <= 10));
        fig2 << orbit.reach();
        fig2.write("test_discrete_time_evolver_ty");
    }

    Void test_subdivide() {
        Real a = 1.3_dec;
        Real b = 0.3_dec;
        RealVariable x("x"), y("y");
        DifferenceAlgebraicEquation henon({next(x) = a - x * x + b * y, next(y) = x});

        RealExpressionBoundedConstraintSet initial_set({0.4_dec <= x <= 0.7_dec, 0.9_dec <= y <= 1.1_dec});

        DifferenceAlgebraicEquation::TimeType time(6);

        DiscreteTimeEvolver<ValidatedTag> evolver(henon);
        evolver.configuration().set_maximum_enclosure_radius(0.1);
        evolver.configuration().set_enable_subdivisions(true);

        auto orbit = evolver.orbit(initial_set, time, Semantics::UPPER);

        ARIADNE_TESTING_PROVE(orbit.final().size() > 1);

        LabelledFigure fig(Axes2d(-10 <= x <= 10, -10 <= y <= 10));
        fig << line_style(true) << fill_colour(cyan) << orbit.reach();
        fig << fill_colour(yellow) << orbit.final();
        fig.write("test_discrete_time_evolver_xy_subdivisions");
    }
};

int main(int argc, const char* argv[]) {
    ARIADNE_TESTING_SUITE("test_discrete_time_evolver");
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeEvolver, test_stopping_policies);
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeEvolver, test_lower_witness_and_upper_subdivision);
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeEvolver, test_lower_constrained_source_matches_continuous_admission);
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeEvolver, test_upper_initial_auxiliary_outputs_and_resume);
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeEvolver, test_evolve);
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeEvolver, test_subdivide);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
