/***************************************************************************
 *            test_continuous_time_approximate_evolver.cpp
 *
 *  Copyright  2006-20  Luca Geretti, Mirko Albanese
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

/** @file tests/test_continuous_time_approximate_evolver.cpp
 *  @brief Test vector field evolver behavior and numerical invariants.
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include "algebra/algebra.hpp"
#include "algebra/matrix.hpp"
#include "algebra/vector.hpp"
#include "evolution/continuous_time_approximate_evolver.hpp"
#include "evolution/enclosure.hpp"
#include "evolution/orbit.hpp"
#include "function/constraint.hpp"
#include "function/function.hpp"
#include "geometry/box.hpp"
#include "integration/integrator.hpp"
#include "labelled/figure.hpp"
#include "log_verbosity_scope.hpp"
#include "logging/test_runner.hpp"
#include "numeric/rational.hpp"
#include "set/list_set.hpp"
#include "utility/tuple.hpp"

#include "ariadne-testing.hpp"
#include "expression-set/expression_set.hpp"

using namespace Ariadne;
using namespace std;

class TestContinuousTimeApproximateEvolver : public Ariadne::Testing::ContextFixture {
  public:
    using Ariadne::Testing::ContextFixture::ContextFixture;

  public:
    Void test_linear_flow_and_input_labels() {
        RealVariable x("x"), y("y"), z("z"), missing("missing");
        OrdinaryDifferentialAlgebraicEquation field({dot(x) = 1, dot(y) = 0}, {let(z) = x + 2 * y});
        ContinuousTimeEvolver<ApproximateTag> evolver(field);
        evolver.configuration().set_step_size(0.25);
        using PointType = ContinuousTimeEvolver<ApproximateTag>::ApproximatePointType;
        using Points = ContinuousTimeEvolver<ApproximateTag>::ApproximateListPointType;
        PointType start(RealSpace({x, y}), FloatDPApproximationPoint({0.0_x, 2.0_x}, dp));
        PointType reordered(RealSpace({y, x}), FloatDPApproximationPoint({2.0_x, 0.0_x}, dp));
        for (unsigned verbosity: {0u, 100u}) {
            LogVerbosityScope logging(verbosity);
            auto orbit = evolver.orbit(Points{start, reordered}, Real(1));
            ARIADNE_TESTING_EXPECT_EQUAL(orbit.curves().size(), 2u);
            for (const auto& curve: orbit.curves()) {
                const auto& endpoint = curve.end();
                ARIADNE_TESTING_EXPECT_EQUAL(endpoint->first.get_d(), 1.0);
                ARIADNE_TESTING_PROVE(std::abs(endpoint->second[curve.space()[x]].raw().get_d() - 1.0) < 1e-12);
                ARIADNE_TESTING_PROVE(std::abs(endpoint->second[curve.space()[y]].raw().get_d() - 2.0) < 1e-12);
                ARIADNE_TESTING_PROVE(std::abs(endpoint->second[curve.space()[z]].raw().get_d() - 5.0) < 1e-12);
            }
        }
        ARIADNE_TESTING_EXPECT_EQUAL(evolver.orbit(Points{}, Real(1)).curves().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(evolver.orbit(Points{start}, Real(0)).curve(0).size(), 1u);
        PointType incompatible(RealSpace({missing, y}), FloatDPApproximationPoint({0.0_x, 2.0_x}, dp));
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.orbit(Points{incompatible}, Real(1)));
    }

    Void test_configuration_and_sampling() {
        RealVariable x("x"), y("y"), unknown("unknown");
        OrdinaryDifferentialAlgebraicEquation field({dot(x) = 0, dot(y) = 0});
        ContinuousTimeEvolver<ApproximateTag> evolver(field);
        auto& config = evolver.configuration();
        ARIADNE_TESTING_EXPECT_EXCEPTION(config.set_step_size(0.0));
        ARIADNE_TESTING_EXPECT_EXCEPTION(config.set_step_size(-1.0));
        ARIADNE_TESTING_EXPECT_EXCEPTION(config.set_sampling_count(x, 0u));
        config.set_sampling_count(x, 2u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(config.set_sampling_count(x, 3u));
        ARIADNE_TESTING_EXPECT_EXCEPTION(config.set_sampling_counts(SamplingCounts{
            {x, 0u}
        }));
        ARIADNE_TESTING_EXPECT_EQUAL(config.sampling_counts().size(), 1u);
        config.clear_sampling_counts();
        ARIADNE_TESTING_PROVE(config.sampling_counts().empty());
        config.set_sampling_counts(SamplingCounts{});
        ARIADNE_TESTING_PROVE(config.sampling_counts().empty());
        config.set_sampling_counts(SamplingCounts{
            {x, 2u}
        });
        const auto& configured = evolver;
        ARIADNE_TESTING_EXPECT_EQUAL(configured.configuration().sampling_counts().size(), 1u);
        UpperBoxType box({
            {0, 1},
            {0, 1}
        });
        auto divided = evolver.orbit(box, Real(0));
        ARIADNE_TESTING_EXPECT_EQUAL(divided.curves().size(), 2u);
        config.set_sampling_counts(SamplingCounts{
            {unknown, 2u}
        });
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.orbit(box, Real(0)));
        config.clear_sampling_counts();
        auto sampled = evolver.orbit(box, Real(0));
        ARIADNE_TESTING_EXPECT_EQUAL(sampled.curves().size(), 1u);
        std::ostringstream description;
        const ConfigurationInterface& base = config;
        description << base;
        ARIADNE_TESTING_PROVE(description.str().find("sampling_counts") != std::string::npos);
    }

    Void test_multiple_trajectories_x() {
        typedef ContinuousTimeEvolver<ApproximateTag>::ApproximateListPointType ListPointType;

        Real mu = Dyadic(0.5_x);
        RealVariable x("x"), y("y"), z("z");

        OrdinaryDifferentialAlgebraicEquation vanderpol({dot(x) = y, dot(y) = mu * (1 - x * x) * y - x}, {let(z) = sqrt(sqr(x) + sqr(y))});
        ARIADNE_TESTING_NOTE(vanderpol);

        RealExpressionBoundedConstraintSet initial_set({-2.0_dec <= x <= -1.5_dec, 0.5_dec <= y <= 1});

        Real time = 1.0_dec;

        ContinuousTimeEvolver<ApproximateTag> evolver(vanderpol);
        evolver.configuration().set_step_size(0.05);

        Nat trajectory_number = 4;

        for (Nat t = 1; t <= trajectory_number; t++) {
            evolver.configuration().set_sampling_count(x, t);
            ARIADNE_TESTING_NOTE(evolver.configuration());

            Orbit<ListPointType> orbit = evolver.orbit(initial_set, time);

            List<LabelledInterpolatedCurve> curves = orbit.curves();

            ARIADNE_TESTING_EXPECT_EQUAL(curves.size(), t);

            evolver.configuration().clear_sampling_counts();
        }
    }

    Void test_multiple_trajectories_xy() {
        typedef ContinuousTimeEvolver<ApproximateTag>::ApproximateListPointType ListPointType;

        Real mu = Dyadic(0.5_x);
        RealVariable x("x"), y("y"), z("z");

        OrdinaryDifferentialAlgebraicEquation vanderpol({dot(x) = y, dot(y) = mu * (1 - x * x) * y - x}, {let(z) = sqrt(sqr(x) + sqr(y))});
        ARIADNE_TESTING_NOTE(vanderpol);

        RealExpressionBoundedConstraintSet initial_set({-2.0_dec <= x <= -1.5_dec, 0.5_dec <= y <= 1});

        Real time = 1.0_dec;

        ContinuousTimeEvolver<ApproximateTag> evolver(vanderpol);
        evolver.configuration().set_step_size(0.05);

        Nat trajectory_number = 4;

        for (Nat first = 1; first <= trajectory_number; first++) {
            for (Nat second = 1; second <= trajectory_number; second++) {
                evolver.configuration().set_sampling_count(x, first);
                evolver.configuration().set_sampling_count(y, second);
                ARIADNE_TESTING_NOTE(evolver.configuration());

                Orbit<ListPointType> orbit = evolver.orbit(initial_set, time);

                List<LabelledInterpolatedCurve> curves = orbit.curves();

                ARIADNE_TESTING_EXPECT_EQUAL(curves.size(), first * second);

                evolver.configuration().clear_sampling_counts();
            }
        }
    }

    Void test_multiple_real_expr_vdp_cycle() {
        typedef ContinuousTimeEvolver<ApproximateTag>::ApproximateListPointType ListPointType;

        Real mu = Dyadic(0.5_x);
        RealVariable x("x"), y("y"), z("z");

        OrdinaryDifferentialAlgebraicEquation vanderpol({dot(x) = y, dot(y) = mu * (1 - x * x) * y - x}, {let(z) = sqrt(sqr(x) + sqr(y))});
        ARIADNE_TESTING_NOTE(vanderpol);

        RealExpressionBoundedConstraintSet initial_set({-2.0_dec <= x <= -1.5_dec, 0.5_dec <= y <= 1});

        Real time = 6.3_dec;

        ContinuousTimeEvolver<ApproximateTag> evolver(vanderpol);
        evolver.configuration().set_step_size(0.05);
        evolver.configuration().set_sampling_count(x, 1);
        evolver.configuration().set_sampling_count(y, 1);

        ARIADNE_TESTING_NOTE(evolver.configuration());

        Orbit<ListPointType> orbit = evolver.orbit(initial_set, time);

        LabelledFigure fig({-2.5_dec <= x <= 2.5_dec, -2.5_dec <= y <= 2.5_dec});
        fig << orbit;
        fig.write("test_real_expr_continuous_time_evolver_xy");

        LabelledFigure fig2({0 <= TimeVariable() <= 6.3_dec, 0 <= z <= 2.5_dec});
        fig2 << orbit;
        fig2.write("test_real_expr_continuous_time_evolver_tz");
    }

    Void test_multiple_real_box_vdp_cycle() {
        typedef ContinuousTimeEvolver<ApproximateTag>::ApproximateListPointType ListPointType;

        Real mu = Dyadic(0.5_x);
        RealVariable x("x"), y("y"), z("z");

        OrdinaryDifferentialAlgebraicEquation vanderpol({dot(x) = y, dot(y) = mu * (1 - x * x) * y - x}, {let(z) = sqrt(sqr(x) + sqr(y))});
        ARIADNE_TESTING_NOTE(vanderpol);

        RealVariablesBox initial_set({-2.0_dec <= x <= -1.5_dec, 0.5_dec <= y <= 1});

        Real time = 6.3_dec;

        ContinuousTimeEvolver<ApproximateTag> evolver(vanderpol);
        evolver.configuration().set_step_size(0.05);
        evolver.configuration().set_sampling_count(x, 1);
        evolver.configuration().set_sampling_count(y, 1);

        ARIADNE_TESTING_NOTE(evolver.configuration());

        Orbit<ListPointType> orbit = evolver.orbit(initial_set, time);

        LabelledFigure fig({-2.5_dec <= x <= 2.5_dec, -3 <= y <= 3});
        fig << orbit;
        fig.write("test_real_box_list_continuous_time_evolver_xy");

        LabelledFigure fig2({0 <= TimeVariable() <= 6.3_dec, 0 <= z <= 2.5_dec});
        fig2 << orbit;
        fig2.write("test_real_box_list_continuous_time_evolver_tz");
    }
};

int main(int argc, const char* argv[]) {
    ARIADNE_TESTING_SUITE("test_continuous_time_approximate_evolver");
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeApproximateEvolver, test_linear_flow_and_input_labels);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeApproximateEvolver, test_configuration_and_sampling);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeApproximateEvolver, test_multiple_real_expr_vdp_cycle);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeApproximateEvolver, test_multiple_real_box_vdp_cycle);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeApproximateEvolver, test_multiple_trajectories_x);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeApproximateEvolver, test_multiple_trajectories_xy);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
