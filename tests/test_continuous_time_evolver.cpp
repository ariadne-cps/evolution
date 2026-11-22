/***************************************************************************
 *            test_continuous_time_validated_evolver.cpp
 *
 *  Copyright  2006-21  Pieter Collins
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

/** @file tests/test_continuous_time_validated_evolver.cpp
 *  @brief Test vector field evolver behavior and numerical invariants.
 */

#include <fstream>
#include <iostream>

#include "../src/evolution_time.hpp"
#include "algebra/algebra.hpp"
#include "algebra/matrix.hpp"
#include "algebra/vector.hpp"
#include "evolution/continuous_time_evolver.hpp"
#include "evolution/enclosure.hpp"
#include "evolution/orbit.hpp"
#include "function/constraint.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
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

template<class IVL1, class IVL2>
decltype(auto) subset(LabelledBox<IVL1> bx1, LabelledBox<IVL2> bx2) {
    return subset(bx1.euclidean_set(), bx2.euclidean_set());
}

class TestContinuousTimeEvolver : public Ariadne::Testing::ContextFixture {
  public:
    Void test_autonomous_constant_parameters_preserve_the_endpoint_map() {
        RealVariable x("x"), y("y"), p("p"), q("q");
        const OrdinaryDifferentialAlgebraicEquation equation(
            {dot(x) = p, dot(y) = q * x},
            RealVariablesBox({-1 <= p <= 1, 0.5_dec <= q <= 1}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-8);
        ContinuousTimeEvolver<ValidatedTag> evolver(equation, integrator);
        evolver.configuration().set_enable_reconditioning(false);
        evolver.configuration().set_maximum_step_size(0.0625);
        const auto result = evolver.orbit(RealVariablesBox({x == 0, y == 0}),
                                          Real(0.09375_x), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(result.reach().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.final().size(), 1u);
        const auto& final = result.final()[0u];
        ARIADNE_TESTING_EXPECT_EQUAL(final.state_dimension(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(final.number_of_parameters(), 4u);
        const auto P = EffectiveScalarMultivariateFunction::coordinate(4u, 2u);
        const auto Q = EffectiveScalarMultivariateFunction::coordinate(4u, 3u);
        const EffectiveVectorMultivariateFunction expected = {
            P * Real(0.09375_x), P * Q * Real(0.00439453125_x)};
        ARIADNE_TESTING_EXPECT_LESS(norm(final.state_function() - expected), 1e-7_pr);
        // Evaluate opposite parameter choices in the SAME map. A box restart
        // or a swapped parameter order would lose this dependence.
        const Vector<FloatDPBounds> positive({FloatDPBounds(0, dp), FloatDPBounds(0, dp),
                                              FloatDPBounds(1, dp), FloatDPBounds(0.75_x, dp)});
        const Vector<FloatDPBounds> negative({FloatDPBounds(0, dp), FloatDPBounds(0, dp),
                                              FloatDPBounds(-1, dp), FloatDPBounds(0.75_x, dp)});
        const auto delta = final.state_function()(positive) - final.state_function()(negative);
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.1875_x, dp), delta[0u]));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.006591796875_x, dp), delta[1u]));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.09375_x, dp), final.time_function()(positive)));
    }

    Void test_constant_parameters_survive_reconditioning_and_subdivision() {
        // Integral_0^H p*(2*t-H) dt = 0 for every FIXED p. Reselecting p
        // on the second step yields a much wider interval and fails this test.
        RealVariable x("x"), t("t"), p("p");
        const OrdinaryDifferentialAlgebraicEquation equation(
            {dot(x) = p * (2 * t - Real(0.09375_x))}, RealVariablesBox({-1 <= p <= 1}), t);
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-8);
        ContinuousTimeEvolver<ValidatedTag> evolver(equation, integrator);
        evolver.configuration().set_maximum_step_size(0.0625);
        evolver.configuration().set_maximum_spacial_error(0x1p-18);
        evolver.configuration().set_enable_subdivisions(true);
        evolver.configuration().set_maximum_enclosure_radius(0.0005);
        const ExactBoxType domain({
            {0.0_x, 0.0_x}
        });
        auto model = ValidatedVectorMultivariateTaylorFunctionModelDP::identity(
            domain, ThresholdSweeper<FloatDP>(dp, 1e-10));
        model[0u].set_error(FloatDPError(0x1p-16_x, dp));
        EnclosureConfiguration config(integrator.function_factory(), 1u);
        const LabelledEnclosure initial(Enclosure(domain, model, config), RealSpace({x}));
        const auto result = evolver.orbit(initial, Real(0.09375_x), Semantics::UPPER);
        ARIADNE_TESTING_PROVE(result.final().size() > 1u);
        const auto final_box = cast_exact_box(result.final().bounding_box().euclidean_set());
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {-0x1p-16_x, 0x1p-16_x}
        }),
                                     final_box));
        ARIADNE_TESTING_PROVE(subset(final_box, ExactBoxType({
                                                    {-0x1p-15_x, 0x1p-15_x}
        })));
        for (const auto& final: result.final()) {
            ARIADNE_TESTING_EXPECT_EQUAL(final.state_dimension(), 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(final.variable_kinds()[0u], EnclosureVariableKind::PARAMETER);
            ARIADNE_TESTING_EXPECT_EQUAL(final.variable_kinds()[1u], EnclosureVariableKind::ERROR);
            ARIADNE_TESTING_EXPECT_EQUAL(final.number_of_parameters(), 2u);
        }
    }

    using Ariadne::Testing::ContextFixture::ContextFixture;

  public:
    Void test_non_autonomous_two_step_orbit() {
        // The same numerical time coordinate must advance across consecutive
        // steps: x'=t, x(0)=0 has x(1/8)=1/128.
        RealVariable x("x"), t("t");
        const OrdinaryDifferentialAlgebraicEquation equation(
            {dot(x) = t}, RealVariablesBox{}, t);
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
        ContinuousTimeEvolver<ValidatedTag> evolver(equation, integrator);
        evolver.configuration().set_enable_reconditioning(false);
        evolver.configuration().set_maximum_step_size(0.0625);
        const auto result = evolver.orbit(RealVariablesBox({x == 0}),
                                          Real(0.125_x), Semantics::UPPER);
        ARIADNE_TESTING_PROVE(result.reach().size() >= 2u);
        ARIADNE_TESTING_PROVE(result.final().size() > 0u);
        const auto final_box = result.final().bounding_box().euclidean_set();
        // Test the certified box's outward-rounded binary endpoints as exact
        // values. Containment in an UpperBox is only upper semidecidable and
        // cannot prove this inclusion, even when the numerical bounds fit.
        const ExactBoxType exact_endpoint({
            {0.0078125_x, 0.0078125_x}
        });
        ARIADNE_TESTING_PROVE(subset(exact_endpoint, cast_exact_box(final_box)));
        const RealBox expected({
            {0.007_dec, 0.009_dec}
        });
        if (not definitely(subset(final_box, expected))) {
            // Diagnostic only: keep the tight numerical assertion unchanged.
            // Expose the certified range of each step to distinguish a
            // widened validated enclosure from a time-coordinate error.
            std::cerr << "[NonautonomousTwoStepOrbit] final=" << final_box
                      << " reach_count=" << result.reach().size()
                      << " intermediate_count=" << result.intermediate().size()
                      << std::endl;
            SizeType i = 0u;
            for (const auto& stage: result.intermediate()) {
                if (i >= 8u) {
                    break;
                }
                std::cerr << "[NonautonomousTwoStepOrbit] endpoint=" << i
                          << " box=" << stage.euclidean_set().bounding_box()
                          << " absolute_time=" << stage.time_function().range()
                          << std::endl;
                ++i;
            }
        }
        ARIADNE_TESTING_PROVE(subset(final_box, expected));
    }

    Void test_binary64_step_limits() {
        using EvolutionDetail::limit_step_to_remaining_time;
        const Dyadic quarter(0.25_x), half(0.5_x);
        ARIADNE_TESTING_EXPECT_EQUAL(limit_step_to_remaining_time(half, quarter), quarter);
        ARIADNE_TESTING_EXPECT_EQUAL(limit_step_to_remaining_time(quarter, half), quarter);
        const Dyadic smallest(Integer(1), std::size_t(1074u));
        ARIADNE_TESTING_EXPECT_EQUAL(limit_step_to_remaining_time(half, smallest), smallest);
        const Dyadic unrepresentable(Integer(1), std::size_t(1100u));
        ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error,
                                      limit_step_to_remaining_time(half, unrepresentable));
        const Dyadic between(Integer(9007199254740993ULL), std::size_t(54u));
        ARIADNE_TESTING_EXPECT_EQUAL(limit_step_to_remaining_time(Dyadic(1u), between), half);
    }

    Void test_lower_termination_diagnostics() {
        LogVerbosityScope logging(100u);
        test_subdivide_along_evolution();
    }
    Void test_constrained_initial_subdivision_and_tiny_horizon() {
        LogVerbosityScope logging(100u);
        RealVariable x("x");
        OrdinaryDifferentialAlgebraicEquation field(List<DottedRealAssignment>{dot(x) = 0});
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
        ContinuousTimeEvolver<ValidatedTag> evolver(field, integrator);
        evolver.configuration().set_enable_reconditioning(false);
        evolver.configuration().set_maximum_enclosure_radius(0.75);
        for (const auto& constraint: List<ContinuousPredicate>{x <= 0.25_dec, x >= 1.75_dec}) {
            RealExpressionBoundedConstraintSet initial({0 <= x <= 2}, {constraint});
            auto result = evolver.orbit(initial, Real(0), Semantics::UPPER);
            ARIADNE_TESTING_EXPECT_EQUAL(result.final().size(), 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(result.reach().size(), 0u);
        }
        auto result = evolver.orbit(RealVariablesBox({0 <= x <= 0}),
                                    rec(pow(Real(2), 1100u)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(result.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.reach().size(), 0u);
    }
    Void test_components_and_preconditioned_orbits() {
        LogVerbosityScope logging(100u);
        RealVariable x("x"), other("other");
        OrdinaryDifferentialAlgebraicEquation field(List<DottedRealAssignment>{dot(x) = 1});
        for (auto mode: {TaylorSeriesPreconditioning::IDENTITY, TaylorSeriesPreconditioning::QR}) {
            PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error = 1e-6);
            integrator.set_preconditioning(mode);
            integrator.set_diagnostics(true);
            ContinuousTimeEvolver<ValidatedTag> evolver(field, integrator);
            evolver.configuration().set_maximum_step_size(0.03125);
            evolver.configuration().set_enable_reconditioning(false);
            auto initial = evolver.enclosure(ExactBoxType({
                {-0.03125_x, 0.03125_x}
            }));
            std::unique_ptr<ContinuousTimeEvolver<ValidatedTag>> clone(evolver.clone());
            const ContinuousTimeEvolver<ValidatedTag>& immutable = *clone;
            ARIADNE_TESTING_EXPECT_EQUAL(immutable.configuration().maximum_step_size(), 0.03125_x);
            auto result = immutable.orbit(initial, Real(1), Semantics::UPPER);
            ARIADNE_TESTING_PROVE(result.final().size() > 0u);
            auto box = result.final().bounding_box().euclidean_set();
            ARIADNE_TESTING_PROVE(subset(box, RealBox({
                                                  {0.96_dec, 1.04_dec}
            })));
            auto wrong = LabelledEnclosure(ExactBoxType({
                                               {-0.03125_x, 0.03125_x}
            }),
                                           RealSpace({other}), EnclosureConfiguration(integrator.function_factory()));
            ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.orbit(wrong, Real(1), Semantics::UPPER));
        }
    }

    Void test_carried_expansion_diagnostics_checkpoints() {
        RealVariable x("x");
        OrdinaryDifferentialAlgebraicEquation field(List<DottedRealAssignment>{dot(x) = 0});
        PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error = 1e-6);
        integrator.set_carried_expansion_diagnostics(true);
        ContinuousTimeEvolver<ValidatedTag> evolver(field, integrator);
        evolver.configuration().set_maximum_step_size(1.0 / 2048.0);
        auto result = evolver.orbit(RealVariablesBox({0.9_dec <= x <= 1.1_dec}), Real(1), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(result.reach().size(), 2048u);
        ARIADNE_TESTING_PROVE(subset(result.final().bounding_box().euclidean_set(),
                                     RealBox({
                                         {0.89_dec, 1.11_dec}
        })));
    }

    Void test_single_trajectory() {
        typedef ContinuousTimeEvolver<ValidatedTag>::EnclosureType EnclosureType;

        // Set up the evolution parameters and grid
        Real time(2.0_dec);
        ExactDouble step_size(0.5_x);
        ExactDouble enclosure_radius(0.25_x);

        ThresholdSweeper<FloatDP> sweeper(DoublePrecision(), 1e-8_pr);

        // Set up the evaluators
        TaylorPicardIntegrator picard_integrator(step_maximum_error = 1e-7_pr, sweeper, lipschitz_tolerance = 0.5_x,
                                                 minimum_temporal_order = 0, maximum_temporal_order = 8);
        // Set up the evaluators
        GradedTaylorSeriesIntegrator series_integrator(step_maximum_error = 1e-7_pr, sweeper, lipschitz_tolerance = 0.5_x,
                                                       minimum_spacial_order = 1, minimum_temporal_order = 4,
                                                       maximum_spacial_order = 3, maximum_temporal_order = 8);

        IntegratorInterface& integrator = picard_integrator;

        ARIADNE_TESTING_NOTE(integrator);

        // Set up the vector field
        Real mu = Dyadic(0.5_x);
        RealVariable x("x"), v("v");

        OrdinaryDifferentialAlgebraicEquation vanderpol({dot(x) = v, dot(v) = mu * (1 - x * x) * v - x});
        ARIADNE_TESTING_NOTE(vanderpol);

        Semantics semantics = Semantics::LOWER;

        ContinuousTimeEvolver<ValidatedTag> evolver(vanderpol, integrator);
        evolver.configuration().set_maximum_enclosure_radius(enclosure_radius);
        evolver.configuration().set_maximum_step_size(step_size);
        ARIADNE_TESTING_NOTE(evolver.configuration());

        LabelledFigure fig(Axes2d(-2.25_dec <= x <= +2.25_dec, -2.25_dec <= v <= +2.25_dec));

        // Define the initial set
        RealExpressionBoundedConstraintSet initial_set({0.95_dec <= x <= 1.05_dec, 0.45_dec <= v <= 0.55_dec});
        Orbit<EnclosureType> flow_tube = evolver.orbit(initial_set, time, semantics);
        ARIADNE_TESTING_NOTE(flow_tube);

        // The requested horizon is t=2. Reference integration of the initial
        // rectangle gives approximately x in [-0.15902,0.00873] and
        // v in [-1.51529,-1.39682]. Check tight validated enclosures,
        // not the obsolete bounds from evolution beyond t=2.
        LabelledBox<RealInterval> expected_final_flow_tube_bounding_box({
                                                                            x, v
        },
                                                                        {{-0.165_dec, 0.015_dec}, {-1.522_dec, -1.390_dec}});
        ARIADNE_TESTING_PROVE(subset(flow_tube.final().bounding_box(), expected_final_flow_tube_bounding_box));

        fig << line_style(true) << fill_colour(cyan) << flow_tube.reach();
        fig << fill_colour(magenta) << flow_tube.intermediate();
        fig << fill_colour(red) << flow_tube.final();

        // Define the initial point
        RealExpressionBoundedConstraintSet initial_point({1.0_dec <= x <= 1.0_dec, 0.5_dec <= v <= 0.5_dec});
        Orbit<EnclosureType> trajectory = evolver.orbit(initial_point, time, semantics);

        // For the point (1,0.5), independent integration to t=2 yields
        // approximately (-0.076261,-1.456835); keep a narrow regression box.
        LabelledBox<RealInterval> expected_final_trajectory_bounding_box({
                                                                             x, v
        },
                                                                         {{-0.081_dec, -0.071_dec}, {-1.462_dec, -1.452_dec}});
        ARIADNE_TESTING_PROVE(subset(trajectory.final().bounding_box(), expected_final_trajectory_bounding_box));

        fig << line_style(true) << fill_colour(blue) << trajectory.reach();

        fig.write("test_continuous_time_evolver-vdp");
    }

    Void test_failure() {
        // The system in this test is stiff and is expected to fail.

        RealVariable x("x"), y("y");

        // Set up the evolution parameters and grid
        Real time(0.5_dec);
        ExactDouble maximum_step_size(0.015625_pr);
        ExactDouble enclosure_radius(0.25_x);

        ThresholdSweeper<FloatDP> sweeper(DoublePrecision(), 1e-8_pr);

        // Set up the evaluators
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6_pr, sweeper, lipschitz_tolerance = 0.5_x,
                                          minimum_temporal_order = 0, maximum_temporal_order = 8);

        OrdinaryDifferentialAlgebraicEquation fail_vf({dot(x) = 1, dot(y) = y * y * 100});
        ContinuousTimeEvolver<ValidatedTag> evolver(fail_vf, integrator);
        evolver.configuration().set_maximum_enclosure_radius(enclosure_radius);
        evolver.configuration().set_maximum_step_size(maximum_step_size);
        //        evolver.configuration().set_minimum_step_size(minimum_step_size);

        RealVariablesBox initial_box({x == 0, y == 1});

        time = 1.5_dec;

        // Compute the reachable sets; expect failure as required step size becomes too small
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.orbit(initial_box, time, Semantics::UPPER));
    }

    Void test_subdivide_initially() {
        RealConstant mu("mu", 1.0_dec);
        RealVariable x("x"), y("y");

        OrdinaryDifferentialAlgebraicEquation dynamics({dot(x) = y, dot(y) = mu * y * (1 - sqr(x)) - x});

        StepMaximumError max_err = 1e-8;

        TaylorPicardIntegrator integrator(max_err);

        ContinuousTimeEvolver<ValidatedTag> evolver(dynamics, integrator);
        evolver.configuration().set_maximum_enclosure_radius(0.1);
        evolver.configuration().set_maximum_step_size(0.1);
        evolver.configuration().set_maximum_spacial_error(1e-5);
        evolver.configuration().set_enable_subdivisions(true);

        Real x0(1.40_dec);
        Real y0(2.40_dec);
        Real eps_x0 = 15 / 100_q;
        Real eps_y0 = 5 / 100_q;

        RealExpressionBoundedConstraintSet initial_set({x0 - eps_x0 <= x <= x0 + eps_x0, y0 - eps_y0 <= y <= y0 + eps_y0});

        Real evolution_time(0.2_dec);

        auto evolution = evolver.orbit(initial_set, evolution_time, Semantics::UPPER);

        ARIADNE_TESTING_PROVE(evolution.final().size() > 1);
    }

    Void test_subdivide_along_evolution() {
        RealVariable x("x"), y("y");

        OrdinaryDifferentialAlgebraicEquation dynamics({dot(x) = 1, dot(y) = sqr(x)});

        StepMaximumError max_err = 1e-8;

        TaylorPicardIntegrator integrator(max_err);

        ContinuousTimeEvolver<ValidatedTag> evolver(dynamics, integrator);
        evolver.configuration().set_maximum_enclosure_radius(0.1);
        evolver.configuration().set_maximum_step_size(0.1);
        evolver.configuration().set_maximum_spacial_error(1e-5);
        evolver.configuration().set_enable_subdivisions(true);

        Real x0(1.40_dec);
        Real y0(2.40_dec);
        Real eps_x0 = 5 / 100_q;
        Real eps_y0 = 5 / 100_q;

        RealExpressionBoundedConstraintSet initial_set({x0 - eps_x0 <= x <= x0 + eps_x0, y0 - eps_y0 <= y <= y0 + eps_y0});

        Real evolution_time(3);

        auto orbit_upper = evolver.orbit(initial_set, evolution_time, Semantics::UPPER);

        LabelledFigure fig(Axes2d(1 <= x <= 5, 2 <= y <= 30));
        fig << orbit_upper.reach();
        fig.write("test_continuous_time_evolver-subdivisions");

        auto orbit_lower = evolver.orbit(initial_set, evolution_time, Semantics::LOWER);

        ARIADNE_TESTING_PROVE(orbit_upper.reach().size() > orbit_lower.reach().size());
    }
};

int main(int argc, const char* argv[]) {
    ARIADNE_TESTING_SUITE("test_continuous_time_evolver");
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeEvolver, test_single_trajectory);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeEvolver, test_failure);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeEvolver, test_subdivide_initially);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeEvolver, test_subdivide_along_evolution);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeEvolver, test_components_and_preconditioned_orbits);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeEvolver, test_carried_expansion_diagnostics_checkpoints);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeEvolver, test_constrained_initial_subdivision_and_tiny_horizon);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeEvolver, test_binary64_step_limits);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeEvolver, test_non_autonomous_two_step_orbit);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeEvolver, test_autonomous_constant_parameters_preserve_the_endpoint_map);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeEvolver, test_constant_parameters_survive_reconditioning_and_subdivision);
    ARIADNE_TESTING_ADD_TEST(TestContinuousTimeEvolver, test_lower_termination_diagnostics);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
