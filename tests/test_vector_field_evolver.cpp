/***************************************************************************
 *            test_vector_field_evolver.cpp
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

/** @file tests/test_vector_field_evolver.cpp
 *  @brief Test vector field evolver behavior and numerical invariants.
 */

#include <fstream>
#include <iostream>
#include <sstream>

#include "../src/carried_expansion_diagnostics.hpp"
#include "../src/evolution_time.hpp"
#include "algebra/algebra.hpp"
#include "algebra/matrix.hpp"
#include "algebra/vector.hpp"
#include "evolution/enclosure.hpp"
#include "evolution/orbit.hpp"
#include "evolution/vector_field_evolver.hpp"
#include "function/constraint.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "geometry/box.hpp"
#include "integration/integrator.hpp"
#include "labelled/figure.hpp"
#include "log_verbosity_scope.hpp"
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

class TestVectorFieldEvolver : public Ariadne::Testing::ContextFixture {
  public:
    using Ariadne::Testing::ContextFixture::ContextFixture;

  public:
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
        VectorField field(List<DottedRealAssignment>{dot(x) = 0});
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
        VectorFieldEvolver evolver(field, integrator);
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
    Void test_carried_expansion_profile_buckets() {
        Sweeper<FloatDP> sweeper = ThresholdSweeper<FloatDP>(dp, 0.0);
        auto coordinate = ValidatedTaylorModelDP::coordinate(1u, 0u, sweeper);
        auto model = coordinate + pow(coordinate, 3) * FloatDPBounds(ExactDouble(1e-13), dp) + pow(coordinate, 32) * FloatDPBounds(ExactDouble(1e-14), dp);
        ExactBoxType domain = ExactBoxType::unit_box(1u);
        ValidatedVectorMultivariateTaylorFunctionModelDP mapping(
            domain, Vector<ValidatedTaylorModelDP>{model});
        std::ostringstream text;
        EvolutionDetail::print_carried_expansion_snapshot(text, "test", 42u, mapping, sweeper);
        const auto result = text.str();
        ARIADNE_TESTING_PROVE(result.find("stage=test step=42") != std::string::npos);
        ARIADNE_TESTING_PROVE(result.find("degree=31 below_3e-14_count=1") != std::string::npos);
        ARIADNE_TESTING_PROVE(result.find("degree=3 below_3e-14_count=0 bridge_3e-14_to_1e-12_count=1") != std::string::npos);
        ARIADNE_TESTING_PROVE(result.find("above_1e-12_count=1") != std::string::npos);
        ValidatedVectorMultivariateTaylorFunctionModelDP empty(
            domain, Vector<ValidatedTaylorModelDP>(0u, model));
        std::ostringstream empty_text;
        EvolutionDetail::print_carried_expansion_snapshot(empty_text, "empty", 0u, empty, sweeper);
        ARIADNE_TESTING_PROVE(empty_text.str().empty());
    }
    Void test_components_and_preconditioned_orbits() {
        LogVerbosityScope logging(100u);
        RealVariable x("x"), other("other");
        VectorField field(List<DottedRealAssignment>{dot(x) = 1});
        for (auto mode: {TaylorSeriesPreconditioning::IDENTITY, TaylorSeriesPreconditioning::QR}) {
            PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error = 1e-6);
            integrator.set_preconditioning(mode);
            integrator.set_diagnostics(true);
            VectorFieldEvolver evolver(field, integrator);
            evolver.configuration().set_maximum_step_size(0.03125);
            evolver.configuration().set_enable_reconditioning(false);
            auto initial = evolver.enclosure(ExactBoxType({
                {-0.03125_x, 0.03125_x}
            }));
            std::unique_ptr<VectorFieldEvolver> clone(evolver.clone());
            const VectorFieldEvolver& immutable = *clone;
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
        VectorField field(List<DottedRealAssignment>{dot(x) = 0});
        PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error = 1e-6);
        integrator.set_carried_expansion_diagnostics(true);
        VectorFieldEvolver evolver(field, integrator);
        evolver.configuration().set_maximum_step_size(1.0 / 2048.0);
        auto result = evolver.orbit(RealVariablesBox({0.9_dec <= x <= 1.1_dec}), Real(1), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(result.reach().size(), 2048u);
        ARIADNE_TESTING_PROVE(subset(result.final().bounding_box().euclidean_set(),
                                     RealBox({
                                         {0.89_dec, 1.11_dec}
        })));
    }

    Void test_single_trajectory() {
        typedef VectorFieldEvolver::EnclosureType EnclosureType;

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

        VectorField vanderpol({dot(x) = v, dot(v) = mu * (1 - x * x) * v - x});
        ARIADNE_TESTING_NOTE(vanderpol);

        Semantics semantics = Semantics::LOWER;

        VectorFieldEvolver evolver(vanderpol, integrator);
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

        fig.write("test_vector_field_evolver-vdp");
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

        VectorField fail_vf({dot(x) = 1, dot(y) = y * y * 100});
        VectorFieldEvolver evolver(fail_vf, integrator);
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

        VectorField dynamics({dot(x) = y, dot(y) = mu * y * (1 - sqr(x)) - x});

        StepMaximumError max_err = 1e-8;

        TaylorPicardIntegrator integrator(max_err);

        VectorFieldEvolver evolver(dynamics, integrator);
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

        VectorField dynamics({dot(x) = 1, dot(y) = sqr(x)});

        StepMaximumError max_err = 1e-8;

        TaylorPicardIntegrator integrator(max_err);

        VectorFieldEvolver evolver(dynamics, integrator);
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
        fig.write("test_vector_field_evolver-subdivisions");

        auto orbit_lower = evolver.orbit(initial_set, evolution_time, Semantics::LOWER);

        ARIADNE_TESTING_PROVE(orbit_upper.reach().size() > orbit_lower.reach().size());
    }
};

int main(int argc, const char* argv[]) {
    ARIADNE_TESTING_SUITE("test_vector_field_evolver");
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldEvolver, test_single_trajectory);
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldEvolver, test_failure);
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldEvolver, test_subdivide_initially);
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldEvolver, test_subdivide_along_evolution);
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldEvolver, test_components_and_preconditioned_orbits);
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldEvolver, test_carried_expansion_diagnostics_checkpoints);
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldEvolver, test_carried_expansion_profile_buckets);
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldEvolver, test_constrained_initial_subdivision_and_tiny_horizon);
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldEvolver, test_binary64_step_limits);
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldEvolver, test_lower_termination_diagnostics);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
