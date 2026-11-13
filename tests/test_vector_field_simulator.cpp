/***************************************************************************
 *            test_vector_field_simulator.cpp
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

/** @file tests/test_vector_field_simulator.cpp
 *  @brief Test vector field simulator behavior and numerical invariants.
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include "algebra/algebra.hpp"
#include "algebra/matrix.hpp"
#include "algebra/vector.hpp"
#include "evolution/enclosure.hpp"
#include "evolution/orbit.hpp"
#include "evolution/vector_field_simulator.hpp"
#include "function/constraint.hpp"
#include "function/function.hpp"
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

class TestVectorFieldSimulator : public Ariadne::Testing::ContextFixture {
  public:
    using Ariadne::Testing::ContextFixture::ContextFixture;

  public:
    Void test_linear_flow_and_input_labels() {
        RealVariable x("x"), y("y"), z("z"), missing("missing");
        VectorField field({dot(x) = 1, dot(y) = 0}, {let(z) = x + 2 * y});
        VectorFieldSimulator simulator(field);
        simulator.configuration().set_step_size(0.25);
        using PointType = VectorFieldSimulator::ApproximatePointType;
        using Points = VectorFieldSimulator::ApproximateListPointType;
        PointType start(RealSpace({x, y}), FloatDPApproximationPoint({0.0_x, 2.0_x}, dp));
        PointType reordered(RealSpace({y, x}), FloatDPApproximationPoint({2.0_x, 0.0_x}, dp));
        for (unsigned verbosity: {0u, 100u}) {
            LogVerbosityScope logging(verbosity);
            auto orbit = simulator.orbit(Points{start, reordered}, Real(1));
            ARIADNE_TESTING_EXPECT_EQUAL(orbit.curves().size(), 2u);
            for (const auto& curve: orbit.curves()) {
                const auto& endpoint = curve.end();
                ARIADNE_TESTING_EXPECT_EQUAL(endpoint->first.get_d(), 1.0);
                ARIADNE_TESTING_PROVE(std::abs(endpoint->second[curve.space()[x]].raw().get_d() - 1.0) < 1e-12);
                ARIADNE_TESTING_PROVE(std::abs(endpoint->second[curve.space()[y]].raw().get_d() - 2.0) < 1e-12);
                ARIADNE_TESTING_PROVE(std::abs(endpoint->second[curve.space()[z]].raw().get_d() - 5.0) < 1e-12);
            }
        }
        ARIADNE_TESTING_EXPECT_EQUAL(simulator.orbit(Points{}, Real(1)).curves().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(simulator.orbit(Points{start}, Real(0)).curve(0).size(), 1u);
        PointType incompatible(RealSpace({missing, y}), FloatDPApproximationPoint({0.0_x, 2.0_x}, dp));
        ARIADNE_TESTING_EXPECT_EXCEPTION(simulator.orbit(Points{incompatible}, Real(1)));
    }

    Void test_configuration_and_sampling() {
        RealVariable x("x"), y("y"), unknown("unknown");
        VectorField field({dot(x) = 0, dot(y) = 0});
        VectorFieldSimulator simulator(field);
        auto& config = simulator.configuration();
        ARIADNE_TESTING_EXPECT_EXCEPTION(config.set_step_size(0.0));
        ARIADNE_TESTING_EXPECT_EXCEPTION(config.set_step_size(-1.0));
        ARIADNE_TESTING_EXPECT_EXCEPTION(config.insert_subspace(x, 0u));
        config.insert_subspace(x, 2u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(config.insert_subspace(x, 3u));
        ARIADNE_TESTING_EXPECT_EXCEPTION(config.set_subspace(subspace{
            {x, 0u}
        }));
        ARIADNE_TESTING_EXPECT_EQUAL(config.get_subspace().size(), 1u);
        config.remove_keys();
        ARIADNE_TESTING_PROVE(!config.is_using_subspace());
        config.set_subspace(subspace{});
        ARIADNE_TESTING_PROVE(!config.is_using_subspace());
        config.set_subspace(subspace{
            {x, 2u}
        });
        config.set_num_subdivisions(3u);
        const auto& configured = simulator;
        ARIADNE_TESTING_PROVE(configured.configuration().is_using_subspace());
        UpperBoxType box({
            {0, 1},
            {0, 1}
        });
        auto divided = simulator.orbit(box, Real(0));
        ARIADNE_TESTING_EXPECT_EQUAL(config.num_subdivisions(), 3u);
        ARIADNE_TESTING_PROVE(divided.curves().size() > 0u);
        config.set_subspace(subspace{
            {unknown, 2u}
        });
        ARIADNE_TESTING_EXPECT_EXCEPTION(simulator.orbit(box, Real(0)));
        config.remove_keys();
        config.set_num_subdivisions(0u);
        for (auto type: {DiscretisationType::None, DiscretisationType::Recombine}) {
            config.set_discretisation_type(type);
            auto sampled = simulator.orbit(box, Real(0));
            ARIADNE_TESTING_PROVE(sampled.curves().size() > 0u);
            std::ostringstream description;
            const ConfigurationInterface& base = config;
            description << base;
            ARIADNE_TESTING_PROVE(description.str().find("VectorFieldSimulatorConfiguration") != std::string::npos);
        }
    }

    Void test_multiple_trajectories_x() {
        typedef VectorFieldSimulator::ApproximateListPointType ListPointType;

        Real mu = Dyadic(0.5_x);
        RealVariable x("x"), y("y"), z("z");

        VectorField vanderpol({dot(x) = y, dot(y) = mu * (1 - x * x) * y - x}, {let(z) = sqrt(sqr(x) + sqr(y))});
        ARIADNE_TESTING_NOTE(vanderpol);

        RealExpressionBoundedConstraintSet initial_set({-2.0_dec <= x <= -1.5_dec, 0.5_dec <= y <= 1});

        Real time = 1.0_dec;

        VectorFieldSimulator simulator(vanderpol);
        simulator.configuration().set_step_size(0.05);
        simulator.configuration().set_num_subdivisions(0);
        simulator.configuration().set_discretisation_type(DiscretisationType::Mince);

        Nat trajectory_number = 4;

        for (Nat t = 1; t <= trajectory_number; t++) {
            simulator.configuration().insert_subspace(x, t);
            ARIADNE_TESTING_NOTE(simulator.configuration());

            Orbit<ListPointType> orbit = simulator.orbit(initial_set, time);

            List<LabelledInterpolatedCurve> curves = orbit.curves();

            ARIADNE_TESTING_EXPECT_EQUAL(curves.size(), t);

            simulator.configuration().remove_keys();
        }
    }

    Void test_multiple_trajectories_xy() {
        typedef VectorFieldSimulator::ApproximateListPointType ListPointType;

        Real mu = Dyadic(0.5_x);
        RealVariable x("x"), y("y"), z("z");

        VectorField vanderpol({dot(x) = y, dot(y) = mu * (1 - x * x) * y - x}, {let(z) = sqrt(sqr(x) + sqr(y))});
        ARIADNE_TESTING_NOTE(vanderpol);

        RealExpressionBoundedConstraintSet initial_set({-2.0_dec <= x <= -1.5_dec, 0.5_dec <= y <= 1});

        Real time = 1.0_dec;

        VectorFieldSimulator simulator(vanderpol);
        simulator.configuration().set_step_size(0.05);
        simulator.configuration().set_num_subdivisions(0);
        simulator.configuration().set_discretisation_type(DiscretisationType::Mince);

        Nat trajectory_number = 4;

        for (Nat first = 1; first <= trajectory_number; first++) {
            for (Nat second = 1; second <= trajectory_number; second++) {
                simulator.configuration().insert_subspace(x, first);
                simulator.configuration().insert_subspace(y, second);
                ARIADNE_TESTING_NOTE(simulator.configuration());

                Orbit<ListPointType> orbit = simulator.orbit(initial_set, time);

                List<LabelledInterpolatedCurve> curves = orbit.curves();

                ARIADNE_TESTING_EXPECT_EQUAL(curves.size(), first * second);

                simulator.configuration().remove_keys();
            }
        }
    }

    Void test_multiple_real_expr_vdp_cycle() {
        typedef VectorFieldSimulator::ApproximateListPointType ListPointType;

        Real mu = Dyadic(0.5_x);
        RealVariable x("x"), y("y"), z("z");

        VectorField vanderpol({dot(x) = y, dot(y) = mu * (1 - x * x) * y - x}, {let(z) = sqrt(sqr(x) + sqr(y))});
        ARIADNE_TESTING_NOTE(vanderpol);

        RealExpressionBoundedConstraintSet initial_set({-2.0_dec <= x <= -1.5_dec, 0.5_dec <= y <= 1});

        Real time = 6.3_dec;

        VectorFieldSimulator simulator(vanderpol);
        simulator.configuration().set_step_size(0.05);
        simulator.configuration().set_num_subdivisions(0);
        simulator.configuration().set_discretisation_type(DiscretisationType::Mince);
        simulator.configuration().insert_subspace(x, 1);
        simulator.configuration().insert_subspace(y, 1);

        ARIADNE_TESTING_NOTE(simulator.configuration());

        Orbit<ListPointType> orbit = simulator.orbit(initial_set, time);

        LabelledFigure fig({-2.5_dec <= x <= 2.5_dec, -2.5_dec <= y <= 2.5_dec});
        fig << orbit;
        fig.write("test_real_expr_vector_field_simulator_xy");

        LabelledFigure fig2({0 <= TimeVariable() <= 6.3_dec, 0 <= z <= 2.5_dec});
        fig2 << orbit;
        fig2.write("test_real_expr_vector_field_simulator_tz");
    }

    Void test_multiple_real_box_vdp_cycle() {
        typedef VectorFieldSimulator::ApproximateListPointType ListPointType;

        Real mu = Dyadic(0.5_x);
        RealVariable x("x"), y("y"), z("z");

        VectorField vanderpol({dot(x) = y, dot(y) = mu * (1 - x * x) * y - x}, {let(z) = sqrt(sqr(x) + sqr(y))});
        ARIADNE_TESTING_NOTE(vanderpol);

        RealVariablesBox initial_set({-2.0_dec <= x <= -1.5_dec, 0.5_dec <= y <= 1});

        Real time = 6.3_dec;

        VectorFieldSimulator simulator(vanderpol);
        simulator.configuration().set_step_size(0.05);
        simulator.configuration().set_num_subdivisions(0);
        simulator.configuration().set_discretisation_type(DiscretisationType::Mince);
        simulator.configuration().insert_subspace(x, 1);
        simulator.configuration().insert_subspace(y, 1);

        ARIADNE_TESTING_NOTE(simulator.configuration());

        Orbit<ListPointType> orbit = simulator.orbit(initial_set, time);

        LabelledFigure fig({-2.5_dec <= x <= 2.5_dec, -3 <= y <= 3});
        fig << orbit;
        fig.write("test_real_box_list_vector_field_simulator_xy");

        LabelledFigure fig2({0 <= TimeVariable() <= 6.3_dec, 0 <= z <= 2.5_dec});
        fig2 << orbit;
        fig2.write("test_real_box_list_vector_field_simulator_tz");
    }
};

int main(int argc, const char* argv[]) {
    ARIADNE_TESTING_SUITE("test_vector_field_simulator");
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldSimulator, test_linear_flow_and_input_labels);
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldSimulator, test_configuration_and_sampling);
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldSimulator, test_multiple_real_expr_vdp_cycle);
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldSimulator, test_multiple_real_box_vdp_cycle);
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldSimulator, test_multiple_trajectories_x);
    ARIADNE_TESTING_ADD_TEST(TestVectorFieldSimulator, test_multiple_trajectories_xy);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
