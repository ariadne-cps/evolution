/***************************************************************************
 *            test_ordinary_dynamics_contract.cpp
 *
 *  Copyright  2026  Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/

#include <stdexcept>
#include <type_traits>

#include "evolution/differential_inclusion_evolver.hpp"
#include "evolution/vector_field_evolver.hpp"
#include "evolution/vector_field_simulator.hpp"
#include "function/taylor_model.hpp"
#include "integration/integrator.hpp"
#include "logging/test_runner.hpp"
#include "system/odae.hpp"
#include "system/odai.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

static_assert(std::is_same_v<VectorFieldEvolver::SystemType, OrdinaryDifferentialAlgebraicEquation>);
static_assert(std::is_same_v<VectorFieldSimulator::SystemType, OrdinaryDifferentialAlgebraicEquation>);
static_assert(std::is_same_v<DifferentialInclusionEvolver::SystemType, OrdinaryDifferentialAlgebraicInclusion>);

class TestOrdinaryDynamicsContract : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_equation_consumers_reject_unsupported_arguments() {
        RealVariable x("x"), t("t"), p("p");
        const OrdinaryDifferentialAlgebraicEquation timed({dot(x) = t - x}, RealVariablesBox{}, t);
        const OrdinaryDifferentialAlgebraicEquation parameterised({dot(x) = p - x}, RealVariablesBox({0 <= p <= 1}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-5);
        for (const auto& equation: {timed, parameterised}) {
            ARIADNE_TESTING_EXPECT_EQUAL(equation.dimension(), 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(equation.function().argument_size(), 2u);
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, VectorFieldEvolver(equation, integrator));
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, (void)VectorFieldSimulator(equation));
        }
    }

    Void test_inclusion_consumer_rejects_unsupported_arguments() {
        RealVariable x("x"), t("t"), p("p"), u("u");
        const RealVariablesBox inputs({-1 <= u <= 1});
        const OrdinaryDifferentialAlgebraicInclusion timed({dot(x) = t + u}, inputs, RealVariablesBox{}, t);
        const OrdinaryDifferentialAlgebraicInclusion parameterised({dot(x) = p + u}, inputs, RealVariablesBox({0 <= p <= 1}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-5);
        LohnerReconditioner reconditioner(1u, 1u, 1u, 1.0_x);
        for (const auto& inclusion: {timed, parameterised}) {
            ARIADNE_TESTING_EXPECT_EQUAL(inclusion.dimension(), 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(inclusion.function().argument_size(), 3u);
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, DifferentialInclusionEvolver(inclusion, integrator, reconditioner));
        }
    }

    Void test_interval_inclusion_is_not_treated_as_an_input_function() {
        RealVariable x("x");
        const BoxDomainType domain({
            {-1.0_x, 1.0_x}
        });
        const ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        const auto uncertain = ValidatedIntervalTaylorFunctionModel<FloatDP>::constant(
            domain, FloatDPUpperInterval(-1, 1, dp), sweeper);
        const ValidatedVectorIntervalTaylorFunctionModel<FloatDP> field({uncertain});
        const OrdinaryDifferentialAlgebraicInclusion inclusion(RealSpace({x}), field);
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-5);
        LohnerReconditioner reconditioner(1u, 1u, 1u, 1.0_x);
        ARIADNE_TESTING_PROVE(!inclusion.has_input_function());
        ARIADNE_TESTING_EXPECT_EQUAL(inclusion.number_of_inputs(), 0u);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, DifferentialInclusionEvolver(inclusion, integrator, reconditioner));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_ordinary_dynamics_contract");
    ARIADNE_TESTING_ADD_TEST(TestOrdinaryDynamicsContract, test_equation_consumers_reject_unsupported_arguments);
    ARIADNE_TESTING_ADD_TEST(TestOrdinaryDynamicsContract, test_inclusion_consumer_rejects_unsupported_arguments);
    ARIADNE_TESTING_ADD_TEST(TestOrdinaryDynamicsContract, test_interval_inclusion_is_not_treated_as_an_input_function);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
