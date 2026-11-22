/***************************************************************************
 *            test_ordinary_dynamics_contract.cpp
 *
 *  Copyright  2026  Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/

#include <stdexcept>
#include <type_traits>

#include "evolution/continuous_time_approximate_evolver.hpp"
#include "evolution/continuous_time_evolver.hpp"
#include "evolution/traits.hpp"
#include "function/taylor_model.hpp"
#include "integration/integrator.hpp"
#include "logging/test_runner.hpp"
#include "system/odae.hpp"
#include "system/odai.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

static_assert(std::is_same_v<ContinuousTimeEvolver<ValidatedTag>::SystemType, ContinuousDynamics>);
static_assert(std::is_same_v<ContinuousTimeEvolver<ApproximateTag>::SystemType, OrdinaryDifferentialAlgebraicEquation>);
static_assert(std::is_same_v<EvolutionTraits<ValidatedTag, OrdinaryDifferentialAlgebraicEquation>::EvolverType,
                             EvolutionTraits<ValidatedTag, OrdinaryDifferentialAlgebraicInclusion>::EvolverType>);
static_assert(std::is_same_v<EvolutionTraits<ValidatedTag, ContinuousTimeEvolver<ValidatedTag>::SystemType>::EnclosureType,
                             ContinuousTimeEvolver<ValidatedTag>::EnclosureType>);

class TestOrdinaryDynamicsContract : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_equation_consumers_enforce_numerical_capabilities() {
        RealVariable x("x"), t("t"), p("p");
        const OrdinaryDifferentialAlgebraicEquation timed({dot(x) = t - x}, RealVariablesBox{}, t);
        const OrdinaryDifferentialAlgebraicEquation parameterised({dot(x) = p - x}, RealVariablesBox({0 <= p <= 1}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-5);
        // Both explicit time and trajectory-constant parameters are supported.
        ARIADNE_TESTING_EXPECT_EQUAL(timed.dimension(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(timed.function().argument_size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(parameterised.function().argument_size(), 2u);
        ContinuousTimeEvolver<ValidatedTag> timed_evolver(timed, integrator);
        ARIADNE_TESTING_PROVE(timed_evolver.system().has_explicit_time());
        ARIADNE_TESTING_EXPECT_EQUAL(timed.auxiliary_space().size(), 0u);
        // Time-aware auxiliary outputs use the separate clock argument.
        RealVariable output("output");
        const OrdinaryDifferentialAlgebraicEquation timed_auxiliary(
            {dot(x) = t - x}, RealVariablesBox{}, t, {let(output) = t + x});
        ContinuousTimeEvolver<ValidatedTag> timed_output(timed_auxiliary, integrator);
        ARIADNE_TESTING_PROVE(timed_output.system().has_explicit_time());
        ContinuousTimeEvolver<ValidatedTag> parameterised_evolver(parameterised, integrator);
        ARIADNE_TESTING_EXPECT_EQUAL(parameterised_evolver.system().parameter_space().dimension(), 1u);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      parameterised_evolver.orbit(RealVariablesBox({x == 0}), Real(0.125_x), Semantics::LOWER));
        const OrdinaryDifferentialAlgebraicEquation parameterised_auxiliary(
            {dot(x) = p - x}, RealVariablesBox({0 <= p <= 1}), std::nullopt, {let(output) = p + x});
        ContinuousTimeEvolver<ValidatedTag> parameterised_output(parameterised_auxiliary, integrator);
        ARIADNE_TESTING_EXPECT_EQUAL(parameterised_output.system().auxiliary_space(), RealSpace({output}));
        PreconditionedGradedTaylorSeriesIntegrator preconditioned(step_maximum_error = 1e-5);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      ContinuousTimeEvolver<ValidatedTag>(parameterised, preconditioned));
        ContinuousTimeEvolver<ApproximateTag> timed_simulator(timed);
        ARIADNE_TESTING_PROVE(timed_simulator.configuration().step_size().raw().get_d() > 0);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, (void)ContinuousTimeEvolver<ApproximateTag>(parameterised));
    }

    Void test_inclusion_consumer_rejects_unsupported_arguments() {
        RealVariable x("x"), t("t"), p("p"), u("u");
        const RealVariablesBox inputs({-1 <= u <= 1});
        const OrdinaryDifferentialAlgebraicInclusion timed({dot(x) = t + u}, inputs, RealVariablesBox{}, t);
        const OrdinaryDifferentialAlgebraicInclusion parameterised({dot(x) = p + u}, inputs, RealVariablesBox({0 <= p <= 1}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-5);
        for (const auto& inclusion: {timed, parameterised}) {
            ARIADNE_TESTING_EXPECT_EQUAL(inclusion.dimension(), 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(inclusion.function().argument_size(), 3u);
            ContinuousTimeEvolver<ValidatedTag> evolver(inclusion, integrator, ZeroApproximation(), ConstantApproximation());
            ARIADNE_TESTING_EXPECT_EQUAL(evolver.system().dimension(), 1u);
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                          evolver.orbit(RealVariablesBox({x == 0}), Real(0.0625_x), Semantics::LOWER));
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
        ARIADNE_TESTING_PROVE(!inclusion.has_input_function());
        ARIADNE_TESTING_EXPECT_EQUAL(inclusion.number_of_inputs(), 0u);
        ContinuousTimeEvolver<ValidatedTag> evolver(inclusion, integrator, ZeroApproximation(), ConstantApproximation());
        ARIADNE_TESTING_EXPECT_EQUAL(evolver.system().state_space(), RealSpace({x}));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_ordinary_dynamics_contract");
    ARIADNE_TESTING_ADD_TEST(TestOrdinaryDynamicsContract, test_equation_consumers_enforce_numerical_capabilities);
    ARIADNE_TESTING_ADD_TEST(TestOrdinaryDynamicsContract, test_inclusion_consumer_rejects_unsupported_arguments);
    ARIADNE_TESTING_ADD_TEST(TestOrdinaryDynamicsContract, test_interval_inclusion_is_not_treated_as_an_input_function);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
