/***************************************************************************
 *            test_ordinary_dynamics_contract.cpp
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

/** @file tests/test_ordinary_dynamics_contract.cpp
 *  @brief Verify ordinary dynamics contract behavior.
 */

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
static_assert(std::is_same_v<ContinuousTimeEvolver<ApproximateTag>::SystemType, ContinuousDynamics>);
static_assert(std::is_same_v<EvolutionTraits<ApproximateTag, OrdinaryDifferentialAlgebraicEquation>::EvolverType,
                             EvolutionTraits<ApproximateTag, OrdinaryDifferentialAlgebraicInclusion>::EvolverType>);
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
        ContinuousTimeEvolver<ValidatedTag> validated_timed_continuous_evolver(timed, integrator);
        ARIADNE_TESTING_PROVE(validated_timed_continuous_evolver.system().has_explicit_time());
        ARIADNE_TESTING_EXPECT_EQUAL(timed.auxiliary_space().size(), 0u);
        // Time-aware auxiliary outputs use the separate clock argument.
        RealVariable output("output");
        const OrdinaryDifferentialAlgebraicEquation timed_auxiliary(
            {dot(x) = t - x}, RealVariablesBox{}, t, {let(output) = t + x});
        ContinuousTimeEvolver<ValidatedTag> timed_output(timed_auxiliary, integrator);
        ARIADNE_TESTING_PROVE(timed_output.system().has_explicit_time());
        // Resolve the analytic residual below more tightly than the broad
        // capability checks, which use a 1e-5 local integration tolerance.
        ContinuousTimeEvolver<ValidatedTag> parameterised_evolver(
            parameterised, TaylorPicardIntegrator(step_maximum_error = 1e-7));
        ARIADNE_TESTING_EXPECT_EQUAL(parameterised_evolver.system().parameter_space().dimension(), 1u);
        parameterised_evolver.configuration().set_maximum_step_size(0.03125);
        const auto parameterised_lower = parameterised_evolver.orbit(
            RealVariablesBox({x == 0}), Real(0.125_x), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(parameterised_lower.final().size() == 1u);
        const auto& parameterised_final = parameterised_lower.final()[0u];
        const auto physical = parameterised_final.constant_parameter_function(1u).get(0u);
        ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0, 1), cast_exact_interval(physical.range())));
        // x'=p-x from x(0)=0 has x(H)=p*(1-exp(-H)) for the same p.
        const auto gain = 1 - exp(FloatDPBounds(-0.125_x, dp));
        const auto residual = parameterised_final.state_function()[0u] - physical * gain;
        ARIADNE_TESTING_NOTE(residual.range());
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(residual.range()),
                                     ExactIntervalType(-0x1p-16_x, 0x1p-16_x)));
        const OrdinaryDifferentialAlgebraicEquation parameterised_auxiliary(
            {dot(x) = p - x}, RealVariablesBox({0 <= p <= 1}), std::nullopt, {let(output) = p + x});
        ContinuousTimeEvolver<ValidatedTag> parameterised_output(parameterised_auxiliary, integrator);
        ARIADNE_TESTING_EXPECT_EQUAL(parameterised_output.system().auxiliary_space(), RealSpace({output}));
        PreconditionedGradedTaylorSeriesIntegrator preconditioned(step_maximum_error = 1e-5);
        for (const auto& equation: {timed, parameterised, parameterised_auxiliary}) {
            ContinuousTimeEvolver<ValidatedTag> consumer(equation, preconditioned);
            consumer.configuration().set_maximum_step_size(0.03125);
            for (auto semantics: {Semantics::UPPER, Semantics::LOWER}) {
                const auto orbit = consumer.orbit(RealVariablesBox({x == 0}), Real(0.125_x), semantics);
                ARIADNE_TESTING_REQUIRE(orbit.final().size() == 1u);
                const auto& final = orbit.final()[0u];
                ARIADNE_TESTING_EXPECT_EQUAL(final.state_dimension(), 1u);
                ARIADNE_TESTING_EXPECT_EQUAL(final.auxiliary_space(), equation.auxiliary_space());
                const FloatDPBounds horizon(0.125_x, dp);
                const auto preconditioned_gain = 1 - exp(-horizon);
                if (equation.parameter_space().dimension() != 0u) {
                    const auto physical_parameter = final.constant_parameter_function(1u).get(0u);
                    const auto analytic_residual = final.state_function()[0u] - physical_parameter * preconditioned_gain;
                    ARIADNE_TESTING_PROVE(subset(cast_exact_interval(analytic_residual.range()), ExactIntervalType(-0x1p-16_x, 0x1p-16_x)));
                } else {
                    const auto analytic_residual = final.state_function()[0u] - (horizon - preconditioned_gain);
                    ARIADNE_TESTING_PROVE(subset(cast_exact_interval(analytic_residual.range()), ExactIntervalType(-0x1p-16_x, 0x1p-16_x)));
                }
            }
        }
        ContinuousTimeEvolver<ApproximateTag> approximate_timed_continuous_evolver(timed);
        ARIADNE_TESTING_PROVE(approximate_timed_continuous_evolver.configuration().step_size().raw().get_d() > 0);
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
        }

        // The selected STRUCTURED integrator now supports timed
        // inclusion LOWER: same physical field and all admissible inputs.
        // x'=t+u, x(0)=0, u in [-1,1], H=1/16 has the exact image
        // [-H+H^2/2, H+H^2/2]=[-31/512,33/512].
        ContinuousTimeEvolver<ValidatedTag> timed_evolver(
            timed, integrator, ZeroApproximation(), ConstantApproximation());
        timed_evolver.configuration().set_maximum_step_size(0.0625);
        timed_evolver.configuration().set_enable_reconditioning(false);
        const auto timed_lower = timed_evolver.orbit(
            RealVariablesBox({x == 0}), Real(0.0625_x), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(timed_lower.final().size() == 1u);
        ARIADNE_TESTING_REQUIRE(timed_lower.reach().size() > 0u);
        const auto& timed_final = timed_lower.final()[0u];
        ARIADNE_TESTING_EXPECT_EQUAL(timed_final.state_dimension(), 1u);
        ARIADNE_TESTING_PROVE(timed_final.number_of_parameters() >= 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(timed_final.number_of_constraints(), 0u);
        const auto timed_box = cast_exact_box(timed_final.euclidean_set().bounding_box());
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {-0.060546875_x, 0.064453125_x}
        }),
                                     timed_box));
        const auto timed_zero = timed_evolver.orbit(
            RealVariablesBox({x == 0}), Real(0), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(timed_zero.final().size(), 1u);

        // A constant physical p remains correlated with the state while
        // independent complete input boxes are added at each actual step.
        ContinuousTimeEvolver<ValidatedTag> parameterised_evolver(
            parameterised, integrator, ZeroApproximation(), ConstantApproximation());
        parameterised_evolver.configuration().set_maximum_step_size(0.03125);
        const auto parameterised_lower = parameterised_evolver.orbit(
            RealVariablesBox({x == 0}), Real(0.0625_x), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(parameterised_lower.final().size() == 1u);
        const auto& parameterised_final = parameterised_lower.final()[0u];
        ARIADNE_TESTING_EXPECT_EQUAL(parameterised_final.constant_parameter_count(), 1u);
        // The one retained constraint is the declared physical p-domain,
        // rather than an artificial equality on a free input coordinate.
        ARIADNE_TESTING_EXPECT_EQUAL(parameterised_final.number_of_constraints(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(parameterised_final.constraint_bounds(), ExactBoxType({
                                                                                  {0, 1}
        }));
        const auto parameter_constraint = parameterised_final.constraint_function()[0u] -
                                          parameterised_final.constant_parameter_function(1u).get(0u);
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(parameter_constraint.range()),
                                     ExactIntervalType(-0x1p-40_x, 0x1p-40_x)));
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {-0.0625_x, 0.125_x}
        }),
                                     cast_exact_box(parameterised_final.state_function().range())));
        ARIADNE_TESTING_PROVE(subset(cast_exact_box(parameterised_final.state_function().range()),
                                     ExactBoxType({
                                         {-0.0703125_x, 0.1328125_x}
        })));
        const auto parameterised_zero = parameterised_evolver.orbit(
            RealVariablesBox({x == 0}), Real(0), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(parameterised_zero.final().size() == 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(parameterised_zero.final()[0u].constant_parameter_count(), 1u);
        ContinuousTimeEvolver<ValidatedTag> boxes(parameterised, FirstOrderMultifunctionIntegrator());
        boxes.configuration().set_maximum_step_size(0.03125);
        const auto boxed_lower = boxes.orbit(
            RealVariablesBox({x == 0}), Real(0.0625_x), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(boxed_lower.final().size() == 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(boxed_lower.reach().size(), 2u);
        const auto& boxed_final = boxed_lower.final()[0u];
        ARIADNE_TESTING_EXPECT_EQUAL(boxed_final.constant_parameter_space(), RealSpace({p}));
        ARIADNE_TESTING_EXPECT_EQUAL(boxed_final.constant_parameter_count(), 1u);
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {0.0_x, 0.0_x}
        }),
                                     cast_exact_box(boxed_final.bounding_box().euclidean_set())));
        const auto boxed_zero = boxes.orbit(
            RealVariablesBox({x == 0}), Real(0), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(boxed_zero.final().size() == 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(boxed_zero.final()[0u].constant_parameter_count(), 1u);
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
