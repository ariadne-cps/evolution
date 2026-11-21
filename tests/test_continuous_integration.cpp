/***************************************************************************
 *            test_continuous_integration.cpp
 *  Copyright  2026  Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#include "../src/continuous_integration.hpp"
#include "evolution/enclosure.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "function/taylor_model.hpp"
#include "integration/validated_step_integrator.hpp"
#include "logging/test_runner.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;
using namespace Ariadne::EvolutionDetail;

class TestContinuousIntegration : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    const ExactBoxType domain{
        {-1.0_x, 1.0_x}
    };
    const StepSizeType h{0.125_x};
    TaylorPicardIntegrator integrator{StepMaximumError(1e-6)};
    EnclosureConfiguration config{integrator.function_factory()};

    Void test_local_step_preserves_correlations_and_absolute_time() {
        const auto s = EffectiveScalarMultivariateFunction::coordinate(1u, 0u);
        Enclosure source(domain, EffectiveVectorMultivariateFunction({s, -s}),
                         EffectiveScalarMultivariateFunction::constant(1u, 2),
                         {ValidatedConstraint(0, s, 1)}, config);
        const ValidatedVectorMultivariateFunction field = {
            ValidatedScalarMultivariateFunction::constant(2u, 1),
            ValidatedScalarMultivariateFunction::constant(2u, -1)};
        const FlowStepIntegratorAdapter solver(integrator);
        const auto step = solver.step(ValidatedIntegrationRequest(field,
                                                                  cast_exact_box(source.bounding_box()), suggest(h), StepSizeType(2u)));
        auto reach = source;
        auto evolved = source;
        apply_integration_step(reach, evolved, step, IntegrationStepCoordinates::STATE);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.parameter_domain(), domain);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.number_of_constraints(), source.number_of_constraints());
        ARIADNE_TESTING_PROVE(reach.number_of_constraints() >= source.number_of_constraints());
        const auto invariant = evolved.state_function()[0u] + evolved.state_function()[1u];
        // 2^-20 is exactly representable even when long double exceeds double.
        const ExactDouble tolerance = 0x1p-20_x;
        ARIADNE_TESTING_PROVE(definitely(subset(cast_exact_interval(invariant.range()), ExactIntervalType(-tolerance, tolerance))));
        const Vector<FloatDPBounds> point({FloatDPBounds(0.25_x, dp)});
        const auto values = evolved.state_function()(point);
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(StepSizeType(0.25_x) + step.time_step(), dp), values[0u]));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0, dp), values[0u] + values[1u]));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(step.final_time(), dp), evolved.time_function()(point)));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(step.time_step(), dp), evolved.dwell_time_function()(point)));
    }

    Void test_parameterized_step_uses_the_independent_endpoint_certificate() {
        Enclosure source(domain, config);
        source.set_time_function(config.function_factory().create_constant(domain, 2));
        const auto reach_domain = product(domain, IntervalDomainType(StepSizeType(0u), h));
        const auto tube = config.function_factory().create_projection(reach_domain, range(0u, 1u));
        auto endpoint = ValidatedVectorMultivariateTaylorFunctionModelDP::identity(domain, ThresholdSweeper<FloatDP>(dp, 1e-10));
        endpoint[0u].set_error(FloatDPError(0.125_x, dp));
        const ValidatedIntegrationStep step(h, tube.range(), endpoint.range(), {},
                                            {tube}, endpoint, StepSizeType(2u), 1u, false, ExactBoxType(0u), StepSizeType(2u));
        auto reach = source;
        auto evolved = source;
        apply_integration_step(reach, evolved, step, IntegrationStepCoordinates::PARAMETERS);
        const Vector<FloatDPBounds> endpoint_point({FloatDPBounds(0, dp)});
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(-0.125_x, 0.125_x, dp), evolved.state_function()(endpoint_point)[0u]));
        ARIADNE_TESTING_EXPECT_EQUAL(reach.number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.number_of_parameters(), 1u);
        const Vector<FloatDPBounds> point({FloatDPBounds(0, dp)});
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(StepSizeType(2u) + h, dp), evolved.time_function()(point)));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(h, dp), evolved.dwell_time_function()(point)));
    }

    Void test_unsupported_certificates_do_not_change_enclosures() {
        Enclosure reach(domain, config), evolved(domain, config);
        const ValidatedIntegrationStep box_only(h, UpperBoxType(domain), UpperBoxType(domain));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      apply_integration_step(reach, evolved, box_only, IntegrationStepCoordinates::STATE));
        const auto d1 = product(domain, IntervalDomainType(StepSizeType(0u), hlf(h)));
        const auto d2 = product(domain, IntervalDomainType(hlf(h), h));
        const auto& factory = config.function_factory();
        const ValidatedIntegrationStep split(h, UpperBoxType(domain), UpperBoxType(domain), {},
                                             {factory.create_projection(d1, range(0u, 1u)), factory.create_projection(d2, range(0u, 1u))},
                                             factory.create_identity(domain), StepSizeType(0u), 1u);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      apply_integration_step(reach, evolved, split, IntegrationStepCoordinates::PARAMETERS));
        ARIADNE_TESTING_EXPECT_EQUAL(reach.parameter_domain(), domain);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.parameter_domain(), domain);
    }

    Void test_argument_maps_reject_missing_or_out_of_domain_parameters() {
        Enclosure source(domain, config);
        const ExactBoxType parameters({
            {0.5_x, 0.75_x}
        });
        source.adjoin_constant_parameters(parameters);
        const auto t = EffectiveScalarMultivariateFunction::coordinate(3u, 1u);
        const auto p = EffectiveScalarMultivariateFunction::coordinate(3u, 2u);
        const EffectiveVectorMultivariateFunction field = {t + p};
        TaylorSeriesIntegrator taylor(step_maximum_error = 1e-6, order = 8);
        const FlowStepIntegratorAdapter solver(taylor);
        const auto step = solver.step(ValidatedIntegrationRequest(
            field, domain, StepSizeType(0u), h, parameters, UpperBoxType({
                                                                {-2.0_x, 2.0_x}
        })));
        auto reach = source;
        auto evolved = source;
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      apply_integration_step(reach, evolved, step, source.state_function()));
        const auto fixed = config.function_factory().create_coordinate(source.parameter_domain(), 1u);
        const auto outside = join(source.state_function(), fixed + Real(1));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      apply_integration_step(reach, evolved, step, outside));
        ARIADNE_TESTING_EXPECT_EQUAL(reach.parameter_domain(), source.parameter_domain());
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.parameter_domain(), source.parameter_domain());
        ARIADNE_TESTING_EXPECT_EQUAL(cast_exact_box(reach.state_function().range()), cast_exact_box(source.state_function().range()));
        ARIADNE_TESTING_EXPECT_EQUAL(cast_exact_box(evolved.state_function().range()), cast_exact_box(source.state_function().range()));
    }

    Void test_explicit_time_state_step_rebases_reach_without_losing_correlations() {
        // The numerical reach uses absolute time [1/2,5/8], whereas Enclosure
        // must gain a local temporal parameter [0,1/8].
        const auto s = EffectiveScalarMultivariateFunction::coordinate(1u, 0u);
        Enclosure source(domain, EffectiveVectorMultivariateFunction({s, -s}),
                         EffectiveScalarMultivariateFunction::constant(1u, 0.5_x), {}, config);
        const auto t = EffectiveScalarMultivariateFunction::coordinate(3u, 2u);
        const EffectiveVectorMultivariateFunction field = {t, -t};
        const UpperBoxType bound({
            {-2.0_x, 2.0_x},
            {-2.0_x, 2.0_x}
        });
        TaylorSeriesIntegrator taylor(step_maximum_error = 1e-6, order = 8);
        const FlowStepIntegratorAdapter solver(taylor);
        const StepSizeType t0(0.5_x);
        const auto step = solver.step(ValidatedIntegrationRequest(
            field, cast_exact_box(source.bounding_box()), t0, h, ExactBoxType(0u), bound));
        ARIADNE_TESTING_PROVE(step.has_explicit_time());
        ARIADNE_TESTING_EXPECT_EQUAL(step.reach_time_offset(), StepSizeType(0u));
        auto reach = source;
        auto evolved = source;
        apply_integration_step(reach, evolved, step, IntegrationStepCoordinates::STATE);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.parameter_domain(), domain);
        const auto invariant = evolved.state_function()[0u] + evolved.state_function()[1u];
        const ExactDouble tolerance = 0x1p-20_x;
        ARIADNE_TESTING_PROVE(definitely(subset(cast_exact_interval(invariant.range()),
                                                ExactIntervalType(-tolerance, tolerance))));
        const Vector<FloatDPBounds> point({FloatDPBounds(0.25_x, dp)});
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(t0 + h, dp), evolved.time_function()(point)));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(h, dp), evolved.dwell_time_function()(point)));
        ARIADNE_TESTING_EXPECT_EQUAL(reach.number_of_parameters(), source.number_of_parameters() + 1u);
    }

    Void test_explicit_time_parameters_are_preserved_across_two_steps() {
        // x'=t+p; p is constant throughout both steps, not reselected.
        const auto t = EffectiveScalarMultivariateFunction::coordinate(3u, 1u);
        const auto p = EffectiveScalarMultivariateFunction::coordinate(3u, 2u);
        const EffectiveVectorMultivariateFunction field = {t + p};
        const ExactBoxType initial({
            {0.0_x, 0.125_x}
        });
        const ExactBoxType parameters({
            {0.5_x, 0.75_x}
        });
        const ExactBoxType parameter_domain = product(initial, parameters);
        const auto x0 = EffectiveScalarMultivariateFunction::coordinate(2u, 0u);
        Enclosure source(parameter_domain, EffectiveVectorMultivariateFunction({x0}),
                         EffectiveScalarMultivariateFunction::constant(2u, 0.5_x), {}, config);
        const UpperBoxType bound({
            {-1.0_x, 2.0_x}
        });
        const StepSizeType t0(0.5_x);
        const StepSizeType h1(0.0625_x);
        const StepSizeType h2(0.03125_x);
        TaylorSeriesIntegrator taylor(step_maximum_error = 1e-6, order = 8);
        const FlowStepIntegratorAdapter solver(taylor);
        const auto first = solver.step(ValidatedIntegrationRequest(field, initial, t0, h1, parameters, bound));
        const auto second = solver.step(ValidatedIntegrationRequest(field, first, h2, bound));
        // First-step fields use [x0, absolute time, p], while continuation
        // retains [x0, p, absolute time]. Both must keep the same original p.
        ARIADNE_TESTING_EXPECT_EQUAL(first.reach_time_index().value(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(second.reach_time_index().value(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(first.constant_parameter_domain(), parameters);
        ARIADNE_TESTING_EXPECT_EQUAL(second.constant_parameter_domain(), parameters);
        auto first_reach = source;
        auto first_evolved = source;
        apply_integration_step(first_reach, first_evolved, first, IntegrationStepCoordinates::PARAMETERS);
        auto second_reach = first_evolved;
        auto second_evolved = first_evolved;
        apply_integration_step(second_reach, second_evolved, second, IntegrationStepCoordinates::PARAMETERS);
        ARIADNE_TESTING_EXPECT_EQUAL(first_evolved.parameter_domain(), parameter_domain);
        ARIADNE_TESTING_EXPECT_EQUAL(second_evolved.parameter_domain(), parameter_domain);
        ARIADNE_TESTING_EXPECT_EQUAL(second_evolved.number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(second_reach.number_of_parameters(), 3u);
        const auto p0 = EffectiveScalarMultivariateFunction::coordinate(2u, 1u);
        const EffectiveVectorMultivariateFunction expected = {
            x0 + p0 * Real(0.09375_x) + Real(0.05126953125_x)};
        ARIADNE_TESTING_EXPECT_LESS(norm(second_evolved.state_function() - expected), 1e-4_pr);
        const Vector<FloatDPBounds> point({FloatDPBounds(0, dp), FloatDPBounds(0.625_x, dp)});
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(t0 + h1 + h2, dp), second_evolved.time_function()(point)));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(h1 + h2, dp), second_evolved.dwell_time_function()(point)));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_continuous_integration");
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_local_step_preserves_correlations_and_absolute_time);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_parameterized_step_uses_the_independent_endpoint_certificate);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_unsupported_certificates_do_not_change_enclosures);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_argument_maps_reject_missing_or_out_of_domain_parameters);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_explicit_time_state_step_rebases_reach_without_losing_correlations);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_explicit_time_parameters_are_preserved_across_two_steps);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
