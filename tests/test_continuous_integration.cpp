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
        ARIADNE_TESTING_PROVE(definitely(subset(cast_exact_interval(invariant.range()), ExactIntervalType(-0.000001_x, 0.000001_x))));
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
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_continuous_integration");
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_local_step_preserves_correlations_and_absolute_time);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_parameterized_step_uses_the_independent_endpoint_certificate);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_unsupported_certificates_do_not_change_enclosures);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
