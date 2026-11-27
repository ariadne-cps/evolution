/***************************************************************************
 *            test_continuous_integration.cpp
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

/** @file tests/test_continuous_integration.cpp
 *  @brief Verify continuous integration behavior.
 */

#include "../src/continuous_integration.hpp"
#include "evolution/enclosure.hpp"
#include "function/function.hpp"
#include "function/symbolic_function.hpp"
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

    Void test_multipatch_reach_keeps_local_intervals_and_independent_inputs() {
        Enclosure source(domain, config);
        const auto s = EffectiveScalarMultivariateFunction::coordinate(1u, 0u);
        source.new_positive_state_constraint(s);
        source.set_time_function(config.function_factory().create_constant(domain, 2));
        const auto& factory = config.function_factory();
        const auto patch = [&](const ExactIntervalType& time, const ExactIntervalType& input) {
            const auto d = product(product(domain, time), input);
            return ValidatedVectorMultivariateFunctionPatch(1u,
                                                            factory.create_coordinate(d, 0u) + factory.create_coordinate(d, 2u) *
                                                                                                   (factory.create_coordinate(d, 1u) - Real(2)));
        };
        const auto first = patch(ExactIntervalType(2.0_x, 2.0625_x), ExactIntervalType(-0.5_x, 0.5_x));
        const auto second = patch(ExactIntervalType(2.0625_x, 2.125_x), ExactIntervalType(-1.0_x, 1.0_x));
        const auto endpoint_domain = product(domain, ExactIntervalType(-0.125_x, 0.125_x));
        const auto endpoint = ValidatedVectorMultivariateFunctionPatch(1u,
                                                                       factory.create_coordinate(endpoint_domain, 0u) + factory.create_coordinate(endpoint_domain, 1u));
        const UpperBoxType bound({
            {-2.0_x, 2.0_x}
        });
        // Deliberately unsorted, with different reach/endpoint input domains.
        const ValidatedIntegrationStep step(h, bound, bound, {}, {second, first}, endpoint,
                                            StepSizeType(2u), 1u, true);
        const auto applied = apply_integration_step(source, step, IntegrationStepCoordinates::PARAMETERS, Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(applied.reach.size(), 2u);
        ARIADNE_TESTING_PROVE(!applied.used_box_enclosure);
        for (SizeType i = 0u; i != applied.reach.size(); ++i) {
            const auto& reach = applied.reach[i];
            const ExactIntervalType expected_time = i == 0u
                                                      ? ExactIntervalType(0.0625_x, 0.125_x)
                                                      : ExactIntervalType(0.0_x, 0.0625_x);
            ARIADNE_TESTING_EXPECT_EQUAL(reach.parameter_domain()[2u], expected_time);
            ARIADNE_TESTING_EXPECT_EQUAL(reach.variable_kinds()[1u], EnclosureVariableKind::INPUT);
            ARIADNE_TESTING_EXPECT_EQUAL(reach.variable_kinds()[2u], EnclosureVariableKind::TEMPORAL);
            ARIADNE_TESTING_EXPECT_EQUAL(reach.number_of_constraints(), source.number_of_constraints());
            const ExactDouble tau = i == 0u ? 0.09375_x : 0.03125_x;
            const Vector<FloatDPBounds> point({FloatDPBounds(0.25_x, dp), FloatDPBounds(0.5_x, dp), FloatDPBounds(tau, dp)});
            ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.25_x, dp) + FloatDPBounds(0.5_x, dp) * FloatDPBounds(tau, dp),
                                          reach.state_function()(point)[0u]));
            ARIADNE_TESTING_PROVE(refines(FloatDPBounds(2, dp) + FloatDPBounds(tau, dp), reach.time_function()(point)));
            ARIADNE_TESTING_PROVE(refines(FloatDPBounds(tau, dp), reach.dwell_time_function()(point)));
        }
        ARIADNE_TESTING_EXPECT_EQUAL(applied.reach[0u].parameter_domain()[1u], ExactIntervalType(-1.0_x, 1.0_x));
        ARIADNE_TESTING_EXPECT_EQUAL(applied.reach[1u].parameter_domain()[1u], ExactIntervalType(-0.5_x, 0.5_x));
        ARIADNE_TESTING_EXPECT_EQUAL(applied.evolved.parameter_domain(), endpoint_domain);
        const Vector<FloatDPBounds> endpoint_point({FloatDPBounds(0, dp), FloatDPBounds(0.125_x, dp)});
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.125_x, dp), applied.evolved.state_function()(endpoint_point)[0u]));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(2.125_x, dp), applied.evolved.time_function()(endpoint_point)));

        const ValidatedIntegrationStep boxed_endpoint(h, bound, bound, {}, {second, first}, std::nullopt,
                                                      StepSizeType(2u), 1u, true);
        const auto mixed = apply_integration_step(source, boxed_endpoint, IntegrationStepCoordinates::PARAMETERS,
                                                  Semantics::UPPER, IntegrationStepRepresentation::ALLOW_BOX_UPPER_ENCLOSURES);
        ARIADNE_TESTING_PROVE(mixed.used_box_enclosure);
        ARIADNE_TESTING_EXPECT_EQUAL(mixed.reach.size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(mixed.reach[0u].parameter_domain(), applied.reach[0u].parameter_domain());
        ARIADNE_TESTING_EXPECT_EQUAL(mixed.evolved.variable_kinds()[1u], EnclosureVariableKind::ERROR);
        const ValidatedIntegrationStep boxed_reach(h, bound, bound, {}, {}, endpoint, StepSizeType(2u));
        const auto other = apply_integration_step(source, boxed_reach, IntegrationStepCoordinates::PARAMETERS,
                                                  Semantics::UPPER, IntegrationStepRepresentation::ALLOW_BOX_UPPER_ENCLOSURES);
        ARIADNE_TESTING_EXPECT_EQUAL(other.evolved.parameter_domain(), endpoint_domain);
        ARIADNE_TESTING_EXPECT_EQUAL(other.evolved.variable_kinds()[1u], EnclosureVariableKind::INPUT);
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.125_x, dp), other.evolved.state_function()(endpoint_point)[0u]));
    }

    Void test_box_only_inclusion_continuation_uses_outer_enclosures() {
        const auto u = EffectiveScalarMultivariateFunction::coordinate(2u, 1u);
        const ValidatedInputFieldBoxAdapterDP field(
            EffectiveVectorMultivariateFunction({u}), domain, domain);
        const ExactBoxType initial({
            {-0.125_x, 0.125_x}
        });
        Enclosure source(initial, config);
        const FirstOrderMultifunctionIntegrator numerical;
        const MultiflowStepIntegratorAdapter solver(numerical);
        const UpperBoxType bound(domain);
        const StepSizeType h1(0.03125_x), h2(0.015625_x);
        const auto first = solver.step(ValidatedIntegrationRequest(field, initial, h1, bound));
        const auto first_applied = apply_integration_step(source, first, IntegrationStepCoordinates::STATE,
                                                          Semantics::LOWER, IntegrationStepRepresentation::ALLOW_BOX_LOWER_ENCLOSURES);
        const auto second = solver.step(ValidatedIntegrationRequest(field, first, h2, bound));
        const auto applied = apply_integration_step(first_applied.evolved, second, IntegrationStepCoordinates::STATE,
                                                    Semantics::LOWER, IntegrationStepRepresentation::ALLOW_BOX_LOWER_ENCLOSURES);
        ARIADNE_TESTING_PROVE(first_applied.used_box_enclosure && applied.used_box_enclosure);
        ARIADNE_TESTING_EXPECT_EQUAL(applied.evolved.state_dimension(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(applied.evolved.number_of_parameters(), 3u);
        const auto final = cast_exact_box(applied.evolved.bounding_box());
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {-0.171875_x, 0.171875_x}
        }),
                                     final));
        ARIADNE_TESTING_PROVE(subset(final, ExactBoxType({
                                                {-0.1728515625_x, 0.1728515625_x}
        })));
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {0.046875_x, 0.046875_x}
        }),
                                     final));
        const Vector<FloatDPBounds> point({FloatDPBounds(0, dp), FloatDPBounds(0, dp), FloatDPBounds(0, dp)});
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.046875_x, dp), applied.evolved.time_function()(point)));
        ARIADNE_TESTING_EXPECT_EQUAL(first_applied.evolved.variable_kinds()[1u], EnclosureVariableKind::ERROR);
        ARIADNE_TESTING_EXPECT_EQUAL(applied.evolved.variable_kinds()[2u], EnclosureVariableKind::ERROR);
        ARIADNE_TESTING_PROVE(!second.has_evolve_model());
        ARIADNE_TESTING_PROVE(!second.has_reach_models());
        ARIADNE_TESTING_EXPECT_EQUAL(source.number_of_parameters(), 1u);
    }

    Void test_structured_inclusion_continuation_appends_inputs_independently() {
        const auto u = EffectiveFormula::coordinate(1u);
        const EffectiveVectorMultivariateFunction function = EffectiveVectorMultivariateFormulaFunction(
            2u, Vector<EffectiveFormula>{u});
        const ValidatedInputFieldBoxAdapterDP field(function, domain, domain);
        const ExactBoxType initial({
            {-0.125_x, 0.125_x}
        });
        Enclosure source(initial, config);
        TaylorSeriesIntegrator base(ThresholdSweeper<FloatDP>(dp, 1e-10), lipschitz_tolerance = 0.5_x, order = 8);
        const SelectedParametricMultifunctionIntegrator selected(
            base, InputApproximation(ZeroApproximation()), InputApproximation(ConstantApproximation()));
        const StructuredInclusionStepIntegratorAdapter solver(selected);
        const StepSizeType h1(0.03125_x), h2(0.015625_x);
        const UpperBoxType bound(domain);
        const auto first = solver.step(ValidatedIntegrationRequest(field, initial, h1, bound));
        const auto first_applied = apply_integration_step(source, first, IntegrationStepCoordinates::PARAMETERS, Semantics::UPPER);
        const auto second = solver.step(ValidatedIntegrationRequest(field, first, h2, bound));
        const auto applied = apply_integration_step(first_applied.evolved, second, IntegrationStepCoordinates::PARAMETERS, Semantics::UPPER);
        ARIADNE_TESTING_PROVE(!applied.used_box_enclosure);
        ARIADNE_TESTING_EXPECT_EQUAL(applied.evolved.parameter_domain(), second.evolve_model().domain());
        ARIADNE_TESTING_PROVE(applied.evolved.number_of_parameters() > first_applied.evolved.number_of_parameters());
        ARIADNE_TESTING_EXPECT_EQUAL(applied.reach[0u].parameter_domain()[applied.reach[0u].number_of_parameters() - 1u],
                                     ExactIntervalType(0.0_x, 0.015625_x));
        const auto final = cast_exact_box(applied.evolved.bounding_box());
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {-0.171875_x, 0.171875_x}
        }),
                                     final));
        ARIADNE_TESTING_PROVE(subset(final, ExactBoxType({
                                                {-0.1728515625_x, 0.1728515625_x}
        })));
    }

    Void test_box_enclosures_retain_fixed_parameters_and_constraints() {
        Enclosure source(domain, config);
        const auto s = EffectiveScalarMultivariateFunction::coordinate(1u, 0u);
        source.new_positive_state_constraint(s);
        const ExactBoxType parameters({
            {-0.5_x, 0.5_x}
        });
        source.adjoin_parameters(parameters, EnclosureVariableKind::PARAMETER);
        const auto p = config.function_factory().create_coordinate(source.parameter_domain(), 1u);
        const auto arguments = join(source.state_function(), p);
        const UpperBoxType bound({
            {-2.0_x, 2.0_x}
        });
        const ValidatedIntegrationStep step(h, bound, bound, {}, {}, std::nullopt,
                                            StepSizeType(0u), std::nullopt, true, parameters);
        const auto applied = apply_integration_step(source, step, arguments, Semantics::LOWER,
                                                    IntegrationStepRepresentation::ALLOW_BOX_LOWER_ENCLOSURES);
        ARIADNE_TESTING_EXPECT_EQUAL(applied.evolved.variable_kinds()[1u], EnclosureVariableKind::PARAMETER);
        ARIADNE_TESTING_EXPECT_EQUAL(applied.evolved.variable_kinds()[2u], EnclosureVariableKind::ERROR);
        ARIADNE_TESTING_EXPECT_EQUAL(applied.evolved.parameter_domain()[1u], parameters[0u]);
        ARIADNE_TESTING_EXPECT_EQUAL(applied.evolved.number_of_constraints(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(applied.reach[0u].number_of_constraints(), 1u);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      apply_integration_step(source, step, IntegrationStepCoordinates::PARAMETERS, Semantics::UPPER,
                                                             IntegrationStepRepresentation::ALLOW_BOX_UPPER_ENCLOSURES));
    }

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
        const auto applied = apply_integration_step(source, step, IntegrationStepCoordinates::STATE, Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(applied.reach.size(), 1u);
        const auto& reach = applied.reach[0u];
        const auto& evolved = applied.evolved;
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
        const auto applied = apply_integration_step(source, step, IntegrationStepCoordinates::PARAMETERS, Semantics::UPPER);
        const auto& reach = applied.reach[0u];
        const auto& evolved = applied.evolved;
        const Vector<FloatDPBounds> endpoint_point({FloatDPBounds(0, dp)});
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(-0.125_x, 0.125_x, dp), evolved.state_function()(endpoint_point)[0u]));
        ARIADNE_TESTING_EXPECT_EQUAL(reach.number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.number_of_parameters(), 1u);
        const Vector<FloatDPBounds> point({FloatDPBounds(0, dp)});
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(StepSizeType(2u) + h, dp), evolved.time_function()(point)));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(h, dp), evolved.dwell_time_function()(point)));
    }

    Void test_missing_certificates_require_explicit_box_policy() {
        Enclosure source(domain, config);
        const ValidatedIntegrationStep box_only(h, UpperBoxType(domain), UpperBoxType(domain));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      apply_integration_step(source, box_only, IntegrationStepCoordinates::STATE, Semantics::UPPER));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      apply_integration_step(source, box_only, IntegrationStepCoordinates::STATE, Semantics::LOWER,
                                                             IntegrationStepRepresentation::ALLOW_BOX_UPPER_ENCLOSURES));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      apply_integration_step(source, box_only, IntegrationStepCoordinates::STATE, Semantics::UPPER,
                                                             IntegrationStepRepresentation::ALLOW_BOX_LOWER_ENCLOSURES));
        const auto lower = apply_integration_step(source, box_only, IntegrationStepCoordinates::STATE, Semantics::LOWER,
                                                  IntegrationStepRepresentation::ALLOW_BOX_LOWER_ENCLOSURES);
        ARIADNE_TESTING_PROVE(lower.used_box_enclosure);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.evolved.variable_kinds()[1u], EnclosureVariableKind::ERROR);
        ARIADNE_TESTING_EXPECT_EQUAL(source.parameter_domain(), domain);
    }

    Void test_argument_maps_reject_missing_or_out_of_domain_parameters() {
        Enclosure source(domain, config);
        const ExactBoxType parameters({
            {0.5_x, 0.75_x}
        });
        source.adjoin_parameters(parameters, EnclosureVariableKind::PARAMETER);
        const auto t = EffectiveScalarMultivariateFunction::coordinate(3u, 1u);
        const auto p = EffectiveScalarMultivariateFunction::coordinate(3u, 2u);
        const EffectiveVectorMultivariateFunction field = {t + p};
        TaylorSeriesIntegrator taylor(step_maximum_error = 1e-6, order = 8);
        const FlowStepIntegratorAdapter solver(taylor);
        const auto step = solver.step(ValidatedIntegrationRequest(
            field, domain, StepSizeType(0u), h, parameters, UpperBoxType({
                                                                {-2.0_x, 2.0_x}
        })));
        const auto before = cast_exact_box(source.state_function().range());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      apply_integration_step(source, step, source.state_function(), Semantics::UPPER));
        const auto fixed = config.function_factory().create_coordinate(source.parameter_domain(), 1u);
        const auto outside = join(source.state_function(), fixed + Real(1));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      apply_integration_step(source, step, outside, Semantics::UPPER));
        ARIADNE_TESTING_EXPECT_EQUAL(cast_exact_box(source.state_function().range()), before);
        ARIADNE_TESTING_EXPECT_EQUAL(source.parameter_domain(), product(domain, parameters));
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
        const auto applied = apply_integration_step(source, step, IntegrationStepCoordinates::STATE, Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(applied.reach.size(), 1u);
        const auto& reach = applied.reach[0u];
        const auto& evolved = applied.evolved;
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
        const auto first_applied = apply_integration_step(source, first, IntegrationStepCoordinates::PARAMETERS, Semantics::UPPER);
        const auto& first_evolved = first_applied.evolved;
        const auto second_applied = apply_integration_step(first_evolved, second, IntegrationStepCoordinates::PARAMETERS, Semantics::UPPER);
        const auto& second_reach = second_applied.reach[0u];
        const auto& second_evolved = second_applied.evolved;
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
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_multipatch_reach_keeps_local_intervals_and_independent_inputs);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_box_only_inclusion_continuation_uses_outer_enclosures);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_structured_inclusion_continuation_appends_inputs_independently);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_box_enclosures_retain_fixed_parameters_and_constraints);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_local_step_preserves_correlations_and_absolute_time);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_parameterized_step_uses_the_independent_endpoint_certificate);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_missing_certificates_require_explicit_box_policy);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_argument_maps_reject_missing_or_out_of_domain_parameters);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_explicit_time_state_step_rebases_reach_without_losing_correlations);
    ARIADNE_TESTING_ADD_TEST(TestContinuousIntegration, test_explicit_time_parameters_are_preserved_across_two_steps);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
