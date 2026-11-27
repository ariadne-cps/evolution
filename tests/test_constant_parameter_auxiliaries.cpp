/***************************************************************************
 *            test_constant_parameter_auxiliaries.cpp
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

/** @file tests/test_constant_parameter_auxiliaries.cpp
 *  @brief Verify correlated auxiliary outputs with trajectory-constant parameters.
 */

#include <optional>
#include <stdexcept>

#include "logging/test_runner.hpp"
#include "paving/paver.hpp"

#include "ariadne-evolution.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;

class TestConstantParameterAuxiliaries : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    const RealVariable x{"x"}, p{"p"}, r{"r"}, a{"a"}, b{"b"};
    const TimeVariable t;

    void check_auxiliary(const Enclosure& set, const OrdinaryDifferentialAlgebraicEquation& equation, bool timed) {
        const auto parameters = set.constant_parameter_function(2u);
        const auto pv = parameters[equation.parameter_space()[p]];
        const auto rv = parameters[equation.parameter_space()[r]];
        const auto auxiliary = set.auxiliary_function();
        auto expected = set.state_function()[0u] + 2 * pv + 3 * rv;
        if (timed)
            expected = expected + set.time_function();
        const ExactIntervalType tolerance(-0x1p-28_x, 0x1p-28_x);
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval((auxiliary[equation.auxiliary_space()[a]] - expected).range()), tolerance));
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval((auxiliary[equation.auxiliary_space()[b]] - pv + rv).range()), tolerance));
        ARIADNE_TESTING_EXPECT_EQUAL(set.auxiliary_parameter_count(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(set.state_dimension(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(set.auxiliary_arguments(), timed ? EnclosureFunctionArguments::STATE_TIME_PARAMETERS : EnclosureFunctionArguments::STATE_PARAMETERS);
    }

    void test_analytic_parameter_outputs_and_preserved_bindings() {
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-7);
        const Real h(0.125_x);
        for (bool timed: {false, true}) {
            const RealExpression clock = timed ? RealExpression(t) : RealExpression(0);
            const OrdinaryDifferentialAlgebraicEquation equation(
                {dot(x) = p - r + clock}, RealVariablesBox({2 <= p <= 3, 5 <= r <= 6}),
                timed ? std::optional<RealVariable>(t) : std::nullopt,
                {let(a) = x + 2 * p + 3 * r + clock, let(b) = p - r});
            ContinuousTimeEvolver<ValidatedTag> evolver(equation, integrator);
            evolver.configuration().set_maximum_step_size(0.03125);
            for (Semantics semantics: {Semantics::UPPER, Semantics::LOWER})
                for (Nat initial_kind = 0u; initial_kind != 3u; ++initial_kind) {
                    const Real origin(initial_kind == 2u ? 2 : 0);
                    const auto result = [&] {
                        if (initial_kind == 0u)
                            return evolver.orbit(RealVariablesBox({x == 0}), h, semantics);
                        if (initial_kind == 1u)
                            return evolver.orbit(RealExpressionBoundedConstraintSet({0 <= x <= 0}, {}), h, semantics);
                        auto initial = evolver.enclosure(ExactBoxType({
                            {0, 0}
                        }));
                        initial.set_time_function(integrator.function_factory().create_constant(initial.parameter_domain(), origin));
                        return evolver.orbit(initial, h, semantics);
                    }();
                    ARIADNE_TESTING_REQUIRE(result.final().size() > 0u);
                    for (const auto& final: result.final()) {
                        check_auxiliary(final, equation, timed);
                        const auto parameters = final.constant_parameter_function(2u);
                        auto expected = (parameters[equation.parameter_space()[p]] - parameters[equation.parameter_space()[r]]) * h;
                        if (timed)
                            expected = expected + origin * h + h * h / 2;
                        ARIADNE_TESTING_PROVE(subset(cast_exact_interval((final.state_function()[0u] - expected).range()), ExactIntervalType(-0x1p-28_x, 0x1p-28_x)));
                        for (const auto& reach: result.reach())
                            check_auxiliary(reach, equation, timed);
                        auto conditioned = final;
                        conditioned.adjoin_parameters(ExactBoxType(10u, ExactIntervalType(-1, 1)), EnclosureVariableKind::INPUT);
                        conditioned.recondition();
                        ARIADNE_TESTING_PROVE(conditioned.number_of_parameters() < final.number_of_parameters() + 10u);
                        check_auxiliary(conditioned, equation, timed);
                        const auto split = conditioned.split();
                        check_auxiliary(split.first, equation, timed);
                        check_auxiliary(split.second, equation, timed);
                    }
                }
        }
    }

    void test_parameter_used_only_by_auxiliary() {
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-7);
        const OrdinaryDifferentialAlgebraicEquation equation(
            {dot(x) = 0}, RealVariablesBox({2 <= p <= 3}), std::nullopt, {let(a) = p});
        ContinuousTimeEvolver<ValidatedTag> evolver(equation, integrator);
        evolver.configuration().set_maximum_step_size(0.03125);
        for (const Real& horizon: {Real(0), Real(0.125_x)}) {
            const auto result = evolver.orbit(RealVariablesBox({x == 0}), horizon, Semantics::UPPER);
            ARIADNE_TESTING_REQUIRE(result.final().size() > 0u);
            for (auto final: result.final()) {
                final.adjoin_parameters(ExactBoxType(10u, ExactIntervalType(-1, 1)), EnclosureVariableKind::INPUT);
                final.recondition();
                ARIADNE_TESTING_PROVE(subset(ExactIntervalType(2, 3), cast_exact_interval(final.auxiliary_function()[0u].range())));
                ARIADNE_TESTING_PROVE(subset(cast_exact_interval((final.auxiliary_function()[0u] - final.constant_parameter_function(1u)[0u]).range()), ExactIntervalType(-0x1p-28_x, 0x1p-28_x)));
            }
        }
    }

    void test_invalid_bindings_and_unsupported_conversions() {
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-7);
        const EnclosureConfiguration configuration(integrator.function_factory());
        Enclosure set(ExactBoxType({
                          {0, 0}
        }),
                      configuration);
        const auto mapping = make_function(RealSpace({x, p}), Vector<RealExpression>{x + p});
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, set.set_auxiliary_mapping(mapping, EnclosureFunctionArguments::STATE_PARAMETERS));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, set.argument_function(EnclosureFunctionArguments::STATE, 1u));
        ARIADNE_TESTING_EXPECT_THROWS(std::logic_error, set.set_auxiliary_mapping(mapping, EnclosureFunctionArguments::STATE_PARAMETERS, 1u));
        ARIADNE_TESTING_EXPECT_EQUAL(set.auxiliary_parameter_count(), 0u);
        set.adjoin_parameters(ExactBoxType({
                                  {2, 3}
        }),
                              EnclosureVariableKind::PARAMETER);
        set.set_auxiliary_mapping(mapping, EnclosureFunctionArguments::STATE_PARAMETERS, 1u);
        Paver paver(new AffinePaver());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, outer_approximation(set, Grid(1u), 1u, paver));
        const List<Enclosure> sources{set, set};
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, recombine(sources));
        const auto identity = recombine(List<Enclosure>{set});
        ARIADNE_TESTING_EXPECT_EQUAL(identity.auxiliary_parameter_count(), 1u);
        ARIADNE_TESTING_PROVE(subset(ExactIntervalType(2, 3), cast_exact_interval(identity.auxiliary_function()[0u].range())));
        auto mapped = set;
        mapped.apply_map(EffectiveVectorMultivariateFunction::identity(1u));
        ARIADNE_TESTING_EXPECT_EQUAL(mapped.auxiliary_parameter_count(), 1u);
        mapped.apply_map(EffectiveVectorMultivariateFunction::identity(1u), EffectiveVectorMultivariateFunction::identity(1u));
        ARIADNE_TESTING_EXPECT_EQUAL(mapped.auxiliary_parameter_count(), 0u);
        // Replacing the mapping also replaces its parameter binding.
        set.set_auxiliary_mapping(EffectiveVectorMultivariateFunction::identity(1u));
        ARIADNE_TESTING_EXPECT_EQUAL(set.auxiliary_parameter_count(), 0u);
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_constant_parameter_auxiliaries");
    ARIADNE_TESTING_ADD_TEST(TestConstantParameterAuxiliaries, test_analytic_parameter_outputs_and_preserved_bindings);
    ARIADNE_TESTING_ADD_TEST(TestConstantParameterAuxiliaries, test_parameter_used_only_by_auxiliary);
    ARIADNE_TESTING_ADD_TEST(TestConstantParameterAuxiliaries, test_invalid_bindings_and_unsupported_conversions);
    return Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
