/***************************************************************************
 *            test_constant_parameter_maps.cpp
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

/** @file tests/test_constant_parameter_maps.cpp
 *  @brief Verify physical parameter maps through normalization and continuation.
 */

#include <algorithm>
#include <cstddef>
#include <stdexcept>

#include "logging/test_runner.hpp"
#include "threading/thread_manager.hpp"

#include "ariadne-evolution.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;

class TestConstantParameterMaps : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    const RealVariable x{"x"}, p{"p"}, q{"q"}, a{"a"};
    const TimeVariable t;

    ValidatedFunctionPatchFactory factory() const {
        return TaylorPicardIntegrator(step_maximum_error = 1e-7).function_factory();
    }

    Enclosure source(const ExactIntervalType& physical) const {
        const auto maker = factory();
        Enclosure result(ExactBoxType({
                             {0, 0}
        }),
                         EnclosureConfiguration(maker));
        result.adjoin_constant_parameters(RealSpace({p}), ExactBoxType({physical}));
        result.set_time_function(maker.create_constant(result.parameter_domain(), Real(2)));
        result.set_auxiliary_mapping(make_function(RealSpace({x, t, p}), Vector<RealExpression>{x + t + p}),
                                     EnclosureFunctionArguments::STATE_TIME_PARAMETERS, 1u);
        return result;
    }

    void check_timed_solution(const Enclosure& set, const Real& duration) {
        const auto parameters = set.constant_parameter_function(1u);
        const auto expected = parameters[0u] * duration + 2 * duration + duration * duration / 2;
        const ExactIntervalType tolerance(-0x1p-26_x, 0x1p-26_x);
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval((set.state_function()[0u] - expected).range()), tolerance));
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval((set.auxiliary_function()[0u] - set.state_function()[0u] - set.time_function() - parameters[0u]).range()), tolerance));
        ARIADNE_TESTING_PROVE(definitely(parameters[0u].error() <= 0));
        ARIADNE_TESTING_EXPECT_EQUAL(set.constant_parameter_space(), RealSpace({p}));
        ARIADNE_TESTING_EXPECT_EQUAL(set.constant_parameter_count(), 1u);
    }

    void test_recombination_and_timed_continuation() {
        struct ConcurrencyScope {
            std::size_t previous = ThreadManager::instance().concurrency();
            ~ConcurrencyScope() { ThreadManager::instance().set_concurrency(previous); }
        } restore;
        const List<Enclosure> inputs{source(ExactIntervalType(2, 3)), source(ExactIntervalType(5, 6)),
                                     source(ExactIntervalType(2, 3)), source(ExactIntervalType(5, 6))};
        const OrdinaryDifferentialAlgebraicEquation equation(
            {dot(x) = p + t}, RealVariablesBox({2 <= p <= 6}), t, {let(a) = x + t + p});
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-7);
        ContinuousTimeEvolver<ValidatedTag> evolver(equation, integrator);
        evolver.configuration().set_maximum_step_size(0.03125);
        const auto held = evolver.orbit(LabelledEnclosure(inputs[0u], RealSpace({x}), RealSpace({a})), Real(0), Semantics::UPPER);
        ARIADNE_TESTING_REQUIRE(held.final().size() == 1u);
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(held.final()[0u].constant_parameter_function(1u).get(0u).range()), ExactIntervalType(2, 3)));
        for (auto execution: {EnclosureRecombinationExecution::SEQUENTIAL, EnclosureRecombinationExecution::PARALLEL}) {
            ThreadManager::instance().set_concurrency(std::min<std::size_t>(2u, ThreadManager::instance().maximum_concurrency()));
            const auto batch = ReferenceEnclosureRecombiner(execution).recombine(inputs);
            ARIADNE_TESTING_REQUIRE(batch.outputs.size() == 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(batch.outputs[0u].source_indices.size(), 4u);
            auto merged = batch.outputs[0u].enclosure;
            ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                             {2, 6}
            }),
                                         merged.constant_parameter_domain()));
            ARIADNE_TESTING_PROVE(subset(ExactIntervalType(2, 6), cast_exact_interval(merged.constant_parameter_function(1u).get(0u).range())));
            ARIADNE_TESTING_PROVE(merged.number_of_protected_parameters() >= 2u);
            // Explicit feasible witnesses retain points from both source domains.
            ARIADNE_TESTING_REQUIRE(merged.number_of_parameters() == 3u);
            for (Nat branch = 0u; branch != 2u; ++branch) {
                Vector<FloatDPBounds> witness(3u, FloatDPBounds(0, dp));
                witness[1u] = FloatDPBounds(0.5_x, dp);
                witness[2u] = FloatDPBounds(branch == 0u ? 0.0_x : 3.0_x, dp);
                const ExactDouble expected = branch == 0u ? 2.5_x : 5.5_x;
                const auto physical = merged.constant_parameter_function(1u).get(0u)(witness);
                const FloatDPBounds target(expected, dp);
                ARIADNE_TESTING_PROVE(physical.lower_raw() <= target.lower_raw() && physical.upper_raw() >= target.upper_raw());
                for (const auto& constraint: merged.constraint_models()) {
                    const auto value = constraint.function()(witness);
                    const FloatDPBounds lower(constraint.lower_bound(), dp), upper(constraint.upper_bound(), dp);
                    ARIADNE_TESTING_PROVE(value.lower_raw() >= lower.upper_raw());
                    ARIADNE_TESTING_PROVE(value.upper_raw() <= upper.lower_raw());
                }
            }
            merged.adjoin_parameters(ExactBoxType(10u, ExactIntervalType(-1, 1)), EnclosureVariableKind::INPUT);
            merged.recondition();
            ARIADNE_TESTING_PROVE(merged.number_of_parameters() < batch.outputs[0u].enclosure.number_of_parameters() + 10u);
            auto current = LabelledEnclosure(merged, RealSpace({x}), RealSpace({a}));
            const auto unchanged = evolver.orbit(current, Real(0), Semantics::UPPER);
            ARIADNE_TESTING_REQUIRE(unchanged.final().size() == 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(unchanged.final()[0u].number_of_parameters(), current.number_of_parameters());
            for (Nat step = 1u; step != 3u; ++step) {
                const auto orbit = evolver.orbit(current, Real(0.125_x), Semantics::UPPER);
                ARIADNE_TESTING_REQUIRE(orbit.final().size() == 1u);
                current = orbit.final()[0u];
                check_timed_solution(current, Real(step) / 8);
            }
            const auto pieces = current.split();
            check_timed_solution(pieces.first, Real(0.25_x));
            check_timed_solution(pieces.second, Real(0.25_x));
            const auto identity = recombine(List<Enclosure>{current});
            check_timed_solution(identity, Real(0.25_x));
        }
    }

    void test_correlated_binding_and_coordinate_transports() {
        const auto maker = factory();
        Enclosure initial(ExactBoxType({
                              {-1, 1}
        }),
                          EnclosureConfiguration(maker));
        auto parameters = maker.create_zeros(1u, initial.parameter_domain());
        parameters.set(0u, initial.state_function()[0u] + 3);
        initial.bind_constant_parameters(RealSpace({
                                             p
        }),
                                         ExactBoxType({{2, 4}}), parameters);
        ARIADNE_TESTING_EXPECT_EQUAL(initial.number_of_protected_parameters(), 1u);
        initial.adjoin_parameters(ExactBoxType(10u, ExactIntervalType(-1, 1)), EnclosureVariableKind::INPUT);
        initial.recondition();
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval((initial.constant_parameter_function(1u).get(0u) - initial.state_function()[0u] - 3).range()), ExactIntervalType(-0x1p-26_x, 0x1p-26_x)));
        auto narrowed = initial;
        const auto domains = split(initial.parameter_domain(), 0u);
        narrowed.restrict(domains.first);
        ARIADNE_TESTING_PROVE(narrowed.has_constant_parameter_binding());
        const auto extended = product(narrowed, ExactIntervalType(-2, 2));
        ARIADNE_TESTING_EXPECT_EQUAL(extended.state_dimension(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(extended.constant_parameter_space(), RealSpace({p}));
        const auto extended_box = product(narrowed, ExactBoxType({
                                                        {-2, 2},
                                                        {-3, 3}
        }));
        ARIADNE_TESTING_EXPECT_EQUAL(extended_box.state_dimension(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(extended_box.constant_parameter_count(), 1u);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, product(narrowed, narrowed));

        const OrdinaryDifferentialAlgebraicEquation equation({dot(x) = p}, RealVariablesBox({2 <= p <= 4}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-7);
        ContinuousTimeEvolver<ValidatedTag> evolver(equation, integrator);
        evolver.configuration().set_maximum_step_size(0.03125);
        auto current = LabelledEnclosure(initial, RealSpace({x}));
        for (Nat step = 1u; step != 3u; ++step) {
            const auto orbit = evolver.orbit(current, Real(0.125_x), Semantics::UPPER);
            ARIADNE_TESTING_REQUIRE(orbit.final().size() == 1u);
            current = orbit.final()[0u];
            const auto expected = (1 + Real(step) / 8) * current.constant_parameter_function(1u).get(0u) - 3;
            ARIADNE_TESTING_PROVE(subset(cast_exact_interval((current.state_function()[0u] - expected).range()), ExactIntervalType(-0x1p-26_x, 0x1p-26_x)));
        }
    }

    void test_rounding_error_becomes_persistent_parameter() {
        const auto maker = factory();
        for (bool tiny: {false, true}) {
            Enclosure initial(ExactBoxType({
                                  {0, 0}
            }),
                              EnclosureConfiguration(maker));
            auto values = maker.create_zeros(1u, initial.parameter_domain());
            values.set(0u, maker.create_constant(initial.parameter_domain(), Real(2)) +
                               (tiny ? FloatDPBounds(-0x1p-50_x, 0x1p-50_x, dp) : FloatDPBounds(-0.125_x, 0.125_x, dp)));
            initial.bind_constant_parameters(RealSpace({
                                                 p
            }),
                                             ExactBoxType({{1.875_x, 2.125_x}}), values);
            ARIADNE_TESTING_PROVE(initial.number_of_parameters() > 1u);
            ARIADNE_TESTING_PROVE(definitely(initial.constant_parameter_function(1u).get(0u).error() <= 0));
            const ExactIntervalType expected = tiny ? ExactIntervalType(0x1.ffffffffffffcp0_x, 0x1.0000000000002p1_x) : ExactIntervalType(1.875_x, 2.125_x);
            ARIADNE_TESTING_PROVE(subset(expected, cast_exact_interval(initial.constant_parameter_function(1u).get(0u).range())));
            const OrdinaryDifferentialAlgebraicEquation equation({dot(x) = p}, RealVariablesBox({Real(1.875_x) <= p <= Real(2.125_x)}));
            TaylorPicardIntegrator integrator(step_maximum_error = 1e-7);
            ContinuousTimeEvolver<ValidatedTag> evolver(equation, integrator);
            evolver.configuration().set_maximum_step_size(0.03125);
            auto current = LabelledEnclosure(initial, RealSpace({x}));
            for (Nat step = 1u; step != 3u; ++step) {
                const auto orbit = evolver.orbit(current, Real(0.125_x), Semantics::UPPER);
                ARIADNE_TESTING_REQUIRE(orbit.final().size() == 1u);
                current = orbit.final()[0u];
                const auto residual = current.state_function()[0u] - current.constant_parameter_function(1u).get(0u) * (Real(step) / 8);
                ARIADNE_TESTING_PROVE(subset(cast_exact_interval(residual.range()), ExactIntervalType(-0x1p-26_x, 0x1p-26_x)));
            }
        }
    }

    void test_inclusion_continuation_retains_parameter_binding() {
        const RealVariable u("u");
        const OrdinaryDifferentialAlgebraicInclusion inclusion(
            {dot(x) = p + u}, RealVariablesBox({Real(-0.125_x) <= u <= Real(0.125_x)}), RealVariablesBox({2 <= p <= 3}));
        ContinuousTimeEvolver<ValidatedTag> evolver(inclusion, FirstOrderMultifunctionIntegrator());
        evolver.configuration().set_maximum_step_size(0.03125);
        auto orbit = evolver.orbit(RealVariablesBox({x == 0}), Real(0.125_x), Semantics::UPPER);
        ARIADNE_TESTING_REQUIRE(orbit.final().size() == 1u);
        auto current = orbit.final()[0u];
        const auto unchanged = evolver.orbit(current, Real(0), Semantics::UPPER);
        ARIADNE_TESTING_REQUIRE(unchanged.final().size() == 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(unchanged.final()[0u].number_of_parameters(), current.number_of_parameters());
        orbit = evolver.orbit(current, Real(0.125_x), Semantics::UPPER);
        ARIADNE_TESTING_REQUIRE(orbit.final().size() == 1u);
        current = orbit.final()[0u];
        ARIADNE_TESTING_EXPECT_EQUAL(current.constant_parameter_space(), RealSpace({p}));
        ARIADNE_TESTING_PROVE(subset(ExactIntervalType(2, 3), cast_exact_interval(current.constant_parameter_function(1u).get(0u).range())));
        ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0.46875_x, 0.78125_x), cast_exact_interval(current.state_function()[0u].range())));
    }

    void test_multiple_parameter_maps_and_order() {
        const auto maker = factory();
        Enclosure first(ExactBoxType({
                            {0, 0}
        }),
                        EnclosureConfiguration(maker));
        Enclosure second = first;
        first.adjoin_constant_parameters(RealSpace({
                                             p, q
        }),
                                         ExactBoxType({{2, 3}, {5, 6}}));
        second.adjoin_constant_parameters(RealSpace({
                                              p, q
        }),
                                          ExactBoxType({{3, 4}, {7, 8}}));
        const List<Enclosure> sources{first, second};
        const auto merged = recombine(sources);
        const auto values = merged.constant_parameter_function(2u);
        ARIADNE_TESTING_EXPECT_EQUAL(merged.constant_parameter_space(), RealSpace({p, q}));
        ARIADNE_TESTING_PROVE(subset(ExactIntervalType(2, 4), cast_exact_interval(values[0u].range())));
        ARIADNE_TESTING_PROVE(subset(ExactIntervalType(5, 8), cast_exact_interval(values[1u].range())));
        ARIADNE_TESTING_PROVE(subset(ExactIntervalType(5, 8), cast_exact_interval(merged.constant_parameter_function(1u).get(0u).range())));
        ARIADNE_TESTING_PROVE(merged.number_of_protected_parameters() >= 4u);
        Enclosure reversed(ExactBoxType({
                               {0, 0}
        }),
                           EnclosureConfiguration(maker));
        reversed.adjoin_constant_parameters(RealSpace({
                                                q, p
        }),
                                            ExactBoxType({{5, 6}, {2, 3}}));
        const List<Enclosure> incompatible{first, reversed};
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, recombine(incompatible));
    }

    void test_incompatible_bindings_are_rejected() {
        const auto maker = factory();
        auto initial = source(ExactIntervalType(2, 3));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, initial.adjoin_constant_parameters(RealSpace({
                                                                                                    p
        }),
                                                                                                ExactBoxType({{2, 3}})));
        ARIADNE_TESTING_EXPECT_THROWS(std::logic_error, initial.constant_parameter_function(2u));
        Enclosure unbound(ExactBoxType({
                              {0, 0}
        }),
                          EnclosureConfiguration(maker));
        ARIADNE_TESTING_EXPECT_THROWS(std::logic_error, unbound.constant_parameter_domain());
        const auto map = maker.create_zeros(1u, unbound.parameter_domain());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, unbound.bind_constant_parameters(RealSpace({
                                                                                                  p, q
        }),
                                                                                              ExactBoxType({{2, 3}}), map));
        ARIADNE_TESTING_PROVE(!unbound.has_constant_parameter_binding());
        auto other = unbound;
        other.adjoin_constant_parameters(RealSpace({
                                             q
        }),
                                         ExactBoxType({{2, 3}}));
        const List<Enclosure> mismatched{initial, other};
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, recombine(mismatched));
        const OrdinaryDifferentialAlgebraicEquation different({dot(x) = q}, RealVariablesBox({2 <= q <= 3}));
        const OrdinaryDifferentialAlgebraicEquation narrower({dot(x) = p}, RealVariablesBox({Real(2.25_x) <= p <= Real(2.75_x)}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-7);
        auto parameter_only = initial;
        parameter_only.set_auxiliary_mapping(EffectiveVectorMultivariateFunction(0u, EuclideanDomain(1u)));
        const LabelledEnclosure labelled(parameter_only, RealSpace({x}));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, ContinuousTimeEvolver<ValidatedTag>(different, integrator).orbit(labelled, Real(0), Semantics::UPPER));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, ContinuousTimeEvolver<ValidatedTag>(narrower, integrator).orbit(labelled, Real(0), Semantics::UPPER));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_constant_parameter_maps");
    ARIADNE_TESTING_ADD_TEST(TestConstantParameterMaps, test_recombination_and_timed_continuation);
    ARIADNE_TESTING_ADD_TEST(TestConstantParameterMaps, test_correlated_binding_and_coordinate_transports);
    ARIADNE_TESTING_ADD_TEST(TestConstantParameterMaps, test_rounding_error_becomes_persistent_parameter);
    ARIADNE_TESTING_ADD_TEST(TestConstantParameterMaps, test_inclusion_continuation_retains_parameter_binding);
    ARIADNE_TESTING_ADD_TEST(TestConstantParameterMaps, test_multiple_parameter_maps_and_order);
    ARIADNE_TESTING_ADD_TEST(TestConstantParameterMaps, test_incompatible_bindings_are_rejected);
    return Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
