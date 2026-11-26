/***************************************************************************
 *            test_global_time.cpp
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

/** @file tests/test_global_time.cpp
 *  @brief Verify global time behavior.
 */

#include "../src/continuous_integration.hpp"
#include "evolution/continuous_time_approximate_evolver.hpp"
#include "evolution/continuous_time_evolver.hpp"
#include "evolution/hybrid_continuous_time_approximate_evolver.hpp"
#include "evolution/hybrid_continuous_time_evolver.hpp"
#include "evolution/hybrid_evolution_tree.hpp"
#include "integration/integrator.hpp"
#include "integration/selected_inclusion_integrator.hpp"
#include "logging/test_runner.hpp"
#include "system/continuous_time_hybrid_automaton.hpp"
#include "system/odae.hpp"
#include "system/odai.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;
using namespace Ariadne::EvolutionDetail;

class TestGlobalTime : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    const RealVariable x{"x"}, a{"a"}, u{"u"};
    const TimeVariable t;

    void test_clock_conditioning_retains_correlated_arrival_and_auxiliary() {
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-7);
        const auto factory = integrator.function_factory();
        Enclosure initial(ExactBoxType({
                              {-0.0625_x, 0.0625_x}
        }),
                          EnclosureConfiguration(factory));
        initial.set_time_function(2 + factory.create_coordinate(initial.parameter_domain(), 0u));
        initial.set_auxiliary_mapping(make_function(RealSpace({x, t}), Vector<RealExpression>{x + t}),
                                      EnclosureFunctionArguments::STATE_TIME);
        const auto numerical = compute_clocked_equation_step(
            make_function(RealSpace({x, t}), Vector<RealExpression>{t}), initial,
            factory.create_zeros(0u, initial.parameter_domain()), Dyadic(0.125_x), integrator);
        const auto applied = apply_integration_step(initial, numerical.step, numerical.arguments, Semantics::UPPER);
        const auto h = Real(numerical.step.time_step());
        const auto expected = initial.state_function()[0u] + initial.time_function() * h + h * h / 2;
        const auto residual = applied.evolved.state_function()[0u] - expected;
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(residual.range()), ExactIntervalType(-0x1p-30_x, 0x1p-30_x)));
        ARIADNE_TESTING_EXPECT_EQUAL(applied.evolved.state_dimension(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(applied.evolved.number_of_parameters(), initial.number_of_parameters());
        ARIADNE_TESTING_EXPECT_EQUAL(applied.evolved.variable_kinds(), initial.variable_kinds());
        auto check_auxiliary = [&](const Enclosure& enclosure) {
            const auto output_residual = enclosure.auxiliary_function()[0u] - enclosure.state_function()[0u] - enclosure.time_function();
            ARIADNE_TESTING_PROVE(subset(cast_exact_interval(output_residual.range()), ExactIntervalType(-0x1p-28_x, 0x1p-28_x)));
            ARIADNE_TESTING_EXPECT_EQUAL(enclosure.auxiliary_arguments(), EnclosureFunctionArguments::STATE_TIME);
        };
        check_auxiliary(applied.evolved);
        auto conditioned = applied.evolved;
        conditioned.recondition();
        check_auxiliary(conditioned);
        const auto pieces = conditioned.split();
        check_auxiliary(pieces.first);
        check_auxiliary(pieces.second);
        ContinuousTimeEvolver<ValidatedTag> ordinary(OrdinaryDifferentialAlgebraicEquation({dot(x) = t}), integrator);
        ordinary.configuration().set_maximum_step_size(0.03125);
        auto labelled = ordinary.enclosure(ExactBoxType({
            {-0.0625_x, 0.0625_x}
        }));
        labelled.set_time_function(2 + factory.create_coordinate(labelled.parameter_domain(), 0u));
        const auto orbit = ordinary.orbit(labelled, Real(0.125_x), Semantics::UPPER);
        ARIADNE_TESTING_REQUIRE(orbit.final().size() > 0u);
        for (const auto& final: orbit.final()) {
            const auto exact = labelled.state_function()[0u] + labelled.time_function() * Real(0.125_x) + Real(0.0078125_x);
            ARIADNE_TESTING_PROVE(subset(cast_exact_interval((final.state_function()[0u] - exact).range()), ExactIntervalType(-0x1p-28_x, 0x1p-28_x)));
        }
    }

    void test_nonautonomous_equations_and_auxiliaries_in_both_evolvers() {
        const OrdinaryDifferentialAlgebraicEquation equation({dot(x) = t}, {let(a) = x + t});
        ARIADNE_TESTING_PROVE(equation.has_explicit_time());
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-7);
        ContinuousTimeEvolver<ValidatedTag> continuous(equation, integrator);
        continuous.configuration().set_maximum_step_size(0.03125);
        const auto ordinary = continuous.orbit(RealVariablesBox({x == 0}), Real(0.125_x), Semantics::UPPER);
        ARIADNE_TESTING_REQUIRE(ordinary.final().size() > 0u);
        for (const auto& final: ordinary.final()) {
            ARIADNE_TESTING_EXPECT_EQUAL(final.state_dimension(), 1u);
            ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0.0078125_x, 0.0078125_x), cast_exact_interval(final.state_function()[0u].range())));
            ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0.1328125_x, 0.1328125_x), cast_exact_interval(final.auxiliary_function()[0u].range())));
        }
        const DiscreteLocation q(1);
        ContinuousTimeHybridAutomaton system;
        system.new_mode(q, equation);
        for (bool tree: {false, true}) {
            HybridContinuousTimeEvolver<ValidatedTag> hybrid(system);
            hybrid.set_integrator(integrator);
            hybrid.configuration().set_maximum_step_size(0.03125);
            if (tree)
                hybrid.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
            const auto initial = hybrid.enclosure(HybridExactBox(q, RealSpace({
                                                                        x
            }),
                                                                 ExactBoxType({{0, 0}})));
            const auto orbit = hybrid.orbit(initial, HybridTime(0.125_x, 1), Semantics::UPPER);
            ARIADNE_TESTING_REQUIRE(orbit.final().size() > 0u);
            for (const auto& final: orbit.final()) {
                const auto box = cast_exact_box(final.state_bounding_box());
                ARIADNE_TESTING_EXPECT_EQUAL(final.state_space(), RealSpace({x}));
                ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                                 {0.0078125_x, 0.0078125_x}
                }),
                                             box));
                ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                                      {0.00781154632568359375_x, 0.00781345367431640625_x}
                })));
                ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0.1328125_x, 0.1328125_x), cast_exact_interval(final.auxiliary_function()[0u].range())));
            }
            // A deterministic, event-free timed equation is admissible under
            // LOWER: its validated numerical enclosure retains a real
            // trajectory witness, without requiring an INNER image.
            const auto lower_zero = hybrid.orbit(initial, HybridTime(0u, 1), Semantics::LOWER);
            ARIADNE_TESTING_REQUIRE(lower_zero.final().size() == 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(lower_zero.final()[0u].location(), q);
            ARIADNE_TESTING_EXPECT_EQUAL(lower_zero.final()[0u].previous_events().size(), 0u);
            ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0, 0),
                                         cast_exact_interval(lower_zero.final()[0u].state_function()[0u].range())));
            ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0, 0),
                                         cast_exact_interval(lower_zero.final()[0u].auxiliary_function()[0u].range())));

            const auto lower = hybrid.orbit(initial, HybridTime(0.125_x, 1), Semantics::LOWER);
            ARIADNE_TESTING_REQUIRE(lower.final().size() == 1u);
            ARIADNE_TESTING_PROVE(lower.reach().size() > 0u);
            const auto& witnessed = lower.final()[0u];
            ARIADNE_TESTING_EXPECT_EQUAL(witnessed.location(), q);
            ARIADNE_TESTING_EXPECT_EQUAL(witnessed.previous_events().size(), 0u);
            // x' = t, x(0) = 0: x(1/8) = 1/128 and a=x+t=17/128.
            ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0.0078125_x, 0.0078125_x),
                                         cast_exact_interval(witnessed.state_function()[0u].range())));
            ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0.1328125_x, 0.1328125_x),
                                         cast_exact_interval(witnessed.auxiliary_function()[0u].range())));
            ARIADNE_TESTING_PROVE(subset(
                cast_exact_interval(witnessed.time_range()), ExactIntervalType(0.125_x, 0.125_x)));
        }
    }

    void test_nonautonomous_inclusion_with_uncertain_arrival() {
        const DiscreteLocation q(1);
        const OrdinaryDifferentialAlgebraicInclusion inclusion({dot(x) = t + u}, RealVariablesBox({-1 <= u <= 1}));
        ContinuousTimeHybridAutomaton system;
        system.new_mode(q, inclusion);
        for (bool boxes: {false, true}) {
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            TaylorPicardIntegrator base(step_maximum_error = 1e-7);
            if (boxes)
                evolver.set_integrator(FirstOrderMultifunctionIntegrator());
            else
                evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
            evolver.configuration().set_maximum_step_size(0.03125);
            auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                                   x
            }),
                                                            ExactBoxType({{0, 0}})));
            auto source = initial.continuous_set();
            source.adjoin_parameters(ExactBoxType({
                                         {1.9375_x, 2.0_x}
            }),
                                     EnclosureVariableKind::INITIAL);
            source.set_time_function(evolver.function_factory().create_coordinate(source.parameter_domain(), 1u));
            initial.replace_continuous_set(source);
            const auto orbit = evolver.orbit(initial, HybridTime(2.125_x, 1), Semantics::UPPER);
            ARIADNE_TESTING_REQUIRE(orbit.final().size() > 0u);
            std::optional<ExactBoxType> states;
            for (const auto& final: orbit.final()) {
                const auto box = cast_exact_box(final.state_bounding_box());
                states = states ? hull(*states, box) : box;
                ARIADNE_TESTING_EXPECT_EQUAL(final.state_dimension(), 1u);
                ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                                      {0, 1}
                })));
                ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()), ExactIntervalType(2.125_x, 2.125_x)));
            }
            ARIADNE_TESTING_REQUIRE(states.has_value());
            // x(T)=1/2(T^2-tau^2)+integral(u). Both extremes use constant
            // admissible inputs and the endpoints of the arrival interval.
            ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                             {0.1328125_x, 0.568359375_x}
            }),
                                         *states));
        }
    }

    void test_autonomous_modes_with_timed_guards_invariants_and_resets() {
        const DiscreteLocation source(1), target(2);
        const DiscreteEvent jump("timer"), bound("bound");
        ContinuousTimeHybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
        system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, {let(a) = x + t}));
        system.new_invariant(source, t <= Real(0.09375_x), bound);
        system.new_transition(source, jump, target, {next(x) = 1 + t}, t >= Real(0.0625_x), EventKind::URGENT);
        system.check_mode(source);
        for (bool tree: {false, true}) {
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            evolver.configuration().set_maximum_step_size(0.03125);
            if (tree)
                evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
            const auto orbit = evolver.orbit(HybridExactBox(source, RealSpace({
                                                                        x
            }),
                                                            ExactBoxType({{0, 0}})),
                                             HybridTime(0.125_x, 2), Semantics::UPPER);
            bool reset_found = false;
            for (const auto& final: orbit.final()) {
                if (final.location() != target)
                    continue;
                reset_found = true;
                ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().back(), jump);
                const auto box = cast_exact_box(final.state_bounding_box());
                ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                                 {1.0625_x, 1.0625_x}
                }),
                                             box));
                // The guard and reset share the source clock. Prove the
                // tight reset value on the constrained set, not its raw box.
                auto above = final, below = final;
                above.new_activation(jump, make_function(RealSpace({x}), x - Real(1.06250095367431640625_x)));
                below.new_activation(jump, make_function(RealSpace({x}), Real(1.06249904632568359375_x) - x));
                ARIADNE_TESTING_PROVE(definitely(above.is_empty()));
                ARIADNE_TESTING_PROVE(definitely(below.is_empty()));
                ARIADNE_TESTING_PROVE(subset(ExactIntervalType(1.1875_x, 1.1875_x), cast_exact_interval(final.auxiliary_function()[0u].range())));
                ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()), ExactIntervalType(0.125_x, 0.125_x)));
            }
            ARIADNE_TESTING_PROVE(reset_found);
        }
        HybridContinuousTimeEvolver<ApproximateTag> approximate(system);
        approximate.configuration().set_step_size(0.03125);
        const HybridContinuousTimeEvolver<ApproximateTag>::HybridApproximatePointType point(source, RealSpace({x}), FloatDPApproximationPoint({0.0_x}, dp));
        auto single = approximate.orbit(point, HybridTime(0.125_x, 2));
        ARIADNE_TESTING_REQUIRE(single.curves().size() == 2u);
        const auto& target_path = single.curves()[1u];
        ARIADNE_TESTING_EXPECT_WITHIN(target_path.euclidean_set().end()->first.get_d(), 0.125, 1e-14);
        // Preserve the existing sampled strict-positive activation rule:
        // first positive sample is 3/32. Reset reads that global time.
        ARIADNE_TESTING_EXPECT_WITHIN(target_path.euclidean_set().end()->second[target_path.space()[x]].raw().get_d(), 1.09375, 1e-14);
        ARIADNE_TESTING_EXPECT_WITHIN(target_path.euclidean_set().end()->second[target_path.space()[a]].raw().get_d(), 1.21875, 1e-14);
        auto multiple = approximate.orbit(HybridContinuousTimeEvolver<ApproximateTag>::HybridApproximateListPointType{point}, HybridTime(0.125_x, 2));
        ARIADNE_TESTING_REQUIRE(multiple.curves(0u).size() == 2u);
        const auto& list_target = multiple.curves(0u)[1u];
        ARIADNE_TESTING_EXPECT_WITHIN(list_target.euclidean_set().end()->second[list_target.space()[x]].raw().get_d(), 1.09375, 1e-14);
    }

    void test_inclusion_timed_events_with_autonomous_and_nonautonomous_fields() {
        const DiscreteLocation source(1), target(2);
        const DiscreteEvent jump("input timer"), bound("input horizon");
        for (bool timed_field: {false, true}) {
            ContinuousTimeHybridAutomaton system;
            const RealExpression drift = timed_field ? RealExpression(t) : RealExpression(0);
            system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion({dot(x) = drift + u}, RealVariablesBox({-1 <= u <= 1})));
            system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, {let(a) = x + t}));
            system.new_invariant(source, t <= Real(0.09375_x), bound);
            system.new_transition(source, jump, target, {next(x) = x + t}, t >= Real(0.0625_x), EventKind::URGENT);
            for (bool boxes: {false, true}) {
                for (bool tree: {false, true}) {
                    HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                    TaylorPicardIntegrator integrator(step_maximum_error = 1e-7);
                    if (boxes)
                        evolver.set_integrator(FirstOrderMultifunctionIntegrator());
                    else
                        evolver.set_integrator(integrator, ZeroApproximation(), ConstantApproximation());
                    if (tree)
                        evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
                    evolver.configuration().set_maximum_step_size(0.03125);
                    const auto orbit = evolver.orbit(HybridExactBox(source, RealSpace({
                                                                                x
                    }),
                                                                    ExactBoxType({{0, 0}})),
                                                     HybridTime(0.125_x, 2), Semantics::UPPER);
                    std::optional<ExactBoxType> states;
                    for (const auto& final: orbit.final()) {
                        if (final.location() != target)
                            continue;
                        ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().back(), jump);
                        ARIADNE_TESTING_EXPECT_EQUAL(final.state_dimension(), 1u);
                        const auto box = cast_exact_box(final.state_bounding_box());
                        states = states ? hull(*states, box) : box;
                        ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                                              {-1, 1}
                        })));
                        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()), ExactIntervalType(0.125_x, 0.125_x)));
                    }
                    ARIADNE_TESTING_REQUIRE(states.has_value());
                    const ExactIntervalType expected = timed_field ? ExactIntervalType(0.001953125_x, 0.126953125_x) : ExactIntervalType(0, 0.125_x);
                    ARIADNE_TESTING_PROVE(subset(ExactBoxType({expected}), *states));
                }
            }
        }
    }

    void test_approximate_equation_evolution_uses_external_time() {
        const OrdinaryDifferentialAlgebraicEquation equation({dot(x) = t}, {let(a) = x + t});
        ContinuousTimeEvolver<ApproximateTag> evolver(equation);
        evolver.configuration().set_step_size(0.03125);
        const ContinuousTimeEvolver<ApproximateTag>::ApproximatePointType initial(RealSpace({x}), FloatDPApproximationPoint({0.0_x}, dp));
        auto orbit = evolver.orbit(ContinuousTimeEvolver<ApproximateTag>::ApproximateListPointType{initial}, Real(0.1_dec));
        const auto& curve = orbit.curve(0u);
        ARIADNE_TESTING_EXPECT_WITHIN(curve.end()->first.get_d(), 0.1, 1e-15);
        ARIADNE_TESTING_EXPECT_WITHIN(curve.end()->second[curve.space()[x]].raw().get_d(), 0.005, 1e-14);
        ARIADNE_TESTING_EXPECT_WITHIN(curve.end()->second[curve.space()[a]].raw().get_d(), 0.105, 1e-14);
        const DiscreteLocation q(1);
        ContinuousTimeHybridAutomaton system;
        system.new_mode(q, equation);
        HybridContinuousTimeEvolver<ApproximateTag> hybrid(system);
        hybrid.configuration().set_step_size(0.03125);
        const HybridContinuousTimeEvolver<ApproximateTag>::HybridApproximatePointType point(q, RealSpace({x}), FloatDPApproximationPoint({0.0_x}, dp));
        const auto single = hybrid.orbit(point, HybridTime(Real(0.1_dec), 2));
        const auto& final_curve = single.curves()[0u];
        ARIADNE_TESTING_EXPECT_WITHIN(final_curve.euclidean_set().end()->first.get_d(), 0.1, 1e-15);
        ARIADNE_TESTING_EXPECT_WITHIN(final_curve.euclidean_set().end()->second[final_curve.space()[x]].raw().get_d(), 0.005, 1e-14);
        ARIADNE_TESTING_EXPECT_WITHIN(final_curve.euclidean_set().end()->second[final_curve.space()[a]].raw().get_d(), 0.105, 1e-14);
        auto multiple = hybrid.orbit(HybridContinuousTimeEvolver<ApproximateTag>::HybridApproximateListPointType{point, point}, HybridTime(Real(0.125_x), 2));
        ARIADNE_TESTING_EXPECT_EQUAL(multiple.number_of_trajectories(), 2u);
        for (Nat i = 0u; i != 2u; ++i) {
            const auto& path = multiple.curves(i)[0u];
            ARIADNE_TESTING_EXPECT_WITHIN(path.euclidean_set().end()->second[path.space()[x]].raw().get_d(), 0.0078125, 1e-14);
        }
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_global_time");
    ARIADNE_TESTING_ADD_TEST(TestGlobalTime, test_clock_conditioning_retains_correlated_arrival_and_auxiliary);
    ARIADNE_TESTING_ADD_TEST(TestGlobalTime, test_nonautonomous_equations_and_auxiliaries_in_both_evolvers);
    ARIADNE_TESTING_ADD_TEST(TestGlobalTime, test_nonautonomous_inclusion_with_uncertain_arrival);
    ARIADNE_TESTING_ADD_TEST(TestGlobalTime, test_autonomous_modes_with_timed_guards_invariants_and_resets);
    ARIADNE_TESTING_ADD_TEST(TestGlobalTime, test_approximate_equation_evolution_uses_external_time);
    ARIADNE_TESTING_ADD_TEST(TestGlobalTime, test_inclusion_timed_events_with_autonomous_and_nonautonomous_fields);
    return Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
