/***************************************************************************
 *            test_hybrid_inclusions.cpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#include <memory>

#include "evolution/hybrid_continuous_time_evolver.hpp"
#include "evolution/hybrid_evolution_tree.hpp"
#include "function/taylor_model.hpp"
#include "integration/integrator.hpp"
#include "integration/selected_inclusion_integrator.hpp"
#include "logging/test_runner.hpp"
#include "system/hybrid_automaton.hpp"
#include "system/odae.hpp"
#include "system/odai.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

class TestHybridInclusions : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_autonomous_inclusion_with_both_explicit_solver_choices() {
        const DiscreteLocation q(1);
        const RealVariable x("x"), u("u");
        HybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion({dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        for (bool box_method: {false, true}) {
            HybridContinuousTimeEvolver evolver(system);
            TaylorPicardIntegrator base(step_maximum_error = 1e-6);
            if (box_method)
                evolver.set_integrator(FirstOrderMultifunctionIntegrator());
            else
                evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
            evolver.configuration().set_maximum_step_size(0.03125);
            auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                                   x
            }),
                                                            ExactBoxType({{-0.125_x, 0.125_x}})));
            initial.set_time_function(evolver.function_factory().create_constant(initial.parameter_domain(), Real(2)));
            const auto orbit = evolver.orbit(initial, HybridTime(2.125_x, 1), Semantics::UPPER);
            ARIADNE_TESTING_REQUIRE(orbit.final().size() > 0u);
            ARIADNE_TESTING_PROVE(orbit.reach().size() >= 4u);
            for (const auto& final: orbit.final()) {
                const auto box = cast_exact_box(final.state_bounding_box());
                ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                                 {-0.25_x, 0.25_x}
                }),
                                             box));
                ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                                      {-0.2509765625_x, 0.2509765625_x}
                })));
                ARIADNE_TESTING_EXPECT_EQUAL(final.location(), q);
                ARIADNE_TESTING_EXPECT_EQUAL(final.state_dimension(), 1u);
                ARIADNE_TESTING_PROVE(final.previous_events().empty());
                ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()), ExactIntervalType(2.125_x, 2.125_x)));
                ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.dwell_time_range()), ExactIntervalType(0.125_x, 0.125_x)));
            }
        }
    }

    Void test_reconditioning_protects_carried_parameters() {
        const DiscreteLocation q(1);
        const RealVariable x("x"), u("u");
        HybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion({dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        HybridContinuousTimeEvolver evolver(system);
        TaylorPicardIntegrator base(step_maximum_error = 1e-6);
        evolver.set_integrator(base, ZeroApproximation(), AffineApproximation());
        evolver.configuration().set_maximum_step_size(0.03125);
        auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                               x
        }),
                                                        ExactBoxType({{0, 0}})));
        auto source = initial.continuous_set();
        const ExactBoxType passive({
            {-1, 1},
            {-1, 1},
            {-1, 1},
            {-1, 1},
            {-1, 1},
            {-1, 1}
        });
        source.adjoin_parameters(passive, EnclosureVariableKind::PARAMETER);
        initial.replace_continuous_set(source);
        const auto orbit = evolver.orbit(initial, HybridTime(0.5_x, 1), Semantics::UPPER);
        ARIADNE_TESTING_REQUIRE(orbit.final().size() > 0u);
        for (const auto& final: orbit.final()) {
            SizeType protected_axes = 0u;
            for (auto kind: final.continuous_set().variable_kinds())
                if (kind == EnclosureVariableKind::PARAMETER)
                    ++protected_axes;
            ARIADNE_TESTING_EXPECT_EQUAL(protected_axes, 6u);
            // The pressure limit allows state/time error axes. An affine
            // endpoint can append two fresh input axes on the final step.
            const SizeType limit = (source.configuration().reconditioning_num_blocks() + 1u) + 6u + 2u + 2u;
            ARIADNE_TESTING_PROVE(final.number_of_parameters() <= limit);
            const auto box = cast_exact_box(final.state_bounding_box());
            ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                             {-0.5_x, 0.5_x}
            }),
                                         box));
            ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                                  {-0.5009765625_x, 0.5009765625_x}
            })));
        }
    }

    Void test_reset_into_inclusion_with_varying_arrival_time() {
        const DiscreteLocation source(1), target(2);
        const DiscreteEvent enter("enter");
        const RealVariable s("s"), x("x"), y("y"), u("u"), output("output");
        HybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(s) = 1}, {let(output) = 2 * s}));
        system.new_mode(target, OrdinaryDifferentialAlgebraicInclusion(
                                    {dot(x) = u, dot(y) = 0}, RealVariablesBox({-1 <= u <= 1})));
        system.new_transition(source, enter, target, {next(x) = 0, next(y) = 1}, output >= 0.125_x, EventKind::URGENT);
        for (const Dyadic& horizon: {Dyadic(0.125_x), Dyadic(0.5_x)}) {
            const auto radius = horizon - Dyadic(0.03125_x);
            for (bool tree: {false, true}) {
                HybridContinuousTimeEvolver evolver(system);
                TaylorPicardIntegrator base(step_maximum_error = 1e-6);
                evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
                evolver.configuration().set_maximum_step_size(0.03125);
                if (tree)
                    evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
                const auto orbit = evolver.orbit(HybridExactBox(source, RealSpace({
                                                                            s
                }),
                                                                ExactBoxType({{0, 0.03125_x}})),
                                                 HybridTime(Real(horizon), 2), Semantics::UPPER);
                ARIADNE_TESTING_REQUIRE(orbit.final().size() > 0u);
                std::optional<ExactBoxType> union_box;
                for (const auto& final: orbit.final()) {
                    ARIADNE_TESTING_EXPECT_EQUAL(final.location(), target);
                    ARIADNE_TESTING_EXPECT_EQUAL(final.state_dimension(), 2u);
                    ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().size(), 1u);
                    ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().back(), enter);
                    const auto box = cast_exact_box(final.state_bounding_box());
                    union_box = union_box ? hull(*union_box, box) : box;
                    ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({ExactIntervalType(-radius - Dyadic(0.015625_x), radius + Dyadic(0.015625_x)), ExactIntervalType(0.99999904632568359375_x, 1.00000095367431640625_x)})));
                    ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()),
                                                 ExactIntervalType(horizon - Dyadic(0x1p-20_x), horizon + Dyadic(0x1p-20_x))));
                }
                ARIADNE_TESTING_REQUIRE(union_box.has_value());
                // Arrival spans [1/32,1/16]; the largest input duration is horizon-1/32.
                ARIADNE_TESTING_PROVE(subset(ExactBoxType({ExactIntervalType(-radius, radius), ExactIntervalType(1, 1)}), *union_box));
            }
        }
    }

    Void test_interval_taylor_modes_use_the_same_cycle() {
        const DiscreteLocation q(1);
        const RealVariable x("x");
        const ExactBoxType domain({
            {-1, 1}
        });
        const auto field = ValidatedIntervalTaylorFunctionModel<FloatDP>::constant(
            domain, FloatDPUpperInterval(-1, 1, dp), ThresholdSweeper<FloatDP>(dp, 1e-10));
        HybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion(RealSpace({x}),
                                                                  ValidatedVectorIntervalTaylorFunctionModel<FloatDP>({field})));
        for (bool boxes: {false, true}) {
            HybridContinuousTimeEvolver evolver(system);
            TaylorPicardIntegrator base(step_maximum_error = 1e-6);
            if (boxes)
                evolver.set_integrator(FirstOrderMultifunctionIntegrator());
            else
                evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
            evolver.configuration().set_maximum_step_size(0.03125);
            const auto orbit = evolver.orbit(HybridExactBox(q, RealSpace({
                                                                   x
            }),
                                                            ExactBoxType({{0, 0}})),
                                             HybridTime(0.125_x, 1), Semantics::UPPER);
            ARIADNE_TESTING_REQUIRE(orbit.final().size() > 0u);
            for (const auto& final: orbit.final()) {
                const auto box = cast_exact_box(final.state_bounding_box());
                ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                                 {-0.125_x, 0.125_x}
                }),
                                             box));
                ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                                      {-0.1259765625_x, 0.1259765625_x}
                })));
            }
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                          evolver.orbit(HybridExactBox(q, RealSpace({
                                                                              x
            }),
                                                                       ExactBoxType({{2, 2}})),
                                                        HybridTime(0.125_x, 1), Semantics::UPPER));
        }
    }

    Void test_box_horizon_slice_and_copy_ownership() {
        const DiscreteLocation q(1);
        const RealVariable x("x");
        auto make = [&]() {
            HybridAutomaton system;
            system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion({dot(x) = 1}, RealVariablesBox({})));
            HybridContinuousTimeEvolver evolver(system);
            evolver.set_integrator(FirstOrderMultifunctionIntegrator());
            evolver.configuration().set_maximum_step_size(0.03125);
            return evolver;
        };
        auto original = make();
        std::unique_ptr<HybridContinuousTimeEvolver> clone(original.clone());
        original = make();
        const auto orbit = clone->orbit(HybridExactBox(q, RealSpace({
                                                              x
        }),
                                                       ExactBoxType({{0, 0}})),
                                        HybridTime(Real(1 / 10_q), 1), Semantics::UPPER);
        ARIADNE_TESTING_REQUIRE(orbit.final().size() > 0u);
        for (const auto& final: orbit.final()) {
            const auto box = cast_exact_box(final.state_bounding_box());
            const FloatDPBounds analytic(Real(1 / 10_q), dp);
            ARIADNE_TESTING_PROVE(subset(ExactBoxType({ExactIntervalType(analytic.lower_raw(), analytic.upper_raw())}), box));
            ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                                  {0.09375_x, 0.1005859375_x}
            })));
            ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()),
                                         ExactIntervalType(0.0999755859375_x, 0.100006103515625_x)));
        }
    }

    Void test_state_invariants_filter_inclusion_reach_and_endpoints() {
        const DiscreteLocation q(1);
        const RealVariable x("x"), u("u");
        const DiscreteEvent bound("bound");
        HybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion(
                               {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        system.new_invariant(q, x <= 0, bound);
        for (bool box_method: {false, true}) {
            HybridContinuousTimeEvolver evolver(system);
            TaylorPicardIntegrator base(step_maximum_error = 1e-6);
            if (box_method)
                evolver.set_integrator(FirstOrderMultifunctionIntegrator());
            else
                evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
            evolver.configuration().set_maximum_step_size(0.03125);
            auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({x}),
                                                             ExactBoxType({{0, 0}})));
            const auto orbit = evolver.orbit(initial, HybridTime(0.125_x, 1), Semantics::UPPER);
            ARIADNE_TESTING_REQUIRE(orbit.final().size() > 0u);
            ARIADNE_TESTING_REQUIRE(orbit.reach().size() > 0u);
            for (const auto& final: orbit.final()) {
                ARIADNE_TESTING_EXPECT_EQUAL(final.location(), q);
                ARIADNE_TESTING_PROVE(final.previous_events().empty());
                ARIADNE_TESTING_PROVE(final.number_of_constraints() > 0u);
                ARIADNE_TESTING_PROVE(!definitely(final.is_empty()));
                ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()),
                                             ExactIntervalType(0.125_x, 0.125_x)));
            }
            for (const auto& reach: orbit.reach()) {
                ARIADNE_TESTING_EXPECT_EQUAL(reach.location(), q);
                ARIADNE_TESTING_PROVE(reach.number_of_constraints() > 0u);
            }
        }
    }

    Void test_inclusion_invariants_preserve_closed_zero_horizon() {
        const DiscreteLocation q(1);
        const RealVariable x("x"), u("u");
        HybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion(
                               {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        system.new_invariant(q, x <= 0, DiscreteEvent("bound"));
        HybridContinuousTimeEvolver evolver(system);
        evolver.set_integrator(FirstOrderMultifunctionIntegrator());
        const auto crossing = evolver.orbit(HybridExactBox(q, RealSpace({x}),
                                                           ExactBoxType({{-0.125_x, 0.125_x}})),
                                            HybridTime(0u, 1), Semantics::UPPER);
        ARIADNE_TESTING_REQUIRE(crossing.final().size() > 0u);
        for (const auto& final: crossing.final()) {
            ARIADNE_TESTING_PROVE(!definitely(final.is_empty()));
            ARIADNE_TESTING_PROVE(final.number_of_constraints() > 0u);
            ARIADNE_TESTING_PROVE(final.previous_events().empty());
        }
    }

    Void test_unsupported_capabilities_fail_explicitly() {
        const DiscreteLocation q(1);
        const RealVariable x("x"), u("u"), t("t"), p("p");
        const auto initial = HybridExactBox(q, RealSpace({
                                                   x
        }),
                                            ExactBoxType({{0, 0}}));
        HybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion({dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        HybridContinuousTimeEvolver evolver(system);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, evolver.orbit(initial, HybridTime(0u, 1), Semantics::UPPER));
        evolver.set_integrator(FirstOrderMultifunctionIntegrator());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, evolver.orbit(initial, HybridTime(0u, 1), Semantics::LOWER));
        evolver.set_integrator(TaylorPicardIntegrator(step_maximum_error = 1e-6));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, evolver.orbit(initial, HybridTime(0.125_x, 1), Semantics::UPPER));
        // State-only invariants are now supported as necessary outer
        // restrictions. They need neither an entry time nor a Lie derivative.
        system.new_invariant(q, x <= 1, DiscreteEvent("bound"));
        HybridContinuousTimeEvolver guarded(system);
        guarded.set_integrator(FirstOrderMultifunctionIntegrator());
        const auto initially_valid = guarded.orbit(initial, HybridTime(0u, 1), Semantics::UPPER);
        ARIADNE_TESTING_PROVE(initially_valid.final().size() > 0u);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, guarded.orbit(initial, HybridTime(0u, 1), Semantics::LOWER));
        // A state invariant does not establish entry-time, progress or reset
        // capabilities; in particular an ODAI jump remains unsupported.
        const DiscreteLocation target(2);
        HybridAutomaton jumping;
        jumping.new_mode(q, OrdinaryDifferentialAlgebraicInclusion(
                                {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        jumping.new_mode(target, {dot(x) = 0});
        jumping.new_transition(q, DiscreteEvent("jump"), target, {next(x) = x},
                               x >= 0, EventKind::URGENT);
        HybridContinuousTimeEvolver cannot_jump(jumping);
        cannot_jump.set_integrator(FirstOrderMultifunctionIntegrator());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      cannot_jump.orbit(initial, HybridTime(0u, 1), Semantics::UPPER));

        HybridAutomaton timed;
        timed.new_mode(q, OrdinaryDifferentialAlgebraicInclusion({dot(x) = t + p + u}, RealVariablesBox({-1 <= u <= 1}),
                                                                 RealVariablesBox({0 <= p <= 1}), t));
        HybridContinuousTimeEvolver timed_evolver(timed);
        timed_evolver.set_integrator(FirstOrderMultifunctionIntegrator());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, timed_evolver.orbit(initial, HybridTime(0u, 1), Semantics::UPPER));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_inclusions");
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_autonomous_inclusion_with_both_explicit_solver_choices);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_reconditioning_protects_carried_parameters);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_reset_into_inclusion_with_varying_arrival_time);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_interval_taylor_modes_use_the_same_cycle);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_box_horizon_slice_and_copy_ownership);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_state_invariants_filter_inclusion_reach_and_endpoints);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_inclusion_invariants_preserve_closed_zero_horizon);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_unsupported_capabilities_fail_explicitly);
    return Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
