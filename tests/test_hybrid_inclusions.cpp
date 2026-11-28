/***************************************************************************
 *            test_hybrid_inclusions.cpp
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

/** @file tests/test_hybrid_inclusions.cpp
 *  @brief Verify hybrid inclusions behavior.
 */

#include <memory>
#include <optional>

#include "evolution/hybrid_continuous_time_evolver.hpp"
#include "evolution/hybrid_evolution_tree.hpp"
#include "function/taylor_model.hpp"
#include "integration/integrator.hpp"
#include "integration/selected_inclusion_integrator.hpp"
#include "logging/test_runner.hpp"
#include "system/continuous_time_hybrid_automaton.hpp"
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
        ContinuousTimeHybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion({dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        for (bool box_method: {false, true}) {
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
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
        ContinuousTimeHybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion({dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
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
        ContinuousTimeHybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(s) = 1}, {let(output) = 2 * s}));
        system.new_mode(target, OrdinaryDifferentialAlgebraicInclusion(
                                    {dot(x) = u, dot(y) = 0}, RealVariablesBox({-1 <= u <= 1})));
        system.new_transition(source, enter, target, {next(x) = 0, next(y) = 1}, output >= 0.125_x, EventKind::URGENT);
        for (const Dyadic& horizon: {Dyadic(0.125_x), Dyadic(0.5_x)}) {
            const auto radius = horizon - Dyadic(0.03125_x);
            for (bool tree: {false, true}) {
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
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
        ContinuousTimeHybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion(RealSpace({x}),
                                                                  ValidatedVectorIntervalTaylorFunctionModel<FloatDP>({field})));
        for (bool boxes: {false, true}) {
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
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
            ContinuousTimeHybridAutomaton system;
            system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion({dot(x) = 1}, RealVariablesBox({})));
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            evolver.set_integrator(FirstOrderMultifunctionIntegrator());
            evolver.configuration().set_maximum_step_size(0.03125);
            return evolver;
        };
        auto original = make();
        std::unique_ptr<HybridContinuousTimeEvolver<ValidatedTag>> clone(original.clone());
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
        ContinuousTimeHybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion(
                               {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        system.new_invariant(q, x <= 0, bound);
        for (bool box_method: {false, true}) {
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            TaylorPicardIntegrator base(step_maximum_error = 1e-6);
            if (box_method)
                evolver.set_integrator(FirstOrderMultifunctionIntegrator());
            else
                evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
            evolver.configuration().set_maximum_step_size(0.03125);
            auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                                   x
            }),
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
        ContinuousTimeHybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion(
                               {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        system.new_invariant(q, x <= 0, DiscreteEvent("bound"));
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        evolver.set_integrator(FirstOrderMultifunctionIntegrator());
        const auto crossing = evolver.orbit(HybridExactBox(q, RealSpace({
                                                                  x
        }),
                                                           ExactBoxType({{-0.125_x, 0.125_x}})),
                                            HybridTime(0u, 1), Semantics::UPPER);
        ARIADNE_TESTING_REQUIRE(crossing.final().size() > 0u);
        for (const auto& final: crossing.final()) {
            ARIADNE_TESTING_PROVE(!definitely(final.is_empty()));
            ARIADNE_TESTING_PROVE(final.number_of_constraints() > 0u);
            ARIADNE_TESTING_PROVE(final.previous_events().empty());
        }
    }

    Void test_upper_inclusion_transitions_use_guarded_reach() {
        const DiscreteLocation source(1), target(2);
        const RealVariable x("x"), u("u");
        for (EventKind kind: {EventKind::URGENT, EventKind::PERMISSIVE}) {
            ContinuousTimeHybridAutomaton system;
            system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                        {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
            system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
            const DiscreteEvent event(kind == EventKind::URGENT ? "urgent" : "permissive");
            system.new_transition(source, event, target, {next(x) = 0}, x >= 0.03125_x, kind);
            for (bool boxes: {false, true}) {
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                TaylorPicardIntegrator base(step_maximum_error = 1e-6);
                if (boxes)
                    evolver.set_integrator(FirstOrderMultifunctionIntegrator());
                else
                    evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
                evolver.configuration().set_maximum_step_size(0.03125);
                const auto orbit = evolver.orbit(HybridExactBox(source, RealSpace({
                                                                            x
                }),
                                                                ExactBoxType({{0, 0}})),
                                                 HybridTime(0.125_x, 1), Semantics::UPPER);
                // Both kinds admit a feasible crossing witness u=1 at t=1/32,
                // but an outer candidate is not a first-entry certificate.
                Bool jumped = false;
                for (const auto& final: orbit.final()) {
                    if (final.location() == target) {
                        jumped = true;
                        ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().size(), 1u);
                        ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().back(), event);
                        ARIADNE_TESTING_PROVE(definitely(final.inside(
                            HybridExactBox(target, RealSpace({
                                                       x
                        }),
                                           ExactBoxType({{-0x1p-10_x, 0x1p-10_x}})))));
                    }
                }
                ARIADNE_TESTING_PROVE(jumped);
                ARIADNE_TESTING_REQUIRE(orbit.reach().size() > 0u);
                // The ordinary continuous branch remains represented even
                // if some of its points may also spuriously jump.
                Bool source_reach = false;
                for (const auto& reach: orbit.reach())
                    source_reach = source_reach || reach.location() == source;
                ARIADNE_TESTING_PROVE(source_reach);
            }
        }
    }

    Void test_upper_inclusion_progress_blocks_flow_without_reset() {
        const DiscreteLocation q(1);
        const RealVariable x("x"), u("u");
        const DiscreteEvent progress("progress");
        ContinuousTimeHybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion(
                               {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        system.new_guard(q, progress, x >= 0.03125_x, EventKind::PROGRESS);
        for (bool boxes: {false, true}) {
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            TaylorPicardIntegrator base(step_maximum_error = 1e-6);
            if (boxes)
                evolver.set_integrator(FirstOrderMultifunctionIntegrator());
            else
                evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
            evolver.configuration().set_maximum_step_size(0.03125);
            const auto initial = evolver.enclosure(
                HybridExactBox(q, RealSpace({
                                      x
            }),
                               ExactBoxType({{0, 0}})));
            const auto orbit = evolver.orbit(initial, HybridTime(0.125_x, 1), Semantics::UPPER);
            ARIADNE_TESTING_REQUIRE(orbit.reach().size() > 0u);
            for (const auto& reach: orbit.reach()) {
                ARIADNE_TESTING_EXPECT_EQUAL(reach.location(), q);
                ARIADNE_TESTING_PROVE(reach.previous_events().empty());
                ARIADNE_TESTING_PROVE(reach.number_of_constraints() > 0u);
            }
            for (const auto& final: orbit.final()) {
                ARIADNE_TESTING_EXPECT_EQUAL(final.location(), q);
                ARIADNE_TESTING_PROVE(final.previous_events().empty());
                ARIADNE_TESTING_PROVE(final.number_of_constraints() > 0u);
            }
            const auto at_zero = evolver.orbit(initial, HybridTime(0u, 1), Semantics::UPPER);
            ARIADNE_TESTING_REQUIRE(at_zero.final().size() > 0u);
            for (const auto& final: at_zero.final()) {
                ARIADNE_TESTING_PROVE(final.previous_events().empty());
                ARIADNE_TESTING_PROVE(!definitely(final.is_empty()));
            }
        }
    }

    Void test_upper_inclusion_competing_urgent_guards_keep_both_boundary_jumps() {
        const DiscreteLocation source(1), left(2), right(3);
        const RealVariable x("x"), u("u");
        const DiscreteEvent first("first"), second("second");
        ContinuousTimeHybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                    {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        system.new_mode(left, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
        system.new_mode(right, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
        system.new_transition(source, first, left, {next(x) = 0}, x >= 0.03125_x, EventKind::URGENT);
        system.new_transition(source, second, right, {next(x) = 0}, x >= 0.03125_x, EventKind::URGENT);
        for (bool boxes: {false, true}) {
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            TaylorPicardIntegrator base(step_maximum_error = 1e-6);
            if (boxes)
                evolver.set_integrator(FirstOrderMultifunctionIntegrator());
            else
                evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
            evolver.configuration().set_maximum_step_size(0.03125);
            const auto orbit = evolver.orbit(HybridExactBox(source, RealSpace({
                                                                        x
            }),
                                                            ExactBoxType({{0, 0}})),
                                             HybridTime(0.125_x, 1), Semantics::UPPER);
            bool found_first = false, found_second = false;
            for (const auto& final: orbit.final()) {
                if (final.location() == left) {
                    found_first = true;
                    ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().back(), first);
                } else if (final.location() == right) {
                    found_second = true;
                    ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().back(), second);
                }
            }
            ARIADNE_TESTING_PROVE(found_first && found_second);
            ARIADNE_TESTING_REQUIRE(orbit.reach().size() > 0u);
        }
    }

    Void test_upper_inclusion_permissive_guard_does_not_block_flow() {
        const DiscreteLocation source(1), target(2);
        const RealVariable x("x"), u("u");
        const DiscreteEvent event("optional");
        ContinuousTimeHybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                    {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
        system.new_transition(source, event, target, {next(x) = 0},
                              x >= 0.03125_x, EventKind::PERMISSIVE);
        for (bool boxes: {false, true}) {
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            TaylorPicardIntegrator base(step_maximum_error = 1e-6);
            if (boxes)
                evolver.set_integrator(FirstOrderMultifunctionIntegrator());
            else
                evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
            evolver.configuration().set_maximum_step_size(0.03125);
            const auto orbit = evolver.orbit(HybridExactBox(source, RealSpace({
                                                                        x
            }),
                                                            ExactBoxType({{0, 0}})),
                                             HybridTime(0.125_x, 1), Semantics::UPPER);
            bool source_final = false, target_final = false;
            for (const auto& final: orbit.final()) {
                source_final = source_final || final.location() == source;
                target_final = target_final || final.location() == target;
            }
            ARIADNE_TESTING_PROVE(source_final);
            ARIADNE_TESTING_PROVE(target_final);
        }
    }

    Void test_structured_lower_inclusion_preserves_complete_input_flow() {
        const DiscreteLocation q(1);
        const RealVariable x("x"), u("u");
        ContinuousTimeHybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion(
                               {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        TaylorPicardIntegrator base(step_maximum_error = 1e-6);
        evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
        evolver.configuration().set_maximum_step_size(0.03125);
        evolver.configuration().set_enable_reconditioning(false);
        const auto initial = evolver.enclosure(HybridExactBox(
            q, RealSpace({
                   x
        }),
            ExactBoxType({{-0.125_x, 0.125_x}})));
        const auto lower = evolver.orbit(initial, HybridTime(0.125_x, 1), Semantics::LOWER);
        const auto upper = evolver.orbit(initial, HybridTime(0.125_x, 1), Semantics::UPPER);
        ARIADNE_TESTING_REQUIRE(lower.final().size() == 1u);
        ARIADNE_TESTING_REQUIRE(upper.final().size() == 1u);
        ARIADNE_TESTING_PROVE(lower.reach().size() >= 4u);
        const auto& lower_final = lower.final()[0u];
        ARIADNE_TESTING_PROVE(lower_final.location() == q);
        ARIADNE_TESTING_PROVE(lower_final.previous_events().empty());
        ARIADNE_TESTING_EXPECT_EQUAL(lower_final.state_dimension(), 1u);
        const auto lower_box = cast_exact_box(lower_final.state_bounding_box());
        const auto upper_box = cast_exact_box(upper.final()[0u].state_bounding_box());
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {-0.25_x, 0.25_x}
        }),
                                     lower_box));
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {-0.25_x, 0.25_x}
        }),
                                     upper_box));
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(lower_final.time_range()),
                                     ExactIntervalType(0.125_x, 0.125_x)));
        const auto& parameter_kinds = lower_final.continuous_set().variable_kinds();
        const auto& parameters = lower_final.continuous_set().parameter_domain();
        SizeType full_inputs = 0u;
        for (SizeType i = 0u; i != parameter_kinds.size(); ++i) {
            if (parameter_kinds[i] != EnclosureVariableKind::INPUT)
                continue;
            ++full_inputs;
            ARIADNE_TESTING_PROVE(parameters[i].lower_bound().raw().get_d() <= -1.0);
            ARIADNE_TESTING_PROVE(parameters[i].upper_bound().raw().get_d() >= 1.0);
        }
        ARIADNE_TESTING_PROVE(full_inputs >= 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(initial.continuous_set().number_of_parameters(), 1u);

        // A box-only method has not yet been audited for the hybrid LOWER
        // invariant, event and crossing certificates, so it remains rejected
        // independently of ordinary event-free box continuation.
        HybridContinuousTimeEvolver<ValidatedTag> box_evolver(system);
        box_evolver.set_integrator(FirstOrderMultifunctionIntegrator());
        const auto zero = HybridTime(0u, 1);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      box_evolver.orbit(initial, zero, Semantics::LOWER));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      box_evolver.orbit(initial, HybridTime(0.125_x, 1), Semantics::LOWER));
    }

    Void test_lower_inclusion_ignores_optional_jump_but_rejects_blockers() {
        const DiscreteLocation source(1), target(2);
        const RealVariable x("x"), u("u");
        const DiscreteEvent optional("optional");
        ContinuousTimeHybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                    {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
        system.new_transition(source, optional, target, {next(x) = 0},
                              x >= 0.03125_x, EventKind::PERMISSIVE);
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        TaylorPicardIntegrator base(step_maximum_error = 1e-6);
        evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
        evolver.configuration().set_maximum_step_size(0.03125);
        const auto initial = HybridExactBox(
            source, RealSpace({
                        x
        }),
            ExactBoxType({{0, 0}}));

        // Never taking an optional jump is a legal behavior. Even if
        // its guard becomes possible during a step, LOWER can retain
        // uninterrupted flow with the whole input box.
        const auto lower = evolver.orbit(initial, HybridTime(0.125_x, 1), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(lower.final().size() == 1u);
        const auto& final = lower.final()[0u];
        ARIADNE_TESTING_PROVE(final.location() == source);
        ARIADNE_TESTING_PROVE(final.previous_events().empty());
        const auto box = cast_exact_box(final.state_bounding_box());
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {-0.125_x, 0.125_x}
        }),
                                     box));

        // A blocking URGENT transition may invalidate uninterrupted
        // flow and cannot be ignored without an independent event proof.
        ContinuousTimeHybridAutomaton urgent;
        urgent.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                    {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        urgent.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
        urgent.new_transition(source, optional, target, {next(x) = 0},
                              x >= 0, EventKind::URGENT);
        HybridContinuousTimeEvolver<ValidatedTag> blocked(urgent);
        blocked.set_integrator(base, ZeroApproximation(), ConstantApproximation());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      blocked.orbit(initial, HybridTime(0.125_x, 1), Semantics::LOWER));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      blocked.orbit(initial, HybridTime(0u, 1), Semantics::LOWER));

        ContinuousTimeHybridAutomaton constrained;
        constrained.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                         {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        constrained.new_invariant(source, x <= 0, DiscreteEvent("bound"));
        HybridContinuousTimeEvolver<ValidatedTag> invariant(constrained);
        invariant.set_integrator(base, ZeroApproximation(), ConstantApproximation());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      invariant.orbit(initial, HybridTime(0.125_x, 1), Semantics::LOWER));
    }

    Void test_lower_inclusion_certified_whole_tube_invariants() {
        const DiscreteLocation q(1);
        const RealVariable x("x"), u("u");
        const DiscreteEvent bound("bound");
        ContinuousTimeHybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion(
                               {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        system.new_invariant(q, x <= 0, bound);
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        TaylorPicardIntegrator base(step_maximum_error = 1e-6);
        evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
        evolver.configuration().set_maximum_step_size(0.03125);
        evolver.configuration().set_enable_reconditioning(false);

        // The WHOLE x=-1 source and its complete input-dependent reach
        // tube remain strictly inside x<=0 for all times in [0,1/8].
        const auto safe_source = HybridExactBox(
            q, RealSpace({
                   x
        }),
            ExactBoxType({{-1.0_x, -1.0_x}}));
        const auto lower = evolver.orbit(
            safe_source, HybridTime(0.125_x, 1), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(lower.final().size() == 1u);
        ARIADNE_TESTING_PROVE(lower.reach().size() >= 4u);
        const auto& final = lower.final()[0u];
        ARIADNE_TESTING_PROVE(final.location() == q);
        ARIADNE_TESTING_PROVE(final.previous_events().empty());
        const auto bounds = cast_exact_box(final.state_bounding_box());
        ARIADNE_TESTING_PROVE(subset(
            ExactBoxType({
                {-1.125_x, -0.875_x}
        }),
            bounds));
        ARIADNE_TESTING_PROVE(subset(
            cast_exact_interval(final.time_range()),
            ExactIntervalType(0.125_x, 0.125_x)));
        const auto& kinds = final.continuous_set().variable_kinds();
        const auto& parameter_box = final.continuous_set().parameter_domain();
        SizeType free_input_boxes = 0u;
        for (SizeType i = 0u; i < kinds.size(); ++i) {
            if (kinds[i] != EnclosureVariableKind::INPUT)
                continue;
            ++free_input_boxes;
            ARIADNE_TESTING_PROVE(parameter_box[i].lower_bound().raw().get_d() <= -1.0);
            ARIADNE_TESTING_PROVE(parameter_box[i].upper_bound().raw().get_d() >= 1.0);
        }
        ARIADNE_TESTING_PROVE(free_input_boxes >= 4u);

        // Equality at the INITIAL invariant boundary is legal at t=0,
        // but the full input box can subsequently violate x<=0:
        // no fabricated LOWER subbox or input choice is permitted.
        const auto boundary = HybridExactBox(
            q, RealSpace({
                   x
        }),
            ExactBoxType({{0, 0}}));
        const auto zero = evolver.orbit(
            boundary, HybridTime(0u, 1), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(zero.final().size() == 1u);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      evolver.orbit(boundary, HybridTime(0.125_x, 1), Semantics::LOWER));

        // This already-invalid WHOLE source is rejected even at
        // zero horizon: intersecting an outer LOWER enclosure with
        // the invariant does not establish a reachable witness.
        const auto invalid = HybridExactBox(
            q, RealSpace({
                   x
        }),
            ExactBoxType({{0.25_x, 0.25_x}}));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      evolver.orbit(invalid, HybridTime(0u, 1), Semantics::LOWER));
    }

    Void test_lower_inclusion_certainly_inactive_blockers() {
        const DiscreteLocation q(1), target(2);
        const RealVariable x("x"), u("u");
        const auto safe = HybridExactBox(
            q, RealSpace({
                   x
        }),
            ExactBoxType({{-1.0_x, -1.0_x}}));
        const auto near_crossing = HybridExactBox(
            q, RealSpace({
                   x
        }),
            ExactBoxType({{-0.0625_x, -0.0625_x}}));
        const auto active = HybridExactBox(
            q, RealSpace({
                   x
        }),
            ExactBoxType({{0, 0}}));

        // These guards have identical flow-blocking semantics when
        // certainly inactive, despite differing transition/reset rules.
        for (EventKind kind: {
                 EventKind::URGENT, EventKind::PROGRESS, EventKind::IMPACT}) {
            ContinuousTimeHybridAutomaton system;
            system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion(
                                   {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
            const DiscreteEvent event("blocker");
            if (kind == EventKind::PROGRESS) {
                system.new_guard(q, event, x >= 0, kind);
            } else {
                system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
                system.new_transition(q, event, target, {next(x) = x},
                                      x >= 0, kind);
            }
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
            evolver.set_integrator(
                integrator, ZeroApproximation(), ConstantApproximation());
            evolver.configuration().set_maximum_step_size(0.03125);
            evolver.configuration().set_enable_reconditioning(false);

            // All x(t) lie inside [-9/8,-7/8] while all independent
            // numerical inputs range over their COMPLETE [-1,1] box.
            // No blocker can activate, so ordinary flow is legal.
            const auto lower = evolver.orbit(
                safe, HybridTime(0.125_x, 1), Semantics::LOWER);
            ARIADNE_TESTING_REQUIRE(lower.final().size() == 1u);
            ARIADNE_TESTING_PROVE(lower.reach().size() >= 4u);
            const auto& final = lower.final()[0u];
            ARIADNE_TESTING_PROVE(final.location() == q);
            ARIADNE_TESTING_PROVE(final.previous_events().empty());
            ARIADNE_TESTING_EXPECT_EQUAL(
                final.time_range().lower_bound().raw().get_d(), 0.125);
            ARIADNE_TESTING_PROVE(subset(
                ExactBoxType({
                    {-1.125_x, -0.875_x}
            }),
                cast_exact_box(final.state_bounding_box())));
            const auto& kinds = final.continuous_set().variable_kinds();
            const auto& parameters = final.continuous_set().parameter_domain();
            SizeType full_inputs = 0u;
            for (SizeType i = 0u; i < kinds.size(); ++i) {
                if (kinds[i] != EnclosureVariableKind::INPUT)
                    continue;
                ++full_inputs;
                ARIADNE_TESTING_PROVE(
                    parameters[i].lower_bound().raw().get_d() <= -1.0);
                ARIADNE_TESTING_PROVE(
                    parameters[i].upper_bound().raw().get_d() >= 1.0);
            }
            ARIADNE_TESTING_PROVE(full_inputs >= 4u);

            // A source strictly below the guard does NOT prove that
            // the complete following tube stays below it.
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                          evolver.orbit(near_crossing,
                                                        HybridTime(0.125_x, 1), Semantics::LOWER));
            // At an active or undecidable blocker, an unverified flow
            // must not be published, even at zero horizon.
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                          evolver.orbit(active, HybridTime(0u, 1), Semantics::LOWER));
        }
    }

    Void test_lower_inclusion_full_source_urgent_reset_without_flow() {
        const DiscreteLocation source(1), target(2), alternate(3);
        const RealVariable x("x"), u("u");
        const DiscreteEvent first("first"), second("second");
        ContinuousTimeHybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                    {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
        system.new_mode(alternate, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
        system.new_transition(source, first, target, {next(x) = x + 2},
                              x >= 0, EventKind::URGENT);
        system.new_transition(source, second, alternate, {next(x) = x + 4},
                              x >= 0, EventKind::URGENT);
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        TaylorPicardIntegrator base(step_maximum_error = 1e-6);
        evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
        evolver.configuration().set_maximum_step_size(0.03125);
        evolver.configuration().set_enable_reconditioning(false);

        // Both guards hold STRICTLY on the complete original source.
        // Each immediate reset is a real legal behavior, independent
        // of the blocked ordinary progression. Do not order events.
        const auto initial = HybridExactBox(
            source, RealSpace({
                        x
        }),
            ExactBoxType({{1, 1}}));
        const auto orbit = evolver.orbit(
            initial, HybridTime(0.125_x, 2), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(orbit.final().size() == 2u);
        bool found_first = false;
        bool found_second = false;
        for (const auto& final: orbit.final()) {
            ARIADNE_TESTING_REQUIRE(final.previous_events().size() == 1u);
            ARIADNE_TESTING_PROVE(subset(
                cast_exact_interval(final.time_range()),
                ExactIntervalType(0.125_x, 0.125_x)));
            if (final.location() == target) {
                found_first = true;
                ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events()[0u], first);
                ARIADNE_TESTING_PROVE(subset(
                    ExactBoxType({
                        {3, 3}
                }),
                    cast_exact_box(final.state_bounding_box())));
            } else if (final.location() == alternate) {
                found_second = true;
                ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events()[0u], second);
                ARIADNE_TESTING_PROVE(subset(
                    ExactBoxType({
                        {5, 5}
                }),
                    cast_exact_box(final.state_bounding_box())));
            } else {
                ARIADNE_TESTING_REQUIRE(false);
            }
        }
        ARIADNE_TESTING_PROVE(found_first && found_second);

        // An outer guard touching the boundary is UNKNOWN, not a
        // whole-source reset certificate. Do not choose a favourable
        // input or numerical child to make an urgent reset fire.
        const auto boundary = HybridExactBox(
            source, RealSpace({
                        x
        }),
            ExactBoxType({{0, 0}}));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      evolver.orbit(boundary, HybridTime(0.125_x, 2),
                                                    Semantics::LOWER));

        // An IMPACT guard being true is not a direction witness.
        // Reuse the same full-source initial state, but no fabricated
        // impact reset may be emitted without derivative evidence.
        ContinuousTimeHybridAutomaton impact;
        impact.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                    {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        impact.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
        impact.new_transition(source, first, target, {next(x) = x + 2},
                              x >= 0, EventKind::IMPACT);
        HybridContinuousTimeEvolver<ValidatedTag> impact_evolver(impact);
        impact_evolver.set_integrator(
            base, ZeroApproximation(), ConstantApproximation());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      impact_evolver.orbit(initial, HybridTime(0.125_x, 2),
                                                           Semantics::LOWER));
    }

    Void test_lower_inclusion_certified_reset_beside_uncertain_progress() {
        const DiscreteLocation source(1), target(2), alternate(3);
        const RealVariable x("x"), u("u");
        const DiscreteEvent first("first"), second("second"), progress("progress");
        TaylorPicardIntegrator base(step_maximum_error = 1e-6);
        ContinuousTimeHybridAutomaton unblocked;
        unblocked.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                       {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        HybridContinuousTimeEvolver<ValidatedTag> flow(unblocked);
        flow.set_integrator(base, ZeroApproximation(), ConstantApproximation());
        flow.configuration().set_maximum_step_size(0.03125);
        flow.configuration().set_enable_reconditioning(false);
        const auto carried = flow.orbit(
            HybridExactBox(source, RealSpace({
                                       x
        }),
                           ExactBoxType({{0, 0}})),
            HybridTime(0.125_x, 1), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(carried.final().size() == 1u);
        SizeType full_inputs = 0u;
        const auto& carried_set = carried.final()[0u].continuous_set();
        for (SizeType i = 0u; i != carried_set.number_of_parameters(); ++i) {
            if (carried_set.variable_kinds()[i] != EnclosureVariableKind::INPUT)
                continue;
            ++full_inputs;
            ARIADNE_TESTING_PROVE(carried_set.parameter_domain()[i].lower_bound().raw().get_d() <= -1.0);
            ARIADNE_TESTING_PROVE(carried_set.parameter_domain()[i].upper_bound().raw().get_d() >= 1.0);
        }
        ARIADNE_TESTING_PROVE(full_inputs >= 4u);

        for (EventKind reset_kind: {EventKind::PERMISSIVE, EventKind::URGENT}) {
            ContinuousTimeHybridAutomaton system;
            system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                        {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
            system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
            system.new_mode(alternate, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
            system.new_invariant(source, x <= 2, DiscreteEvent("bound"));
            system.new_guard(source, progress, x >= 0, EventKind::PROGRESS);
            system.new_transition(source, first, target, {next(x) = x + 2},
                                  x + 1 >= 0, reset_kind);
            system.new_transition(source, second, alternate, {next(x) = x + 4},
                                  x + 1 >= 0, reset_kind);
            for (bool tree: {false, true}) {
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
                evolver.configuration().set_enable_reconditioning(false);
                if (tree)
                    evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
                const auto active = evolver.enclosure(HybridExactBox(
                    source, RealSpace({
                                x
                }),
                    ExactBoxType({{1, 1}})));
                // The carried enclosure contains a reached subset after full-input
                // propagation. Its PROGRESS guard is UNKNOWN; do not select a
                // numerical child or discard its independent input history.
                for (bool history: {false, true}) {
                    const auto& initial = history ? carried.final()[0u] : active;
                    const auto source_box = cast_exact_box(initial.state_bounding_box());
                    const auto source_domain = initial.continuous_set().parameter_domain();
                    const auto source_kinds = initial.continuous_set().variable_kinds();
                    const auto source_constraints = initial.number_of_constraints();
                    const auto orbit = evolver.orbit(
                        initial, HybridTime(0.25_x, 1), Semantics::LOWER);
                    ARIADNE_TESTING_REQUIRE(orbit.final().size() == 2u);
                    bool found_first = false, found_second = false;
                    for (const auto& final: orbit.final()) {
                        ARIADNE_TESTING_REQUIRE(final.previous_events().size() == 1u);
                        ARIADNE_TESTING_EXPECT_EQUAL(final.continuous_set().parameter_domain(), source_domain);
                        ARIADNE_TESTING_EXPECT_EQUAL(final.continuous_set().variable_kinds(), source_kinds);
                        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()),
                                                     cast_exact_interval(initial.time_range())));
                        const auto box = cast_exact_box(final.state_bounding_box());
                        if (final.location() == target) {
                            found_first = true;
                            ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events()[0u], first);
                            ARIADNE_TESTING_PROVE(subset(
                                history ? ExactBoxType({
                                              {1.875_x, 2.125_x}
                            })
                                        : ExactBoxType({{3, 3}}),
                                box));
                            ARIADNE_TESTING_PROVE(subset(box, history
                                                                  ? ExactBoxType({
                                                                        {1.8740234375_x, 2.1259765625_x}
                            })
                                                                  : ExactBoxType({{2.9990234375_x, 3.0009765625_x}})));
                        } else if (final.location() == alternate) {
                            found_second = true;
                            ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events()[0u], second);
                            ARIADNE_TESTING_PROVE(subset(
                                history ? ExactBoxType({
                                              {3.875_x, 4.125_x}
                            })
                                        : ExactBoxType({{5, 5}}),
                                box));
                            ARIADNE_TESTING_PROVE(subset(box, history
                                                                  ? ExactBoxType({
                                                                        {3.8740234375_x, 4.1259765625_x}
                            })
                                                                  : ExactBoxType({{4.9990234375_x, 5.0009765625_x}})));
                        } else {
                            ARIADNE_TESTING_REQUIRE(false);
                        }
                    }
                    ARIADNE_TESTING_PROVE(found_first && found_second);
                    ARIADNE_TESTING_EXPECT_EQUAL(initial.continuous_set().parameter_domain(), source_domain);
                    ARIADNE_TESTING_EXPECT_EQUAL(cast_exact_box(initial.state_bounding_box()), source_box);
                    ARIADNE_TESTING_EXPECT_EQUAL(initial.number_of_constraints(), source_constraints);
                }
            }
        }
    }

    Void test_lower_inclusion_permissive_reset_requires_whole_source_certificates() {
        const DiscreteLocation source(1), target(2);
        const RealVariable x("x"), u("u");
        const DiscreteEvent jump("jump"), progress("progress");
        const auto initial = HybridExactBox(
            source, RealSpace({
                        x
        }),
            ExactBoxType({{-0.125_x, 0.125_x}}));
        TaylorPicardIntegrator base(step_maximum_error = 1e-6);
        for (bool tree: {false, true}) {
            for (EventKind other: {EventKind::PERMISSIVE, EventKind::INVARIANT,
                                   EventKind::URGENT, EventKind::IMPACT}) {
                ContinuousTimeHybridAutomaton system;
                system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                            {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
                system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
                system.new_guard(source, progress, x >= 0, EventKind::PROGRESS);
                // Even a certain optional reset cannot discharge invariants,
                // uncertain URGENT guards or an IMPACT direction obligation.
                // An uncertain optional reset cannot discharge PROGRESS either.
                if (other == EventKind::PERMISSIVE) {
                    system.new_transition(source, jump, target, {next(x) = x + 2},
                                          x >= 0, EventKind::PERMISSIVE);
                } else {
                    system.new_transition(source, jump, target, {next(x) = x + 2},
                                          x + 1 >= 0, EventKind::PERMISSIVE);
                    if (other == EventKind::INVARIANT)
                        system.new_invariant(source, x <= 0, DiscreteEvent("bound"));
                    else
                        system.new_transition(source, DiscreteEvent("blocker"), target,
                                              {next(x) = x},
                                              other == EventKind::URGENT ? x >= 0 : x + 1 >= 0,
                                              other);
                }
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
                if (tree)
                    evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
                ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                              evolver.orbit(initial, HybridTime(0.25_x, 1), Semantics::LOWER));
            }
        }
    }

    Void test_lower_inclusion_certified_urgent_crossing() {
        const DiscreteLocation source(1), target(2);
        const RealVariable x("x"), u("u"), p("p");
        const TimeVariable t;
        const DiscreteEvent jump("crossing");
        for (Nat form: {0u, 1u, 2u, 3u, 4u}) {
            const bool timed = form == 1u || form >= 3u;
            const bool parameterized = form >= 2u;
            const bool joint_guard = form == 4u;
            const RealExpression drift = parameterized ? RealExpression(p) : RealExpression(1);
            const RealExpression clock = timed ? RealExpression(t) : RealExpression(0);
            const RealVariablesBox parameters = parameterized ? RealVariablesBox({1 <= p <= 2}) : RealVariablesBox();
            ContinuousTimeHybridAutomaton system;
            system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                        {dot(x) = drift + clock + u}, RealVariablesBox({Real(-0.125_x) <= u <= Real(0.125_x)}),
                                        parameters, timed ? std::optional<RealVariable>(t) : std::nullopt));
            system.new_mode(target, OrdinaryDifferentialAlgebraicEquation(
                                        {dot(x) = 0}, parameters, timed ? std::optional<RealVariable>(t) : std::nullopt));
            system.new_invariant(source, x <= 1, DiscreteEvent("bound"));
            system.new_transition(source, jump, target, {next(x) = x + 1},
                                  joint_guard ? x + t >= p / 64 : x >= Real(0.015625_x), EventKind::URGENT);
            for (bool tree: {false, true}) {
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                evolver.set_integrator(TaylorPicardIntegrator(step_maximum_error = 1e-7),
                                       ConstantApproximation(), ConstantApproximation());
                evolver.configuration().set_maximum_step_size(0.03125);
                evolver.configuration().set_enable_reconditioning(false);
                if (tree)
                    evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
                const auto initial = evolver.enclosure(HybridExactBox(
                    source, RealSpace({
                                x
                }),
                    ExactBoxType({{form == 0u ? -0.00390625_x : 0.0_x, 0}})));
                const auto source_domain = initial.continuous_set().parameter_domain();
                const auto source_box = cast_exact_box(initial.state_bounding_box());
                const auto orbit = evolver.orbit(initial, HybridTime(0.125_x, 1), Semantics::LOWER);
                ARIADNE_TESTING_REQUIRE(orbit.final().size() == 1u);
                const auto& final = orbit.final()[0u];
                ARIADNE_TESTING_EXPECT_EQUAL(final.location(), target);
                ARIADNE_TESTING_REQUIRE(final.previous_events().size() == 1u);
                ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events()[0u], jump);
                ARIADNE_TESTING_EXPECT_EQUAL(final.state_dimension(), 1u);
                // For the joint guard, p=1+1/256,u=0 has first contact
                // at t=1/128, x=t+t^2. Every actual source trajectory
                // is nevertheless certified, with full p and u domains.
                const ExactIntervalType witness = joint_guard
                                                    ? ExactIntervalType(1.00787353515625_x, 1.00787353515625_x)
                                                    : ExactIntervalType(1.015625_x, 1.015625_x);
                ARIADNE_TESTING_PROVE(subset(witness, cast_exact_interval(final.state_bounding_box()[0u])));
                ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()), ExactIntervalType(0, 0.03125_x)));
                if (joint_guard)
                    ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0.0078125_x, 0.0078125_x),
                                                 cast_exact_interval(final.time_range())));
                if (parameterized) {
                    ARIADNE_TESTING_EXPECT_EQUAL(final.continuous_set().constant_parameter_space(), RealSpace({p}));
                    ARIADNE_TESTING_PROVE(subset(ExactIntervalType(1, 2),
                                                 cast_exact_interval(final.continuous_set().constant_parameter_function(1u).get(0u).range())));
                }
                SizeType input_axes = 0u;
                const auto& continuous = final.continuous_set();
                for (SizeType i = 0u; i != continuous.number_of_parameters(); ++i) {
                    if (continuous.variable_kinds()[i] != EnclosureVariableKind::INPUT)
                        continue;
                    ++input_axes;
                    ARIADNE_TESTING_PROVE(continuous.parameter_domain()[i].lower_bound().raw().get_d() <= -0.125);
                    ARIADNE_TESTING_PROVE(continuous.parameter_domain()[i].upper_bound().raw().get_d() >= 0.125);
                }
                ARIADNE_TESTING_PROVE(input_axes > 0u);
                ARIADNE_TESTING_EXPECT_EQUAL(initial.continuous_set().parameter_domain(), source_domain);
                ARIADNE_TESTING_EXPECT_EQUAL(cast_exact_box(initial.state_bounding_box()), source_box);
                ARIADNE_TESTING_PROVE(initial.previous_events().empty());
            }
        }
    }

    Void test_lower_inclusion_crossing_requires_universal_certificates() {
        const DiscreteLocation source(1), target(2), alternate(3);
        const RealVariable x("x"), u("u");
        const DiscreteEvent jump("crossing"), other("other");
        for (bool tree: {false, true}) {
            for (EventKind competitor: {EventKind::URGENT, EventKind::PROGRESS,
                                        EventKind::INVARIANT, EventKind::IMPACT}) {
                ContinuousTimeHybridAutomaton system;
                system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                            {dot(x) = 1 + u}, RealVariablesBox({Real(-0.125_x) <= u <= Real(0.125_x)})));
                system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
                system.new_mode(alternate, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
                system.new_transition(source, jump, target, {next(x) = x + 1},
                                      x >= Real(0.015625_x), EventKind::URGENT);
                if (competitor == EventKind::INVARIANT)
                    system.new_invariant(source, x <= Real(0.0078125_x), other);
                else if (competitor == EventKind::PROGRESS)
                    system.new_guard(source, other, x >= Real(0.0078125_x), competitor);
                else
                    system.new_transition(source, other, alternate, {next(x) = x + 2},
                                          x >= Real(0.0078125_x), competitor);
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                evolver.set_integrator(TaylorPicardIntegrator(step_maximum_error = 1e-7),
                                       ConstantApproximation(), ConstantApproximation());
                evolver.configuration().set_maximum_step_size(0.03125);
                if (tree)
                    evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
                const auto initial = HybridExactBox(source, RealSpace({
                                                                x
                }),
                                                    ExactBoxType({{0, 0}}));
                ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                              evolver.orbit(initial, HybridTime(0.125_x, 1), Semantics::LOWER));
            }
            // Some inputs cross while others remain negative: an outer
            // candidate cannot establish a whole-source/input certificate.
            ContinuousTimeHybridAutomaton uncertain;
            uncertain.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                           {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
            uncertain.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
            uncertain.new_transition(source, jump, target, {next(x) = x + 1},
                                     x >= Real(0.015625_x), EventKind::URGENT);
            HybridContinuousTimeEvolver<ValidatedTag> evolver(uncertain);
            evolver.set_integrator(TaylorPicardIntegrator(step_maximum_error = 1e-7),
                                   ConstantApproximation(), ConstantApproximation());
            evolver.configuration().set_maximum_step_size(0.03125);
            if (tree)
                evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                          evolver.orbit(HybridExactBox(source, RealSpace({
                                                                                   x
            }),
                                                                       ExactBoxType({{0, 0}})),
                                                        HybridTime(0.125_x, 1), Semantics::LOWER));
        }
    }

    Void test_lower_urgent_crossing_respects_uncertain_arrival_horizon() {
        const DiscreteLocation source(1), target(2);
        const RealVariable x("x"), u("u");
        const TimeVariable t;
        const DiscreteEvent jump("crossing");
        for (bool inclusion: {false, true})
            for (bool tree: {false, true}) {
                ContinuousTimeHybridAutomaton system;
                if (inclusion)
                    system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                                {dot(x) = 1 + u}, RealVariablesBox({Real(-0.125_x) <= u <= Real(0.125_x)}),
                                                RealVariablesBox(), t));
                else
                    system.new_mode(source, OrdinaryDifferentialAlgebraicEquation(
                                                {dot(x) = 1}, RealVariablesBox(), t));
                system.new_mode(target, OrdinaryDifferentialAlgebraicEquation(
                                            {dot(x) = 0}, RealVariablesBox(), t));
                system.new_transition(source, jump, target, {next(x) = x + 1},
                                      x >= Real(0.03125_x), EventKind::URGENT);
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                TaylorPicardIntegrator base(step_maximum_error = 1e-7);
                if (inclusion)
                    evolver.set_integrator(base, ConstantApproximation(), ConstantApproximation());
                else
                    evolver.set_integrator(base);
                evolver.configuration().set_maximum_step_size(0.0625);
                if (tree)
                    evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
                auto initial = evolver.enclosure(HybridExactBox(source, RealSpace({
                                                                            x
                }),
                                                                ExactBoxType({{0, 0}})));
                // A genuine LOWER subset may live only at the late end of
                // this outer time range. Do not select early outer times.
                initial.set_time_function(evolver.function_factory().create_constant(
                                              initial.parameter_domain(), Real(0.03125_x)) +
                                          FloatDPBounds(-0.03125_x, 0.0_x, dp));
                if (inclusion) {
                    ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                                  evolver.orbit(initial, HybridTime(0.046875_x, 1), Semantics::LOWER));
                } else {
                    const auto uncertified = evolver.orbit(initial, HybridTime(0.046875_x, 1), Semantics::LOWER);
                    ARIADNE_TESTING_EXPECT_EQUAL(uncertified.final().size(), 0u);
                }
                // Uncertain arrival is supported when the complete step
                // endpoint, hence every first contact, is within the horizon.
                const auto certified = evolver.orbit(initial, HybridTime(0.125_x, 1), Semantics::LOWER);
                ARIADNE_TESTING_REQUIRE(certified.final().size() == 1u);
                ARIADNE_TESTING_EXPECT_EQUAL(certified.final()[0u].location(), target);
                ARIADNE_TESTING_REQUIRE(certified.final()[0u].previous_events().size() == 1u);
                ARIADNE_TESTING_EXPECT_EQUAL(certified.final()[0u].previous_events()[0u], jump);
                ARIADNE_TESTING_PROVE(subset(cast_exact_interval(certified.final()[0u].time_range()),
                                             ExactIntervalType(0, 0.125_x)));
            }
    }

    Void test_lower_permissive_endpoints_keep_flow_and_all_certain_resets() {
        const DiscreteLocation source(1), first_target(2), second_target(3);
        const RealVariable x("x"), u("u"), p("p"), output("output");
        const TimeVariable t;
        const DiscreteEvent first("first"), second("second");
        const Dyadic duration(0.03125_x);
        for (Nat form: {0u, 1u, 2u, 3u, 4u, 5u, 6u}) {
            const bool inclusion = form < 4u;
            const bool timed = form == 1u || form == 3u || form == 4u || form == 6u;
            const bool parameterized = form == 2u || form == 3u || form >= 5u;
            const RealExpression drift = parameterized ? RealExpression(p) : RealExpression(1);
            const RealExpression clock = timed ? RealExpression(t) : RealExpression(0);
            const RealVariablesBox parameters = parameterized ? RealVariablesBox({1 <= p <= 2}) : RealVariablesBox();
            const auto time = timed ? std::optional<RealVariable>(t) : std::nullopt;
            ContinuousTimeHybridAutomaton system;
            if (inclusion)
                system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                            {dot(x) = drift + clock + u}, RealVariablesBox({Real(-0.125_x) <= u <= Real(0.125_x)}), parameters, time));
            else
                system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = drift + clock}, parameters, time));
            for (const auto& target: {first_target, second_target}) {
                system.new_mode(target, OrdinaryDifferentialAlgebraicEquation(
                                            {dot(x) = 0}, parameters, time, {let(output) = x + drift + clock}));
                system.new_invariant(target, x <= 10, DiscreteEvent("target_bound"));
            }
            system.new_invariant(source, x <= 1, DiscreteEvent("source_bound"));
            system.new_transition(source, first, first_target, {next(x) = x + 2},
                                  parameterized ? (timed ? x + t >= p / 64 : x >= p / 64) : x >= Real(0.015625_x), EventKind::PERMISSIVE);
            system.new_transition(source, second, second_target, {next(x) = x + 4},
                                  parameterized ? x >= p / 128 : x >= Real(0.0078125_x), EventKind::PERMISSIVE);
            for (bool tree: {false, true}) {
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                TaylorPicardIntegrator base(step_maximum_error = 1e-7);
                if (inclusion)
                    evolver.set_integrator(base, ConstantApproximation(), ConstantApproximation());
                else
                    evolver.set_integrator(base);
                evolver.configuration().set_maximum_step_size(0.03125);
                evolver.configuration().set_enable_reconditioning(false);
                if (tree)
                    evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
                const auto initial = evolver.enclosure(HybridExactBox(source, RealSpace({
                                                                                  x
                }),
                                                                      ExactBoxType({{0, 0}})));
                const auto original_domain = initial.continuous_set().parameter_domain();
                const auto orbit = evolver.orbit(initial, HybridTime(0.03125_x, 1), Semantics::LOWER);
                ARIADNE_TESTING_REQUIRE(orbit.final().size() == 3u);
                const HybridEnclosure* flowing = nullptr;
                for (const auto& final: orbit.final())
                    if (final.location() == source)
                        flowing = &final;
                ARIADNE_TESTING_REQUIRE(flowing != nullptr);
                ARIADNE_TESTING_PROVE(flowing->previous_events().empty());
                const Dyadic input_radius(inclusion ? 0.00390625_x : 0.0_x);
                const Dyadic clock_shift(timed ? 0.00048828125_x : 0.0_x);
                ARIADNE_TESTING_PROVE(subset(ExactIntervalType(duration - input_radius + clock_shift,
                                                               (parameterized ? 2 * duration : duration) + input_radius + clock_shift),
                                             cast_exact_interval(flowing->state_bounding_box()[0u])));
                bool found_first = false, found_second = false;
                for (const auto& final: orbit.final()) {
                    ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()), ExactIntervalType(duration, duration)));
                    ARIADNE_TESTING_EXPECT_EQUAL(final.continuous_set().parameter_domain(), flowing->continuous_set().parameter_domain());
                    ARIADNE_TESTING_EXPECT_EQUAL(final.continuous_set().variable_kinds(), flowing->continuous_set().variable_kinds());
                    if (parameterized) {
                        ARIADNE_TESTING_EXPECT_EQUAL(final.continuous_set().constant_parameter_space(), RealSpace({p}));
                        ARIADNE_TESTING_PROVE(subset(ExactIntervalType(1, 2),
                                                     cast_exact_interval(final.continuous_set().constant_parameter_function(1u).get(0u).range())));
                    }
                    if (final.location() == source)
                        continue;
                    ARIADNE_TESTING_REQUIRE(final.previous_events().size() == 1u);
                    const bool first_branch = final.location() == first_target;
                    ARIADNE_TESTING_PROVE(first_branch || final.location() == second_target);
                    found_first = found_first || first_branch;
                    found_second = found_second || !first_branch;
                    ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events()[0u], first_branch ? first : second);
                    const auto residual = final.state_function()[0u] - flowing->state_function()[0u] - (first_branch ? 2 : 4);
                    ARIADNE_TESTING_PROVE(subset(cast_exact_interval(residual.range()), ExactIntervalType(-0.00390625_x, 0.00390625_x)));
                    const auto physical = parameterized
                                            ? final.continuous_set().constant_parameter_function(1u).get(0u)
                                            : evolver.function_factory().create_constant(final.parameter_domain(), Real(1));
                    const auto auxiliary_residual = final.auxiliary_function()[0u] - final.state_function()[0u] - physical -
                                                    (timed ? final.time_function() : 0 * final.time_function());
                    ARIADNE_TESTING_PROVE(subset(cast_exact_interval(auxiliary_residual.range()), ExactIntervalType(-0.00390625_x, 0.00390625_x)));
                }
                ARIADNE_TESTING_PROVE(found_first && found_second);
                ARIADNE_TESTING_EXPECT_EQUAL(initial.continuous_set().parameter_domain(), original_domain);
                ARIADNE_TESTING_PROVE(initial.previous_events().empty());
            }
        }
    }

    Void test_lower_permissive_endpoint_requires_guard_and_target_certificates() {
        const DiscreteLocation source(1), target(2);
        const RealVariable x("x"), u("u");
        const DiscreteEvent jump("jump");
        for (bool tree: {false, true})
            for (Nat condition: {0u, 1u, 2u}) {
                ContinuousTimeHybridAutomaton system;
                system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                            {dot(x) = condition == 0u ? RealExpression(u) : 1 + u},
                                            RealVariablesBox({Real(-0.125_x) <= u <= Real(0.125_x)})));
                system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
                system.new_transition(source, jump, target, {next(x) = x + 2},
                                      condition == 0u ? x >= 0 : x >= Real(0.015625_x), EventKind::PERMISSIVE);
                if (condition != 0u)
                    system.new_invariant(target, x <= (condition == 1u ? Real(2.03125_x) : Real(1)), DiscreteEvent("bound"));
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                evolver.set_integrator(TaylorPicardIntegrator(step_maximum_error = 1e-7),
                                       ConstantApproximation(), ConstantApproximation());
                evolver.configuration().set_maximum_step_size(0.03125);
                if (tree)
                    evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
                const auto orbit = evolver.orbit(HybridExactBox(source, RealSpace({
                                                                            x
                }),
                                                                ExactBoxType({{0, 0}})),
                                                 HybridTime(0.03125_x, 1), Semantics::LOWER);
                ARIADNE_TESTING_REQUIRE(orbit.final().size() == 1u);
                ARIADNE_TESTING_EXPECT_EQUAL(orbit.final()[0u].location(), source);
                ARIADNE_TESTING_PROVE(orbit.final()[0u].previous_events().empty());
            }
    }

    Void test_lower_permissive_endpoint_respects_prefix_and_horizon() {
        const DiscreteLocation source(1), target(2);
        const RealVariable x("x"), u("u");
        const TimeVariable t;
        const DiscreteEvent jump("jump");
        for (bool tree: {false, true}) {
            for (bool inclusion: {false, true}) {
                ContinuousTimeHybridAutomaton system;
                if (inclusion)
                    system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                                {dot(x) = 1 + u}, RealVariablesBox({Real(-0.125_x) <= u <= Real(0.125_x)}), RealVariablesBox(), t));
                else
                    system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = 1}, RealVariablesBox(), t));
                system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, RealVariablesBox(), t));
                system.new_transition(source, jump, target, {next(x) = x + 2},
                                      x >= Real(0.015625_x), EventKind::PERMISSIVE);
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                TaylorPicardIntegrator base(step_maximum_error = 1e-7);
                if (inclusion)
                    evolver.set_integrator(base, ConstantApproximation(), ConstantApproximation());
                else
                    evolver.set_integrator(base);
                evolver.configuration().set_maximum_step_size(0.0625);
                if (tree)
                    evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
                auto initial = evolver.enclosure(HybridExactBox(source, RealSpace({
                                                                            x
                }),
                                                                ExactBoxType({{0, 0}})));
                initial.set_time_function(evolver.function_factory().create_constant(initial.parameter_domain(), Real(0.03125_x)) +
                                          FloatDPBounds(-0.03125_x, 0.0_x, dp));
                const auto orbit = evolver.orbit(initial, HybridTime(0.046875_x, 1), Semantics::LOWER);
                ARIADNE_TESTING_REQUIRE(orbit.final().size() == 1u);
                ARIADNE_TESTING_EXPECT_EQUAL(orbit.final()[0u].location(), source);
                ARIADNE_TESTING_PROVE(orbit.final()[0u].previous_events().empty());
            }
            // A certain optional reset at the endpoint does not license
            // flowing through an earlier blocking condition.
            ContinuousTimeHybridAutomaton blocked;
            blocked.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                         {dot(x) = 1 + u}, RealVariablesBox({Real(-0.125_x) <= u <= Real(0.125_x)})));
            blocked.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
            blocked.new_guard(source, DiscreteEvent("progress"), x >= Real(0.015625_x), EventKind::PROGRESS);
            blocked.new_transition(source, jump, target, {next(x) = x + 2},
                                   x >= Real(0.0078125_x), EventKind::PERMISSIVE);
            HybridContinuousTimeEvolver<ValidatedTag> evolver(blocked);
            evolver.set_integrator(TaylorPicardIntegrator(step_maximum_error = 1e-7), ConstantApproximation(), ConstantApproximation());
            evolver.configuration().set_maximum_step_size(0.03125);
            if (tree)
                evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                          evolver.orbit(HybridExactBox(source, RealSpace({
                                                                                   x
            }),
                                                                       ExactBoxType({{0, 0}})),
                                                        HybridTime(0.03125_x, 1), Semantics::LOWER));
        }
    }

    Void test_lower_permissive_endpoint_resets_preserve_later_continuations() {
        const DiscreteLocation source(1), target(2);
        const RealVariable x("x"), u("u");
        const DiscreteEvent jump("jump");
        ContinuousTimeHybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                    {dot(x) = 1 + u}, RealVariablesBox({Real(-0.125_x) <= u <= Real(0.125_x)})));
        system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
        system.new_transition(source, jump, target, {next(x) = x + 2},
                              x >= Real(0.015625_x), EventKind::PERMISSIVE);
        for (bool tree: {false, true}) {
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            evolver.set_integrator(TaylorPicardIntegrator(step_maximum_error = 1e-7),
                                   ConstantApproximation(), ConstantApproximation());
            evolver.configuration().set_maximum_step_size(0.03125);
            evolver.configuration().set_enable_reconditioning(false);
            if (tree)
                evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
            const auto orbit = evolver.orbit(HybridExactBox(source, RealSpace({
                                                                        x
            }),
                                                            ExactBoxType({{0, 0}})),
                                             HybridTime(0.0625_x, 2), Semantics::LOWER);
            ARIADNE_TESTING_REQUIRE(orbit.final().size() == 3u);
            bool source_final = false, early_reset = false, late_reset = false;
            for (const auto& final: orbit.final()) {
                ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()), ExactIntervalType(0.0625_x, 0.0625_x)));
                SizeType inputs = 0u;
                const auto& continuous = final.continuous_set();
                for (SizeType i = 0u; i != continuous.number_of_parameters(); ++i) {
                    if (continuous.variable_kinds()[i] != EnclosureVariableKind::INPUT)
                        continue;
                    ++inputs;
                    ARIADNE_TESTING_PROVE(continuous.parameter_domain()[i].lower_bound().raw().get_d() <= -0.125);
                    ARIADNE_TESTING_PROVE(continuous.parameter_domain()[i].upper_bound().raw().get_d() >= 0.125);
                }
                const auto box = cast_exact_interval(final.state_bounding_box()[0u]);
                if (final.location() == source) {
                    source_final = true;
                    ARIADNE_TESTING_PROVE(final.previous_events().empty());
                    ARIADNE_TESTING_EXPECT_EQUAL(inputs, 2u);
                    ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0.0546875_x, 0.0703125_x), box));
                } else {
                    ARIADNE_TESTING_EXPECT_EQUAL(final.location(), target);
                    ARIADNE_TESTING_REQUIRE(final.previous_events().size() == 1u);
                    ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events()[0u], jump);
                    if (inputs == 1u) {
                        early_reset = true;
                        ARIADNE_TESTING_PROVE(subset(ExactIntervalType(2.02734375_x, 2.03515625_x), box));
                        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.dwell_time_range()), ExactIntervalType(0.03125_x, 0.03125_x)));
                    } else {
                        late_reset = true;
                        ARIADNE_TESTING_EXPECT_EQUAL(inputs, 2u);
                        ARIADNE_TESTING_PROVE(subset(ExactIntervalType(2.0546875_x, 2.0703125_x), box));
                        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.dwell_time_range()), ExactIntervalType(0, 0)));
                    }
                }
            }
            ARIADNE_TESTING_PROVE(source_final && early_reset && late_reset);
        }
    }

    Void test_upper_competing_inclusion_events_cover_clock_parameters_and_targets() {
        const DiscreteLocation source(1), first_target(2), second_target(3), optional_target(4);
        const RealVariable x("x"), u("u"), w("w"), p("p"), output("output");
        const TimeVariable t;
        const DiscreteEvent first("first"), second("second"), optional("optional");
        const RealVariablesBox parameters({1 <= p <= 2});
        ContinuousTimeHybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                    {dot(x) = p + t + u}, RealVariablesBox({Real(-0.125_x) <= u <= Real(0.125_x)}), parameters, t));
        system.new_mode(first_target, OrdinaryDifferentialAlgebraicEquation(
                                          {dot(x) = 0}, parameters, t, {let(output) = x - p}));
        system.new_mode(second_target, OrdinaryDifferentialAlgebraicInclusion(
                                           {dot(x) = w}, RealVariablesBox({Real(-0.0625_x) <= w <= Real(0.0625_x)}), parameters, t));
        system.new_mode(optional_target, OrdinaryDifferentialAlgebraicEquation(
                                             {dot(x) = 0}, parameters, t));
        system.new_transition(source, first, first_target, {next(x) = x + p},
                              x >= p / 64, EventKind::URGENT);
        system.new_transition(source, second, second_target, {next(x) = x + 2 * p},
                              x >= Real(0.0234375_x), EventKind::URGENT);
        system.new_transition(source, optional, optional_target, {next(x) = x - p},
                              x >= Real(0.0078125_x), EventKind::PERMISSIVE);
        system.new_invariant(source, x <= p / 8, DiscreteEvent("source_bound"));
        for (const auto& target: {first_target, second_target, optional_target})
            system.new_invariant(target, p <= 2, DiscreteEvent("parameter_bound"));
        for (bool boxes: {false, true})
            for (bool tree: {false, true}) {
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                if (boxes)
                    evolver.set_integrator(FirstOrderMultifunctionIntegrator());
                else
                    evolver.set_integrator(TaylorPicardIntegrator(step_maximum_error = 1e-7),
                                           ConstantApproximation(), ConstantApproximation());
                evolver.configuration().set_maximum_step_size(0.015625);
                if (tree)
                    evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
                auto initial = evolver.enclosure(HybridExactBox(source, RealSpace({
                                                                            x
                }),
                                                                ExactBoxType({{0, 0}})));
                initial.set_time_function(evolver.function_factory().create_constant(initial.parameter_domain(), Real(2)));
                const auto source_domain = initial.continuous_set().parameter_domain();
                const auto orbit = evolver.orbit(initial, HybridTime(2.125_x, 2), Semantics::UPPER);
                std::optional<ExactIntervalType> first_states, second_states, optional_states;
                for (const auto& final: orbit.final()) {
                    ARIADNE_TESTING_EXPECT_EQUAL(final.state_dimension(), 1u);
                    ARIADNE_TESTING_REQUIRE(final.previous_events().size() == 1u);
                    ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()), ExactIntervalType(2.125_x, 2.125_x)));
                    ARIADNE_TESTING_EXPECT_EQUAL(final.continuous_set().constant_parameter_space(), RealSpace({p}));
                    ARIADNE_TESTING_EXPECT_EQUAL(final.continuous_set().constant_parameter_domain(), ExactBoxType({
                                                                                                         {1, 2}
                    }));
                    const auto state = cast_exact_interval(final.state_bounding_box()[0u]);
                    if (final.location() == first_target) {
                        ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events()[0u], first);
                        first_states = first_states ? hull(*first_states, state) : state;
                        const auto physical = final.continuous_set().constant_parameter_function(1u).get(0u);
                        const auto residual = final.auxiliary_function()[0u] - final.state_function()[0u] + physical;
                        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(residual.range()), ExactIntervalType(-0.015625_x, 0.015625_x)));
                    } else if (final.location() == second_target) {
                        ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events()[0u], second);
                        second_states = second_states ? hull(*second_states, state) : state;
                    } else {
                        ARIADNE_TESTING_EXPECT_EQUAL(final.location(), optional_target);
                        ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events()[0u], optional);
                        optional_states = optional_states ? hull(*optional_states, state) : state;
                    }
                }
                ARIADNE_TESTING_REQUIRE(first_states.has_value() && second_states.has_value() && optional_states.has_value());
                // u=0 reaches the first urgent boundary first for p=1,
                // and the second boundary first for p=2. The optional
                // guard is reached before either, for both parameter values.
                ARIADNE_TESTING_PROVE(subset(ExactIntervalType(1.015625_x, 1.015625_x), *first_states));
                ARIADNE_TESTING_PROVE(subset(ExactIntervalType(4.0234375_x, 4.0234375_x), *second_states));
                ARIADNE_TESTING_PROVE(subset(ExactIntervalType(-1.9921875_x, -0.9921875_x), *optional_states));
                ARIADNE_TESTING_EXPECT_EQUAL(initial.continuous_set().parameter_domain(), source_domain);
                ARIADNE_TESTING_PROVE(initial.previous_events().empty());
            }
    }

    Void test_upper_interval_taylor_events_cover_urgent_and_optional_branches() {
        const DiscreteLocation source(1), urgent_target(2), optional_target(3);
        const RealVariable x("x"), w("w");
        const DiscreteEvent urgent("urgent"), optional("optional");
        const ExactBoxType domain({
            {-1, 1}
        });
        const auto field = ValidatedIntervalTaylorFunctionModel<FloatDP>::constant(
            domain, FloatDPUpperInterval(0.875_x, 1.125_x, dp), ThresholdSweeper<FloatDP>(dp, 1e-10));
        ContinuousTimeHybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(RealSpace({x}),
                                                                       ValidatedVectorIntervalTaylorFunctionModel<FloatDP>({field})));
        system.new_mode(urgent_target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
        system.new_mode(optional_target, OrdinaryDifferentialAlgebraicInclusion(
                                             {dot(x) = w}, RealVariablesBox({Real(-0.0625_x) <= w <= Real(0.0625_x)})));
        system.new_transition(source, urgent, urgent_target, {next(x) = x + 1},
                              x >= Real(0.015625_x), EventKind::URGENT);
        system.new_transition(source, optional, optional_target, {next(x) = x - 1},
                              x >= Real(0.0078125_x), EventKind::PERMISSIVE);
        for (bool boxes: {false, true})
            for (bool tree: {false, true}) {
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                if (boxes)
                    evolver.set_integrator(FirstOrderMultifunctionIntegrator());
                else
                    evolver.set_integrator(TaylorPicardIntegrator(step_maximum_error = 1e-7),
                                           ConstantApproximation(), ConstantApproximation());
                evolver.configuration().set_maximum_step_size(0.03125);
                if (tree)
                    evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
                const auto orbit = evolver.orbit(HybridExactBox(source, RealSpace({
                                                                            x
                }),
                                                                ExactBoxType({{0, 0}})),
                                                 HybridTime(0.125_x, 2), Semantics::UPPER);
                std::optional<ExactIntervalType> urgent_states, optional_states;
                for (const auto& final: orbit.final()) {
                    ARIADNE_TESTING_REQUIRE(final.previous_events().size() == 1u);
                    ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()), ExactIntervalType(0.125_x, 0.125_x)));
                    const auto state = cast_exact_interval(final.state_bounding_box()[0u]);
                    if (final.location() == urgent_target) {
                        ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events()[0u], urgent);
                        urgent_states = urgent_states ? hull(*urgent_states, state) : state;
                    } else {
                        ARIADNE_TESTING_EXPECT_EQUAL(final.location(), optional_target);
                        ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events()[0u], optional);
                        optional_states = optional_states ? hull(*optional_states, state) : state;
                    }
                }
                ARIADNE_TESTING_REQUIRE(urgent_states.has_value() && optional_states.has_value());
                ARIADNE_TESTING_PROVE(subset(ExactIntervalType(1.015625_x, 1.015625_x), *urgent_states));
                ARIADNE_TESTING_PROVE(subset(ExactIntervalType(-0.9921875_x, -0.9921875_x), *optional_states));
            }
    }

    Void test_unsupported_capabilities_fail_explicitly() {
        const DiscreteLocation q(1);
        const RealVariable x("x"), u("u"), t("t"), p("p");
        const auto initial = HybridExactBox(q, RealSpace({
                                                   x
        }),
                                            ExactBoxType({{0, 0}}));
        ContinuousTimeHybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicInclusion({dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, evolver.orbit(initial, HybridTime(0u, 1), Semantics::UPPER));
        evolver.set_integrator(FirstOrderMultifunctionIntegrator());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, evolver.orbit(initial, HybridTime(0u, 1), Semantics::LOWER));
        evolver.set_integrator(TaylorPicardIntegrator(step_maximum_error = 1e-6));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, evolver.orbit(initial, HybridTime(0.125_x, 1), Semantics::UPPER));
        // State-only invariants are now supported as necessary outer
        // restrictions. They need neither an entry time nor a Lie derivative.
        system.new_invariant(q, x <= 1, DiscreteEvent("bound"));
        HybridContinuousTimeEvolver<ValidatedTag> guarded(system);
        guarded.set_integrator(FirstOrderMultifunctionIntegrator());
        const auto initially_valid = guarded.orbit(initial, HybridTime(0u, 1), Semantics::UPPER);
        ARIADNE_TESTING_PROVE(initially_valid.final().size() > 0u);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, guarded.orbit(initial, HybridTime(0u, 1), Semantics::LOWER));
        // A state invariant does not establish a certified entry time or
        // LOWER reset capability; possible UPPER jumps use guarded reach.
        const DiscreteLocation target(2);
        ContinuousTimeHybridAutomaton jumping;
        jumping.new_mode(q, OrdinaryDifferentialAlgebraicInclusion(
                                {dot(x) = u}, RealVariablesBox({-1 <= u <= 1})));
        jumping.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
        jumping.new_transition(q, DiscreteEvent("jump"), target, {next(x) = x},
                               x >= 0, EventKind::URGENT);
        HybridContinuousTimeEvolver<ValidatedTag> possible_jump(jumping);
        possible_jump.set_integrator(FirstOrderMultifunctionIntegrator());
        const auto at_zero = possible_jump.orbit(initial, HybridTime(0u, 1), Semantics::UPPER);
        // The initial urgent phase is allowed in UPPER. Its result is a
        // conservative event candidate, never a coverage certificate.
        ARIADNE_TESTING_PROVE(at_zero.final().size() > 0u);

        ContinuousTimeHybridAutomaton timed;
        timed.new_mode(q, OrdinaryDifferentialAlgebraicInclusion({dot(x) = t + p + u}, RealVariablesBox({-1 <= u <= 1}),
                                                                 RealVariablesBox({0 <= p <= 1}), t));
        HybridContinuousTimeEvolver<ValidatedTag> timed_evolver(timed);
        timed_evolver.set_integrator(FirstOrderMultifunctionIntegrator());
        const auto timed_parameters = timed_evolver.orbit(initial, HybridTime(0u, 1), Semantics::UPPER);
        ARIADNE_TESTING_REQUIRE(timed_parameters.final().size() == 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(timed_parameters.final()[0u].continuous_set().constant_parameter_space(), RealSpace({p}));
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(timed_parameters.final()[0u].continuous_set().constant_parameter_function(1u).get(0u).range()), ExactIntervalType(0, 1)));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, timed_evolver.orbit(initial, HybridTime(0u, 1), Semantics::LOWER));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_inclusions");
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_autonomous_inclusion_with_both_explicit_solver_choices);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_reconditioning_protects_carried_parameters);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_reset_into_inclusion_with_varying_arrival_time);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_interval_taylor_modes_use_the_same_cycle);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_box_horizon_slice_and_copy_ownership);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_upper_inclusion_transitions_use_guarded_reach);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_state_invariants_filter_inclusion_reach_and_endpoints);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_inclusion_invariants_preserve_closed_zero_horizon);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_upper_inclusion_progress_blocks_flow_without_reset);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_upper_inclusion_competing_urgent_guards_keep_both_boundary_jumps);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_upper_inclusion_permissive_guard_does_not_block_flow);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_structured_lower_inclusion_preserves_complete_input_flow);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_lower_inclusion_ignores_optional_jump_but_rejects_blockers);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_lower_inclusion_certified_whole_tube_invariants);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_lower_inclusion_certainly_inactive_blockers);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_lower_inclusion_full_source_urgent_reset_without_flow);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_lower_inclusion_certified_reset_beside_uncertain_progress);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_lower_inclusion_permissive_reset_requires_whole_source_certificates);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_lower_inclusion_certified_urgent_crossing);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_lower_inclusion_crossing_requires_universal_certificates);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_lower_urgent_crossing_respects_uncertain_arrival_horizon);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_lower_permissive_endpoints_keep_flow_and_all_certain_resets);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_lower_permissive_endpoint_requires_guard_and_target_certificates);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_lower_permissive_endpoint_respects_prefix_and_horizon);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_lower_permissive_endpoint_resets_preserve_later_continuations);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_upper_competing_inclusion_events_cover_clock_parameters_and_targets);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_upper_interval_taylor_events_cover_urgent_and_optional_branches);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusions, test_unsupported_capabilities_fail_explicitly);
    return Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
