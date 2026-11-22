/***************************************************************************
 *            test_hybrid_step_contracts.cpp
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

/** @file tests/test_hybrid_step_contracts.cpp
 *  @brief Test hybrid step contracts behavior and numerical invariants.
 */

#include "evolution/hybrid_continuous_time_evolver.hpp"
#include "integration/integrator.hpp"
#include "integration/validated_integration_step.hpp"
#include "log_verbosity_scope.hpp"
#include "logging/test_runner.hpp"
#include "system/hybrid_automata.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

class StepProbe : public HybridContinuousTimeEvolver {
  public:
    using HybridContinuousTimeEvolver::_estimate_timing;
    using HybridContinuousTimeEvolver::HybridContinuousTimeEvolver;
    using HybridContinuousTimeEvolverBase::_apply_evolve_step;
    using HybridContinuousTimeEvolverBase::_apply_guard;
    using HybridContinuousTimeEvolverBase::_apply_guard_step;
    using HybridContinuousTimeEvolverBase::_compute_integration_step;
};

// Both certificates enclose x'=1, but the tube deliberately carries a wider
// uniform error. Its endpoint slice must not replace the independent endpoint.
class CertificateProbe : public HybridContinuousTimeEvolver {
  public:
    enum class Shape { STRUCTURED,
                       MULTIPATCH,
                       BOX_ENDPOINT,
                       INPUT_COORDINATE };
    Shape shape = Shape::STRUCTURED;
    mutable std::optional<ExactBoxType> certified_domain;
    using HybridContinuousTimeEvolver::HybridContinuousTimeEvolver;
    using HybridContinuousTimeEvolverBase::_evolution_step;

  protected:
    ValidatedIntegrationStep _compute_integration_step(EffectiveVectorMultivariateFunction,
                                                       const ExactBoxType& initial,
                                                       const StepSizeType& maximum_step) const override {
        certified_domain = initial;
        const auto h = min(maximum_step, StepSizeType(0.25_x));
        auto domain = product(initial, ExactIntervalType(0u, h));
        if (shape == Shape::INPUT_COORDINATE) {
            domain = product(domain, ExactIntervalType(-1, 1));
        }
        const auto factory = function_factory();
        const auto tube = ValidatedVectorMultivariateFunctionPatch(1u,
                                                                   factory.create_coordinate(domain, 0u) + factory.create_coordinate(domain, 1u) +
                                                                       FloatDPBounds(-0.0625_x, 0.0625_x, dp));
        const auto endpoint = ValidatedVectorMultivariateFunctionPatch(1u,
                                                                       factory.create_coordinate(initial, 0u) + Real(h));
        List<ValidatedVectorMultivariateFunctionPatch> reach{tube};
        if (shape == Shape::MULTIPATCH) {
            reach.append(tube);
        }
        return ValidatedIntegrationStep(h, tube.range(), endpoint.range(), {}, reach,
                                        shape == Shape::BOX_ENDPOINT ? std::nullopt : std::optional(endpoint), StepSizeType(0u), 1u);
    }
};

class TestHybridStepContracts : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    void test_full_steps_use_independent_endpoint_and_keep_hybrid_metadata() {
        const DiscreteLocation q(1);
        const RealVariable x("x"), output("output");
        const RealSpace space({x});
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1}, {let(output) = 2 * x});
        CertificateProbe probe(system);
        probe.configuration().set_maximum_step_size(0.25);
        auto initial = probe.enclosure(HybridExactBox(q, space, ExactBoxType({
                                                                    {0, 0}
        })));
        initial.set_time_function(probe.function_factory().create_constant(initial.parameter_domain(), Real(2)));
        const DiscreteEvent prior("prior");
        initial.apply_reset(prior, q, space, make_function(space, Vector<RealExpression>{x}),
                            RealSpace({output}), system.auxiliary_function(q));
        const auto orbit = probe.orbit(initial, HybridTime(2.5_x, 3), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.intermediate().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.reach().size(), 2u);
        for (const auto& final: orbit.final()) {
            const auto box = cast_exact_box(final.state_bounding_box());
            ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                             {0.5_x, 0.5_x}
            }),
                                         box));
            ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                                  {0.4990234375_x, 0.5009765625_x}
            })));
            ARIADNE_TESTING_EXPECT_EQUAL(final.location(), q);
            ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().size(), 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().front(), prior);
            ARIADNE_TESTING_EXPECT_EQUAL(final.state_auxiliary_set().euclidean_set().dimension(), 2u);
            ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()), ExactIntervalType(2.5_x, 2.5_x)));
            ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.dwell_time_range()), ExactIntervalType(0.5_x, 0.5_x)));
        }
    }

    void test_shortened_steps_use_the_tube_at_the_stopping_time() {
        const DiscreteLocation q(1);
        const RealVariable x("x");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1});
        CertificateProbe probe(system);
        probe.configuration().set_maximum_step_size(0.25);
        const auto orbit = probe.orbit(HybridExactBox(q, RealSpace({
                                                             x
        }),
                                                      ExactBoxType({{0, 0}})),
                                       HybridTime(0.125_x, 1), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 1u);
        for (const auto& final: orbit.final()) {
            const auto box = cast_exact_box(final.state_bounding_box());
            ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                             {0.125_x, 0.125_x}
            }),
                                         box));
            ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                                  {0.0625_x, 0.1875_x}
            })));
            ARIADNE_TESTING_PROVE(subset(cast_exact_interval(final.time_range()), ExactIntervalType(0.125_x, 0.125_x)));
        }
    }

    void test_constrained_source_stays_inside_the_certified_composition_domain() {
        const DiscreteLocation q(1);
        const RealVariable x("x");
        const RealSpace space({x});
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1});
        CertificateProbe probe(system);
        probe.configuration().set_maximum_step_size(0.25);
        auto initial = probe.enclosure(HybridExactBox(q, space, ExactBoxType({
                                                                    {-0.125_x, 0.125_x}
        })));
        initial.new_state_constraint(DiscreteEvent("nonnegative"), ValidatedConstraint(0, make_function(space, x), 1));
        const auto orbit = probe.orbit(initial, HybridTime(0.25_x, 1), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 1u);
        ARIADNE_TESTING_REQUIRE(probe.certified_domain.has_value());
        // Even infeasible parameter values remain in the composition domain.
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {-0.125_x, 0.125_x}
        }),
                                     *probe.certified_domain));
        for (const auto& final: orbit.final()) {
            // The original state constraint and the final-time equality remain.
            ARIADNE_TESTING_EXPECT_EQUAL(final.number_of_constraints(), initial.number_of_constraints() + 1u);
            ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                             {0.25_x, 0.375_x}
            }),
                                         cast_exact_box(final.state_bounding_box())));
        }
    }

    void test_unsupported_step_shapes_are_rejected_before_publication() {
        const DiscreteLocation q(1);
        const RealVariable x("x");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1});
        CertificateProbe probe(system);
        const auto initial = probe.enclosure(HybridExactBox(q, RealSpace({
                                                                   x
        }),
                                                            ExactBoxType({{0, 0}})));
        for (auto shape: {CertificateProbe::Shape::MULTIPATCH, CertificateProbe::Shape::BOX_ENDPOINT,
                          CertificateProbe::Shape::INPUT_COORDINATE}) {
            probe.shape = shape;
            HybridEvolutionStepResult result;
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                          probe._evolution_step(initial, Real(1), Semantics::UPPER, result));
            ARIADNE_TESTING_PROVE(result.contributions.empty());
            ARIADNE_TESTING_PROVE(result.successors.empty());
            ARIADNE_TESTING_PROVE(subset(cast_exact_box(initial.state_bounding_box()), ExactBoxType({
                                                                                           {0, 0}
            })));
        }
    }

    void test_constant_and_decreasing_guard_steps() {
        DiscreteLocation q(1);
        DiscreteEvent event("guard");
        RealVariable x("x"), y("y");
        RealSpace space({x, y});
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0});
        StepProbe probe(system);
        auto initial = probe.enclosure(HybridExactBox(q, space, ExactBoxType({
                                                                    {-1, 1},
                                                                    {-1, 1}
        })));
        auto dynamic = system.dynamic_function(q);
        FlowFunctionModel flow = probe._compute_integration_step(dynamic, cast_exact_box(initial.state_bounding_box()), Dyadic(0.5_x)).reach_models()[0];
        Map<DiscreteEvent, CrossingData> crossings;
        Map<DiscreteEvent, TransitionData> transitions;
        Set<DiscreteEvent> active;
        auto timing = probe._estimate_timing(active, Real(2), flow, crossings, transitions, initial);
        for (auto kind: {CrossingKind::NEGATIVE, CrossingKind::POSITIVE, CrossingKind::DECREASING}) {
            TransitionData transition;
            transition.event = event;
            transition.event_kind = EventKind::URGENT;
            RealExpression guard = kind == CrossingKind::DECREASING ? -RealExpression(x)
                                                                    : RealExpression(kind == CrossingKind::POSITIVE ? 1 : -1);
            transition.guard_function = make_function(space, guard);
            auto jump = initial;
            probe._apply_guard_step(jump, dynamic, flow, timing, transition, CrossingData(kind), Semantics::UPPER);
            ARIADNE_TESTING_EXPECT_EQUAL(jump.number_of_parameters(), initial.number_of_parameters());
            ARIADNE_TESTING_EXPECT_EQUAL(definitely(jump.is_empty()), kind == CrossingKind::NEGATIVE);
        }
    }

    void test_invariant_rejects_initial_jump_and_empty_split_halves() {
        DiscreteLocation q(1), target(2);
        DiscreteEvent guard("jump"), invariant("invariant");
        RealVariable x("x"), y("y");
        RealSpace space({x, y});
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0});
        system.new_mode(target, {dot(x) = 0, dot(y) = 0});
        system.new_invariant(q, x <= -0.5_dec, invariant);
        system.new_transition(q, guard, target, {next(x) = x, next(y) = y}, x >= 0, EventKind::URGENT);
        StepProbe probe(system);
        probe.configuration().set_maximum_enclosure_radius(0.25_x);
        probe.configuration().set_enable_subdivisions(true);
        auto orbit = probe.orbit(HybridExactBox(q, space, ExactBoxType({
                                                              {-1, 1},
                                                              { 0, 0}
        })),
                                 HybridTime(0.03125_x, 1), Semantics::UPPER);
        ARIADNE_TESTING_PROVE(orbit.final().size() > 0u);
        for (const auto& set: orbit.reach())
            ARIADNE_TESTING_PROVE(set.location() == q);
        for (const auto& set: orbit.intermediate())
            ARIADNE_TESTING_PROVE(set.location() == q);
        auto excluded = probe.orbit(HybridExactBox(q, space, ExactBoxType({
                                                                 {1, 2},
                                                                 {0, 0}
        })),
                                    HybridTime(0.03125_x, 1), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(excluded.final().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(excluded.intermediate().size(), 0u);
    }

    void test_timing_straddles_final_time_without_unwinding() {
        DiscreteLocation q(1);
        RealVariable x("x"), y("y");
        RealSpace space({x, y});
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0});
        StepProbe probe(system);
        probe.ALLOW_UNWIND = true;
        auto initial = probe.enclosure(HybridExactBox(q, space, ExactBoxType({
                                                                    {-1, 1},
                                                                    {-1, 1}
        })));
        auto factory = probe.function_factory();
        initial.set_time_function(factory.create_coordinate(initial.parameter_domain(), 0u));
        auto dynamic = system.dynamic_function(q);
        FlowFunctionModel flow = probe._compute_integration_step(dynamic, cast_exact_box(initial.state_bounding_box()), Dyadic(0.5_x)).reach_models()[0];
        Set<DiscreteEvent> active;
        Map<DiscreteEvent, CrossingData> crossings;
        Map<DiscreteEvent, TransitionData> transitions;
        auto timing = probe._estimate_timing(active, Real(1 / 8_q), flow, crossings, transitions, initial);
        ARIADNE_TESTING_EXPECT_EQUAL(timing.step_kind, StepKind::CONSTANT_EVOLUTION_TIME);
        ARIADNE_TESTING_EXPECT_EQUAL(timing.finishing_kind, FinishingKind::STRADDLE_FINAL_TIME);
        ARIADNE_TESTING_PROVE(definitely(timing.parameter_dependent_evolution_time.range().lower_bound() > 0));

        // A guard can be active even though the complete set finishes in
        // this step. Creep must preserve the exact final-time termination.
        initial.set_time_function(factory.create_constant(initial.parameter_domain(), Real(0)));
        DiscreteEvent event("crossing");
        CrossingData crossing(CrossingKind::TRANSVERSE);
        crossing.crossing_time = 0.125_x - factory.create_coordinate(flow.space_domain(), 0u);
        crossings.insert(event, crossing);
        TransitionData transition;
        transition.event = event;
        transition.event_kind = EventKind::URGENT;
        transition.guard_function = make_function(space, x - 0.125_dec);
        transitions.insert(event, transition);
        active.insert(event);
        auto finishing = probe._estimate_timing(active, Real(1 / 4_q), flow, crossings, transitions, initial);
        ARIADNE_TESTING_EXPECT_EQUAL(finishing.finishing_kind, FinishingKind::AT_FINAL_TIME);
        ARIADNE_TESTING_EXPECT_EQUAL(finishing.step_kind, StepKind::CONSTANT_FINISHING_TIME);
        initial.set_time_function(factory.create_coordinate(initial.parameter_domain(), 0u));
        auto straddling_crossing = probe._estimate_timing(active, Real(1 / 8_q), flow, crossings, transitions, initial);
        ARIADNE_TESTING_EXPECT_EQUAL(straddling_crossing.finishing_kind, FinishingKind::STRADDLE_FINAL_TIME);

        // Conservative active-event detection may retain a future crossing.
        // It must not shorten a step toward a guard that is not yet reached.
        initial.set_time_function(factory.create_constant(initial.parameter_domain(), Real(0)));
        crossings[event].crossing_time = 3 - factory.create_coordinate(flow.space_domain(), 0u);
        transitions[event].guard_function = make_function(space, x - 3);
        auto future = probe._estimate_timing(active, Real(2), flow, crossings, transitions, initial);
        ARIADNE_TESTING_EXPECT_EQUAL(future.step_kind, StepKind::CONSTANT_EVOLUTION_TIME);
        ARIADNE_TESTING_EXPECT_EQUAL(future.finishing_kind, FinishingKind::BEFORE_FINAL_TIME);
    }

    void test_rejected_step_and_guard_inputs() {
        LogVerbosityScope logging(100u);
        DiscreteLocation q(1);
        DiscreteEvent event("guard");
        RealVariable x("x"), y("y");
        RealSpace space({x, y});
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0});
        StepProbe probe(system);
        auto initial = probe.enclosure(HybridExactBox(q, space, ExactBoxType({
                                                                    {-1, 1},
                                                                    {-1, 1}
        })));
        auto dynamic = system.dynamic_function(q);
        FlowFunctionModel flow = probe._compute_integration_step(dynamic, cast_exact_box(initial.state_bounding_box()), Dyadic(0.5_x)).reach_models()[0];
        Set<DiscreteEvent> active;
        Map<DiscreteEvent, CrossingData> crossings;
        Map<DiscreteEvent, TransitionData> transitions;
        auto timing = probe._estimate_timing(active, Real(2), flow, crossings, transitions, initial);
        auto invalid_timing = timing;
        invalid_timing.step_kind = static_cast<StepKind>(255u);
        auto target = initial;
        ARIADNE_TESTING_EXPECT_EXCEPTION(probe._apply_evolve_step(target, flow, invalid_timing));
        ARIADNE_TESTING_EXPECT_EQUAL(target.number_of_parameters(), initial.number_of_parameters());

        TransitionData transition;
        transition.event = event;
        transition.event_kind = static_cast<EventKind>(255u);
        transition.guard_function = make_function(space, x);
        CrossingData crossing(CrossingKind::DEGENERATE);
        ARIADNE_TESTING_EXPECT_EXCEPTION(probe._apply_guard_step(
            target, dynamic, flow, timing, transition, crossing, Semantics::UPPER));
        transition.event_kind = EventKind::URGENT;
        crossing.crossing_kind = static_cast<CrossingKind>(255u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(probe._apply_guard_step(
            target, dynamic, flow, timing, transition, crossing, Semantics::UPPER));

        auto elapsed = probe.function_factory().create_constant(initial.parameter_domain(), 0.5_x);
        List<HybridEnclosure> sets{initial};
        ARIADNE_TESTING_EXPECT_EXCEPTION(probe._apply_guard(
            sets, elapsed, initial, flow, transition, crossing, Semantics::UPPER));
        crossing.crossing_kind = CrossingKind::DEGENERATE;
        ARIADNE_TESTING_EXPECT_EXCEPTION(probe._apply_guard(
            sets, elapsed, initial, flow, transition, crossing, static_cast<Semantics>(255u)));
        for (auto kind: {CrossingKind::NEGATIVE, CrossingKind::POSITIVE, CrossingKind::DECREASING}) {
            crossing.crossing_kind = kind;
            auto untouched = sets;
            probe._apply_guard(untouched, elapsed, initial, flow, transition, crossing, Semantics::UPPER);
            ARIADNE_TESTING_EXPECT_EQUAL(untouched[0].number_of_constraints(), initial.number_of_constraints());
        }
        auto wrong_elapsed = probe.function_factory().create_constant(ExactBoxType({
                                                                          {-1, 1}
        }),
                                                                      0.5_x);
        ARIADNE_TESTING_EXPECT_EXCEPTION(probe._apply_guard(
            sets, wrong_elapsed, initial, flow, transition, crossing, Semantics::UPPER));
    }

    void test_grazing_whose_maximum_is_in_the_past() {
        LogVerbosityScope logging(100u);
        DiscreteLocation q(1);
        DiscreteEvent event("past maximum");
        RealVariable x("x"), y("y");
        RealSpace space({x, y});
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0});
        StepProbe probe(system);
        auto initial = probe.enclosure(HybridExactBox(q, space,
                                                      ExactBoxType({
                                                          {1.0_x, 1.25_x},
                                                          {0.5_x, 0.75_x}
        })));
        auto factory = probe.function_factory();
        auto dynamic = system.dynamic_function(q);
        FlowFunctionModel flow = probe._compute_integration_step(dynamic, cast_exact_box(initial.state_bounding_box()), Dyadic(0.5_x)).reach_models()[0];
        auto elapsed = factory.create_constant(initial.parameter_domain(), 0.5_x);
        auto final_set = initial;
        final_set.apply_parameter_evolve_step(flow, elapsed);
        List<HybridEnclosure> sets{final_set};
        TransitionData transition;
        transition.event = event;
        transition.event_kind = EventKind::INVARIANT;
        transition.guard_function = make_function(space, y - sqr(x));
        CrossingData crossing(CrossingKind::GRAZING);
        crossing.critical_time = -factory.create_coordinate(flow.space_domain(), 0u);
        crossing.guard_range_at_critical_time = factory.create_coordinate(flow.space_domain(), 1u).range();
        probe._apply_guard(sets, elapsed, initial, flow, transition, crossing, Semantics::UPPER);
        // x only increases from x>=1, so y-x^2 stays strictly negative.
        ARIADNE_TESTING_EXPECT_EQUAL(sets.size(), 1u);
        ARIADNE_TESTING_PROVE(!definitely(sets[0].is_empty()));
        ARIADNE_TESTING_PROVE(definitely(sets[0].inside(HybridExactBox(q, space,
                                                                       ExactBoxType({
                                                                           {1.484375_x, 1.765625_x},
                                                                           {0.484375_x, 0.765625_x}
        })))));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_step_contracts");
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_full_steps_use_independent_endpoint_and_keep_hybrid_metadata);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_shortened_steps_use_the_tube_at_the_stopping_time);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_constrained_source_stays_inside_the_certified_composition_domain);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_unsupported_step_shapes_are_rejected_before_publication);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_constant_and_decreasing_guard_steps);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_invariant_rejects_initial_jump_and_empty_split_halves);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_timing_straddles_final_time_without_unwinding);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_rejected_step_and_guard_inputs);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_grazing_whose_maximum_is_in_the_past);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
