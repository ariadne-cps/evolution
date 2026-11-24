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

#include "../src/hybrid_guard_predicate.hpp"
#include "../src/hybrid_step_results.hpp"
#include "../src/hybrid_transition_handlers.hpp"
#include "evolution/hybrid_continuous_time_evolver.hpp"
#include "integration/integrator.hpp"
#include "integration/validated_integration_step.hpp"
#include "log_verbosity_scope.hpp"
#include "logging/test_runner.hpp"
#include "system/hybrid_automata.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;
using namespace Ariadne::EvolutionDetail;

class StepProbe : public HybridContinuousTimeEvolver<ValidatedTag> {
  public:
    using HybridContinuousTimeEvolver<ValidatedTag>::_estimate_timing;
    using HybridContinuousTimeEvolver<ValidatedTag>::HybridContinuousTimeEvolver;
    using HybridContinuousTimeEvolverBase<ValidatedTag>::_apply_evolve_step;
    using HybridContinuousTimeEvolverBase<ValidatedTag>::_compute_integration_step;
};

// An exact analytic flow x(t)=x0+t, y(t)=y0 isolates guard mathematics
// from adaptive integration and makes elapsed-time witnesses explicit.
struct LinearGuardScenario {
    const RealVariable x{"x"}, y{"y"};
    const RealSpace space{x, y};
    TaylorPicardIntegrator integrator{StepMaximumError(1e-6)};
    ValidatedFunctionPatchFactory factory{integrator.function_factory()};
    HybridEnclosure initial;
    FlowFunctionModel flow;
    EffectiveVectorMultivariateFunction dynamic;
    TimingData timing;

    static FlowFunctionModel make_flow(const ValidatedFunctionPatchFactory& patch_factory,
                                       const ExactBoxType& box, const Dyadic& h) {
        const auto domain = product(box, ExactIntervalType(0u, h));
        auto model = patch_factory.create_zeros(2u, domain);
        model[0u] = patch_factory.create_coordinate(domain, 0u) + patch_factory.create_coordinate(domain, 2u);
        model[1u] = patch_factory.create_coordinate(domain, 1u);
        return model;
    }

    LinearGuardScenario(const ExactBoxType& box, const Dyadic& h)
        : initial(DiscreteLocation(1), space, Enclosure(box, EnclosureConfiguration(factory))),
          flow(make_flow(factory, box, h)), dynamic(make_function(space, Vector<RealExpression>{1, 0})) {
        initial.set_time_function(factory.create_constant(initial.parameter_domain(), Real(2)));
        timing.step_size = h;
        timing.step_kind = StepKind::CONSTANT_EVOLUTION_TIME;
        timing.finishing_kind = FinishingKind::BEFORE_FINAL_TIME;
        timing.final_time = Real(4);
        timing.evolution_time_domain = ExactIntervalType(0u, h);
        timing.evolution_time_coordinate = factory.create_coordinate({timing.evolution_time_domain}, 0u);
        timing.parameter_dependent_evolution_time = factory.create_constant(initial.parameter_domain(), Real(h));
    }

    TransitionData transition(EventKind kind, const RealExpression& guard) const {
        TransitionData result;
        result.event = DiscreteEvent("guard");
        result.event_kind = kind;
        result.guard_function = make_function(space, guard);
        return result;
    }

    HybridGuardContext context(const TransitionData& transition, const CrossingData& crossing, Semantics semantics) const {
        return {initial, flow, transition, crossing, semantics, factory};
    }
};

// Both certificates enclose x'=1, but the tube deliberately carries a wider
// uniform error. Its endpoint slice must not replace the independent endpoint.
class CertificateProbe : public HybridContinuousTimeEvolver<ValidatedTag> {
  public:
    enum class Shape { STRUCTURED,
                       MULTIPATCH,
                       BOX_ENDPOINT,
                       INPUT_COORDINATE };
    Shape shape = Shape::STRUCTURED;
    mutable std::optional<ExactBoxType> certified_domain;
    using HybridContinuousTimeEvolver<ValidatedTag>::HybridContinuousTimeEvolver;
    using HybridContinuousTimeEvolverBase<ValidatedTag>::_evolution_step;

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

    void test_transverse_and_permissive_handlers_keep_distinct_time_maps() {
        LinearGuardScenario scenario(ExactBoxType({
                                         {-1, -0.5_x},
                                         { 0,      0}
        }),
                                     Dyadic(1));
        CrossingData crossing(CrossingGeometry::INCREASING,
                              CrossingTimeCertificate(-scenario.factory.create_coordinate(scenario.initial.parameter_domain(), 0u)));
        const auto first_entry_flow = flow_guard_capability(crossing, Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(first_entry_flow.method, HybridFlowGuardMethod::FIRST_ENTRY_BOUND);
        ARIADNE_TESTING_EXPECT_EQUAL(first_entry_flow.evidence, CrossingEvidenceKind::CROSSING_TIME);
        ARIADNE_TESTING_PROVE(first_entry_flow.uses_timed_flow);
        const auto first_entry_jump = jump_guard_capability(EventKind::URGENT, crossing, Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(first_entry_jump.method, HybridJumpGuardMethod::FIRST_ENTRY);
        ARIADNE_TESTING_EXPECT_EQUAL(first_entry_jump.evidence, CrossingEvidenceKind::CROSSING_TIME);
        ARIADNE_TESTING_PROVE(!first_entry_jump.requires_lie_derivative);
        for (auto semantics: {Semantics::UPPER, Semantics::LOWER}) {
            for (auto kind: {EventKind::URGENT, EventKind::IMPACT}) {
                const auto event = scenario.transition(kind, scenario.x);
                const auto candidate = make_jump_candidate(scenario.context(event, crossing, semantics),
                                                           scenario.dynamic, scenario.timing);
                ARIADNE_TESTING_EXPECT_EQUAL(candidate.elapsed_time.domain(), candidate.enclosure.parameter_domain());
                ARIADNE_TESTING_EXPECT_EQUAL(candidate.enclosure.number_of_parameters(), 2u);
                ARIADNE_TESTING_PROVE(subset(cast_exact_box(candidate.enclosure.state_bounding_box()),
                                             ExactBoxType({
                                                 {-0x1p-40_x, 0x1p-40_x},
                                                 {         0,         0}
                })));
                ARIADNE_TESTING_PROVE(subset(ExactIntervalType(2.5_x, 3.0_x),
                                             cast_exact_interval(candidate.enclosure.time_range())));
                const auto clock_identity = candidate.enclosure.time_function() - candidate.elapsed_time - Real(2);
                ARIADNE_TESTING_PROVE(subset(cast_exact_interval(clock_identity.range()),
                                             ExactIntervalType(-0x1p-40_x, 0x1p-40_x)));
            }
            const auto event = scenario.transition(EventKind::PERMISSIVE, scenario.x);
            const auto candidate = make_jump_candidate(scenario.context(event, crossing, semantics),
                                                       scenario.dynamic, scenario.timing);
            ARIADNE_TESTING_EXPECT_EQUAL(candidate.enclosure.number_of_parameters(), 3u);
            const auto optional = jump_guard_capability(EventKind::PERMISSIVE, crossing, semantics);
            ARIADNE_TESTING_EXPECT_EQUAL(optional.method, HybridJumpGuardMethod::TUBE_ACTIVATION);
            ARIADNE_TESTING_EXPECT_EQUAL(optional.evidence, CrossingEvidenceKind::NONE);
            ARIADNE_TESTING_PROVE(optional.uses_timed_flow);
            ARIADNE_TESTING_PROVE(!optional.requires_lie_derivative);
            ARIADNE_TESTING_EXPECT_EQUAL(candidate.elapsed_time.domain(), candidate.enclosure.parameter_domain());
            // Optional activation includes the interior after entry, not just x=0.
            const auto box = cast_exact_box(candidate.enclosure.state_bounding_box());
            ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                             {0, 0.5_x},
                                             {0,     0}
            }),
                                         box));
            ARIADNE_TESTING_PROVE(box[0u].upper_bound() <= FloatDP(0.50000095367431640625_x, dp));
            const Vector<FloatDPBounds> admitted{FloatDPBounds(-0.75_x, dp), FloatDPBounds(0, dp), FloatDPBounds(1, dp)};
            const Vector<FloatDPBounds> inactive{FloatDPBounds(-1, dp), FloatDPBounds(0, dp), FloatDPBounds(0.5_x, dp)};
            bool rejects_inactive = false;
            for (const auto& constraint: candidate.enclosure.constraints()) {
                const FloatDPBounds lower(constraint.lower_bound(), dp);
                const FloatDPBounds upper(constraint.upper_bound(), dp);
                const auto value = constraint.function()(admitted);
                ARIADNE_TESTING_PROVE(value.lower_raw() >= lower.upper_raw());
                ARIADNE_TESTING_PROVE(value.upper_raw() <= upper.lower_raw());
                const auto rejected = constraint.function()(inactive);
                if (rejected.upper_raw() < lower.lower_raw() || rejected.lower_raw() > upper.upper_raw())
                    rejects_inactive = true;
            }
            ARIADNE_TESTING_PROVE(rejects_inactive);
            const auto state = candidate.enclosure.state_function()(admitted);
            ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.25_x, dp), state[0u]));
        }
    }

    void test_impact_direction_and_source_candidate_elapsed_time() {
        LinearGuardScenario scenario(ExactBoxType({
                                         {-1, -0.5_x},
                                         { 0,      0}
        }),
                                     Dyadic(1));
        const CrossingData active(CrossingGeometry::POSITIVE);
        auto event = scenario.transition(EventKind::URGENT, -scenario.x);
        const auto candidate = make_jump_candidate(scenario.context(event, active, Semantics::UPPER),
                                                   scenario.dynamic, scenario.timing);
        ARIADNE_TESTING_PROVE(!definitely(candidate.enclosure.is_empty()));
        ARIADNE_TESTING_EXPECT_EQUAL(candidate.elapsed_time.domain(), candidate.enclosure.parameter_domain());
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(candidate.elapsed_time.range()), ExactIntervalType(0, 0)));
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(candidate.enclosure.time_range()), ExactIntervalType(2, 2)));
        // The same active source cannot impact against a strictly negative direction.
        event.event_kind = EventKind::IMPACT;
        const auto source_entry = jump_guard_capability(EventKind::URGENT, active, Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(source_entry.method, HybridJumpGuardMethod::SOURCE_ACTIVATION);
        ARIADNE_TESTING_PROVE(!source_entry.requires_lie_derivative);
        const auto impact_entry = jump_guard_capability(EventKind::IMPACT, active, Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(impact_entry.method, HybridJumpGuardMethod::SOURCE_ACTIVATION);
        ARIADNE_TESTING_PROVE(impact_entry.requires_lie_derivative);
        const auto impact = make_jump_candidate(scenario.context(event, active, Semantics::UPPER),
                                                scenario.dynamic, scenario.timing);
        ARIADNE_TESTING_PROVE(definitely(impact.enclosure.is_empty()));
        // A competing blocking guard uses the returned zero-time map directly.
        const auto other = scenario.transition(EventKind::INVARIANT, scenario.x - 2);
        const CrossingData later(CrossingGeometry::INCREASING,
                                 CrossingTimeCertificate(scenario.factory.create_constant(scenario.initial.parameter_domain(), 0.25_x)));
        List<HybridEnclosure> guarded{candidate.enclosure};
        apply_flow_guard(guarded, scenario.context(other, later, Semantics::UPPER), candidate.elapsed_time);
        ARIADNE_TESTING_EXPECT_EQUAL(guarded.size(), 1u);
        ARIADNE_TESTING_PROVE(!definitely(guarded.front().is_empty()));
    }

    void test_initial_guard_activation_classifies_the_complete_source() {
        // On x in [-1/4,1/4], the three guards respectively exclude
        // every source point, enable a strict subset, or enable all points.
        // This is a single-source contract, not collective tree certification.
        LinearGuardScenario scenario(ExactBoxType({
                                         {-0.25_x, 0.25_x},
                                         {      0,      0}
        }),
                                     Dyadic(0.25_x));
        const auto excluded_guard = scenario.transition(EventKind::URGENT, scenario.x - 1);
        const auto partial_guard = scenario.transition(EventKind::URGENT, scenario.x);
        const auto full_guard = scenario.transition(EventKind::URGENT, scenario.x + 1);
        const auto excluded = classify_initial_guard(scenario.initial, excluded_guard);
        const auto partial = classify_initial_guard(scenario.initial, partial_guard);
        const auto full = classify_initial_guard(scenario.initial, full_guard);

        ARIADNE_TESTING_EXPECT_EQUAL(excluded, HybridInitialGuardActivation::EXCLUDED);
        ARIADNE_TESTING_EXPECT_EQUAL(partial, HybridInitialGuardActivation::UNRESOLVED);
        ARIADNE_TESTING_EXPECT_EQUAL(full, HybridInitialGuardActivation::FULL_SOURCE);
        ARIADNE_TESTING_PROVE(!permits_initial_reset(excluded, Semantics::UPPER));
        ARIADNE_TESTING_PROVE(!permits_initial_reset(excluded, Semantics::LOWER));
        ARIADNE_TESTING_PROVE(permits_initial_reset(partial, Semantics::UPPER));
        ARIADNE_TESTING_PROVE(!permits_initial_reset(partial, Semantics::LOWER));
        ARIADNE_TESTING_PROVE(permits_initial_reset(full, Semantics::UPPER));
        ARIADNE_TESTING_PROVE(permits_initial_reset(full, Semantics::LOWER));
        ARIADNE_TESTING_EXPECT_EXCEPTION(
            permits_initial_reset(full, static_cast<Semantics>(255u)));
    }

    void test_shared_source_guard_matches_continuous_activation() {
        LinearGuardScenario scenario(ExactBoxType({
                                         {-0.25_x, 0.25_x},
                                         {      0,      0}
        }), Dyadic(0.25_x));
        const auto partial = scenario.transition(EventKind::URGENT, scenario.x);
        HybridEnclosure common = scenario.initial;
        HybridEnclosure original = scenario.initial;
        restrict_source_guard(common, partial.event, partial.guard_function,
                              partial.source_arguments);
        original.new_activation(partial.event, partial.guard_function,
                                partial.source_arguments);
        ARIADNE_TESTING_EXPECT_EQUAL(common.number_of_constraints(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(common.number_of_constraints(), original.number_of_constraints());
        ARIADNE_TESTING_EXPECT_EQUAL(scenario.initial.number_of_constraints(), 0u);
        ARIADNE_TESTING_PROVE(!definitely(common.is_empty()));
        ARIADNE_TESTING_EXPECT_EQUAL(classify_source_guard(
            scenario.initial, partial.guard_function, partial.source_arguments),
            HybridSourceGuardStatus::UNRESOLVED);
        const auto excluded = scenario.transition(EventKind::URGENT, scenario.x - 1);
        const auto full = scenario.transition(EventKind::URGENT, scenario.x + 1);
        ARIADNE_TESTING_EXPECT_EQUAL(classify_source_guard(
            scenario.initial, excluded.guard_function, excluded.source_arguments),
            HybridSourceGuardStatus::EXCLUDED);
        ARIADNE_TESTING_EXPECT_EQUAL(classify_source_guard(
            scenario.initial, full.guard_function, full.source_arguments),
            HybridSourceGuardStatus::FULL_SOURCE);
        // The continuous hybrid LOWER policy still requires whole-source evidence.
        ARIADNE_TESTING_PROVE(!permits_initial_reset(
            classify_initial_guard(scenario.initial, partial), Semantics::LOWER));
        ARIADNE_TESTING_PROVE(permits_initial_reset(
            classify_initial_guard(scenario.initial, partial), Semantics::UPPER));
    }

    void test_unresolved_geometry_preserves_upper_lower_rules() {
        LinearGuardScenario scenario(ExactBoxType({
                                         {-1, -0.5_x},
                                         { 0,      0}
        }),
                                     Dyadic(1));
        const auto event = scenario.transition(EventKind::URGENT, scenario.x);
        const CrossingData unresolved(CrossingGeometry::UNRESOLVED);
        const auto unresolved_outer = flow_guard_capability(unresolved, Semantics::UPPER);
        const auto unresolved_inner = flow_guard_capability(unresolved, Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(unresolved_outer.method, HybridFlowGuardMethod::TWO_SAMPLES);
        ARIADNE_TESTING_EXPECT_EQUAL(unresolved_inner.method, HybridFlowGuardMethod::DISCARD);
        ARIADNE_TESTING_EXPECT_EQUAL(unresolved_outer.evidence, CrossingEvidenceKind::NONE);
        const auto unresolved_jump = jump_guard_capability(EventKind::URGENT, unresolved, Semantics::UPPER);
        const auto unresolved_discard = jump_guard_capability(EventKind::URGENT, unresolved, Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(unresolved_jump.method, HybridJumpGuardMethod::DIRECTED_BOUNDARY);
        ARIADNE_TESTING_PROVE(unresolved_jump.requires_lie_derivative);
        ARIADNE_TESTING_EXPECT_EQUAL(unresolved_discard.method, HybridJumpGuardMethod::DISCARD);
        const auto upper = make_jump_candidate(scenario.context(event, unresolved, Semantics::UPPER),
                                               scenario.dynamic, scenario.timing);
        const auto lower = make_jump_candidate(scenario.context(event, unresolved, Semantics::LOWER),
                                               scenario.dynamic, scenario.timing);
        ARIADNE_TESTING_PROVE(!definitely(upper.enclosure.is_empty()));
        ARIADNE_TESTING_PROVE(definitely(lower.enclosure.is_empty()));
        auto endpoint = scenario.initial;
        endpoint.apply_parameter_evolve_step(scenario.flow, scenario.timing.parameter_dependent_evolution_time);
        for (auto kind: {CrossingGeometry::UNRESOLVED, CrossingGeometry::CONCAVE}) {
            const CrossingData uncertain(kind);
            List<HybridEnclosure> outer{endpoint}, inner{endpoint};
            apply_flow_guard(outer, scenario.context(event, uncertain, Semantics::UPPER),
                             scenario.timing.parameter_dependent_evolution_time);
            apply_flow_guard(inner, scenario.context(event, uncertain, Semantics::LOWER),
                             scenario.timing.parameter_dependent_evolution_time);
            ARIADNE_TESTING_PROVE(!definitely(outer.front().is_empty()));
            ARIADNE_TESTING_PROVE(definitely(inner.front().is_empty()));
            ARIADNE_TESTING_PROVE(definitely(outer.front().inside(HybridExactBox(
                scenario.initial.location(), scenario.space,
                ExactBoxType({
                    {-0x1p-20_x, 0x1p-20_x},
                    {-0x1p-40_x, 0x1p-40_x}
            })))));
        }
    }

    void test_grazing_split_keeps_primary_and_additional_branch_order() {
        LinearGuardScenario scenario(ExactBoxType({
                                         {      -1,       1},
                                         {-0.125_x, 0.125_x}
        }),
                                     Dyadic(0.5_x));
        const auto event = scenario.transition(EventKind::INVARIANT, scenario.y - sqr(scenario.x));
        scenario.initial.new_invariant(event.event, event.guard_function);
        const auto critical = -scenario.factory.create_coordinate(scenario.initial.parameter_domain(), 0u);
        auto critical_bounds = scenario.initial.state_bounding_box();
        critical_bounds[0u] = UpperIntervalType(0, 0);
        const CrossingData crossing(GuardMaximumCertificate(critical,
                                                            scenario.factory.create_coordinate(scenario.initial.parameter_domain(), 1u).range(), critical_bounds));
        const auto maximum_flow = flow_guard_capability(crossing, Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(maximum_flow.method, HybridFlowGuardMethod::MAXIMUM_SPLIT);
        ARIADNE_TESTING_EXPECT_EQUAL(maximum_flow.evidence, CrossingEvidenceKind::GUARD_MAXIMUM);
        ARIADNE_TESTING_PROVE(maximum_flow.uses_timed_flow);
        const auto maximum_jump = jump_guard_capability(EventKind::URGENT, crossing, Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(maximum_jump.method, HybridJumpGuardMethod::DIRECTED_BOUNDARY);
        ARIADNE_TESTING_EXPECT_EQUAL(maximum_jump.evidence, CrossingEvidenceKind::GUARD_MAXIMUM);
        ARIADNE_TESTING_PROVE(maximum_jump.requires_lie_derivative);
        auto endpoint = scenario.initial;
        endpoint.apply_parameter_evolve_step(scenario.flow, scenario.timing.parameter_dependent_evolution_time);
        List<HybridEnclosure> branches{endpoint, endpoint};
        apply_flow_guard(branches, scenario.context(event, crossing, Semantics::UPPER),
                         scenario.timing.parameter_dependent_evolution_time);
        ARIADNE_TESTING_EXPECT_EQUAL(branches.size(), 6u);
        for (SizeType i = 0u; i != branches.size(); ++i) {
            ARIADNE_TESTING_PROVE(!definitely(branches[i].is_empty()));
            ExactBoxType bound(2u);
            Vector<FloatDPBounds> witness(2u, dp);
            if (i < 2u) {
                bound = ExactBoxType({
                    {-0.50000095367431640625_x, 1.50000095367431640625_x},
                    {-0.12500095367431640625_x,                0x1p-20_x}
                });
                witness = Vector<FloatDPBounds>{FloatDPBounds(-0.75_x, dp), FloatDPBounds(-0.03125_x, dp)};
            } else if (i % 2u == 0u) {
                bound = ExactBoxType({
                    {-0.50000095367431640625_x,                0x1p-20_x},
                    {               -0x1p-20_x, 0.12500095367431640625_x}
                });
                witness = Vector<FloatDPBounds>{FloatDPBounds(-0.75_x, dp), FloatDPBounds(0.03125_x, dp)};
            } else {
                bound = ExactBoxType({
                    {0.49999904632568359375_x, 1.50000095367431640625_x},
                    {              -0x1p-20_x, 0.12500095367431640625_x}
                });
                witness = Vector<FloatDPBounds>{FloatDPBounds(0.75_x, dp), FloatDPBounds(0.03125_x, dp)};
            }
            ARIADNE_TESTING_PROVE(definitely(branches[i].inside(HybridExactBox(
                scenario.initial.location(), scenario.space, bound))));
            // Each ordered branch has a concrete feasible analytic trajectory.
            for (const auto& constraint: branches[i].constraints()) {
                const auto value = constraint.function()(witness);
                const FloatDPBounds lower(constraint.lower_bound(), dp);
                const FloatDPBounds upper(constraint.upper_bound(), dp);
                ARIADNE_TESTING_PROVE(value.lower_raw() >= lower.upper_raw());
                ARIADNE_TESTING_PROVE(value.upper_raw() <= upper.lower_raw());
            }
            const auto state = branches[i].state_function()(witness);
            ARIADNE_TESTING_PROVE(refines(witness[0u] + FloatDPBounds(0.5_x, dp), state[0u]));
            ARIADNE_TESTING_PROVE(refines(witness[1u], state[1u]));
        }
    }

    void test_full_steps_use_independent_endpoint_and_keep_hybrid_metadata() {
        const DiscreteLocation q(1);
        const RealVariable x("x"), output("output");
        const RealSpace space({x});
        ContinuousTimeHybridAutomaton system;
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
        ContinuousTimeHybridAutomaton system;
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
        ContinuousTimeHybridAutomaton system;
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
        ContinuousTimeHybridAutomaton system;
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
        ContinuousTimeHybridAutomaton system;
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
        for (auto kind: {CrossingGeometry::NEGATIVE, CrossingGeometry::POSITIVE, CrossingGeometry::DECREASING}) {
            TransitionData transition;
            transition.event = event;
            transition.event_kind = EventKind::URGENT;
            RealExpression guard = kind == CrossingGeometry::DECREASING ? -RealExpression(x)
                                                                        : RealExpression(kind == CrossingGeometry::POSITIVE ? 1 : -1);
            transition.guard_function = make_function(space, guard);
            const auto candidate = make_jump_candidate(
                {initial, flow, transition, CrossingData(kind), Semantics::UPPER, probe.function_factory()}, dynamic, timing);
            const auto& jump = candidate.enclosure;
            ARIADNE_TESTING_EXPECT_EQUAL(jump.number_of_parameters(), initial.number_of_parameters());
            ARIADNE_TESTING_EXPECT_EQUAL(definitely(jump.is_empty()), kind == CrossingGeometry::NEGATIVE);
        }
    }

    void test_invariant_rejects_initial_jump_and_empty_split_halves() {
        DiscreteLocation q(1), target(2);
        DiscreteEvent guard("jump"), invariant("invariant");
        RealVariable x("x"), y("y");
        RealSpace space({x, y});
        ContinuousTimeHybridAutomaton system;
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
        ContinuousTimeHybridAutomaton system;
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
        const CrossingData crossing(CrossingGeometry::INCREASING, CrossingTimeCertificate(
                                                                      0.125_x - factory.create_coordinate(flow.space_domain(), 0u)));
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
        crossings[event] = CrossingData(CrossingGeometry::INCREASING, CrossingTimeCertificate(
                                                                          3 - factory.create_coordinate(flow.space_domain(), 0u)));
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
        ContinuousTimeHybridAutomaton system;
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
        CrossingData crossing(CrossingGeometry::UNRESOLVED);
        ARIADNE_TESTING_EXPECT_EXCEPTION(make_jump_candidate(
            {target, flow, transition, crossing, Semantics::UPPER, probe.function_factory()}, dynamic, timing));
        transition.event_kind = EventKind::URGENT;
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      CrossingData(static_cast<CrossingGeometry>(255u)));
        ARIADNE_TESTING_EXPECT_EXCEPTION(jump_guard_capability(EventKind::INVARIANT, crossing, Semantics::UPPER));
        ARIADNE_TESTING_EXPECT_EXCEPTION(jump_guard_capability(EventKind::PROGRESS, crossing, Semantics::UPPER));
        ARIADNE_TESTING_EXPECT_EXCEPTION(flow_guard_capability(crossing, static_cast<Semantics>(255u)));
        ARIADNE_TESTING_EXPECT_EXCEPTION(jump_guard_capability(EventKind::URGENT, crossing, static_cast<Semantics>(255u)));

        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      (void)CrossingTimeCertificate(ValidatedScalarMultivariateFunctionPatch{}));

        auto elapsed = probe.function_factory().create_constant(initial.parameter_domain(), 0.5_x);
        List<HybridEnclosure> sets{initial};
        crossing = CrossingData(CrossingGeometry::UNRESOLVED);
        ARIADNE_TESTING_EXPECT_EXCEPTION(apply_flow_guard(sets,
                                                          {initial, flow, transition, crossing, static_cast<Semantics>(255u), probe.function_factory()}, elapsed));
        for (auto kind: {CrossingGeometry::NEGATIVE, CrossingGeometry::POSITIVE, CrossingGeometry::DECREASING}) {
            crossing = CrossingData(kind);
            auto untouched = sets;
            apply_flow_guard(untouched,
                             {initial, flow, transition, crossing, Semantics::UPPER, probe.function_factory()}, elapsed);
            ARIADNE_TESTING_EXPECT_EQUAL(untouched[0].number_of_constraints(), initial.number_of_constraints());
        }
        auto wrong_elapsed = probe.function_factory().create_constant(ExactBoxType({
                                                                          {-1, 1}
        }),
                                                                      0.5_x);
        ARIADNE_TESTING_EXPECT_EXCEPTION(apply_flow_guard(sets,
                                                          {initial, flow, transition, crossing, Semantics::UPPER, probe.function_factory()}, wrong_elapsed));
    }

    void test_reported_endpoint_keeps_its_own_unbounded_continuation() {
        // For x'=1, x0 in [-1,-1/2], the first entry of x>=0
        // occurs at t=-x0 in [1/2,1]. At elapsed time one,
        // the reported endpoint must satisfy this first-entry constraint,
        // but the checkpoint for the next step remains unrestricted.
        LinearGuardScenario scenario(ExactBoxType({
                                         {-1, -0.5_x},
                                         { 0,      0}
        }),
                                     Dyadic(1));
        auto evolved = scenario.initial;
        evolved.apply_parameter_evolve_step(scenario.flow, scenario.timing.parameter_dependent_evolution_time);
        List<HybridEnclosure> endpoints{evolved};
        auto branches = pair_endpoint_branches(endpoints);
        ARIADNE_TESTING_EXPECT_EQUAL(branches.size(), 1u);
        const auto crossing_time = -scenario.factory.create_coordinate(scenario.flow.space_domain(), 0u);
        const CrossingData crossing(CrossingGeometry::INCREASING, CrossingTimeCertificate(crossing_time));
        const auto transition = scenario.transition(EventKind::URGENT, scenario.x);
        restrict_reported_endpoints(branches, scenario.context(transition, crossing, Semantics::UPPER),
                                    scenario.timing.parameter_dependent_evolution_time);
        ARIADNE_TESTING_EXPECT_EQUAL(branches.size(), 1u);
        ARIADNE_TESTING_PROVE(branches[0].reported.number_of_constraints() >
                              branches[0].continuation.number_of_constraints());
        ARIADNE_TESTING_EXPECT_EQUAL(branches[0].continuation.number_of_constraints(), evolved.number_of_constraints());

        List<HybridEnclosure> reach{scenario.initial};
        reach.front().apply_parameter_reach_step(scenario.flow, scenario.timing.parameter_dependent_evolution_time);
        HybridEvolutionStepResult result;
        auto published = publish_step_contributions(reach, branches, scenario.timing, result);
        ARIADNE_TESTING_PROVE(published.progress);
        ARIADNE_TESTING_PROVE(!published.finishing);
        ARIADNE_TESTING_EXPECT_EQUAL(result.successors.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.contributions.size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.contributions[0].kind, HybridEvolutionSetKind::INTERMEDIATE);
        ARIADNE_TESTING_EXPECT_EQUAL(result.contributions[1].kind, HybridEvolutionSetKind::REACH);
        ARIADNE_TESTING_EXPECT_EQUAL(result.successors[0].work.enclosure.number_of_constraints(), evolved.number_of_constraints());
        ARIADNE_TESTING_PROVE(result.contributions[0].enclosure.number_of_constraints() >
                              result.successors[0].work.enclosure.number_of_constraints());
    }

    void test_horizon_publication_keeps_closed_reach_endpoint() {
        LinearGuardScenario scenario(ExactBoxType({
                                         {-1, -0.5_x},
                                         { 0,      0}
        }),
                                     Dyadic(1));
        auto evolved = scenario.initial;
        evolved.apply_parameter_evolve_step(scenario.flow, scenario.timing.parameter_dependent_evolution_time);
        List<HybridEnclosure> endpoints{evolved};
        auto branches = pair_endpoint_branches(endpoints);
        List<HybridEnclosure> reach{scenario.initial};
        reach.front().apply_parameter_reach_step(scenario.flow, scenario.timing.parameter_dependent_evolution_time);

        // At t=5/2 the full-step endpoint is past the horizon, yet the
        // closed reach tube contains the terminal slice at exactly t=5/2.
        auto straddling = scenario.timing;
        straddling.final_time = Real(2.5_x);
        straddling.finishing_kind = FinishingKind::STRADDLE_FINAL_TIME;
        HybridEvolutionStepResult result;
        auto published = publish_step_contributions(reach, branches, straddling, result);
        ARIADNE_TESTING_PROVE(!published.progress);
        ARIADNE_TESTING_PROVE(published.finishing);
        ARIADNE_TESTING_EXPECT_EQUAL(result.successors.size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.contributions.size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.contributions[0].kind, HybridEvolutionSetKind::FINAL);
        ARIADNE_TESTING_EXPECT_EQUAL(result.contributions[1].kind, HybridEvolutionSetKind::REACH);
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(result.contributions[0].enclosure.time_range()),
                                     ExactIntervalType(2.5_x, 2.5_x)));
    }

    void test_grazing_whose_maximum_is_in_the_past() {
        LogVerbosityScope logging(100u);
        DiscreteLocation q(1);
        DiscreteEvent event("past maximum");
        RealVariable x("x"), y("y");
        RealSpace space({x, y});
        ContinuousTimeHybridAutomaton system;
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
        const auto critical = -factory.create_coordinate(flow.space_domain(), 0u);
        auto critical_bounds = initial.state_bounding_box();
        critical_bounds[0u] = UpperIntervalType(0, 0);
        const CrossingData crossing(GuardMaximumCertificate(critical,
                                                            factory.create_coordinate(flow.space_domain(), 1u).range(), critical_bounds));
        apply_flow_guard(sets,
                         {initial, flow, transition, crossing, Semantics::UPPER, probe.function_factory()}, elapsed);
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
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_transverse_and_permissive_handlers_keep_distinct_time_maps);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_impact_direction_and_source_candidate_elapsed_time);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_initial_guard_activation_classifies_the_complete_source);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_shared_source_guard_matches_continuous_activation);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_unresolved_geometry_preserves_upper_lower_rules);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_grazing_split_keeps_primary_and_additional_branch_order);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_full_steps_use_independent_endpoint_and_keep_hybrid_metadata);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_shortened_steps_use_the_tube_at_the_stopping_time);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_constrained_source_stays_inside_the_certified_composition_domain);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_unsupported_step_shapes_are_rejected_before_publication);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_constant_and_decreasing_guard_steps);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_invariant_rejects_initial_jump_and_empty_split_halves);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_timing_straddles_final_time_without_unwinding);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_rejected_step_and_guard_inputs);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_grazing_whose_maximum_is_in_the_past);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_reported_endpoint_keeps_its_own_unbounded_continuation);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts, test_horizon_publication_keeps_closed_reach_endpoint);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
