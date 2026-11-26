/***************************************************************************
 *            test_hybrid_inclusion_impacts.cpp
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

/** @file tests/test_hybrid_inclusion_impacts.cpp
 *  @brief Verify hybrid inclusion impacts behavior.
 */

#include "../src/hybrid_outer_flow.hpp"
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

namespace {
struct ContactScenario {
    const DiscreteLocation source{1}, target{2}, other{3};
    const DiscreteEvent impact{"contact"};
    const RealVariable x{"x"}, y{"y"}, u{"u"}, z{"z"}, output{"output"};
    TaylorPicardIntegrator base{StepMaximumError(1e-6)};
    ValidatedFunctionPatchFactory factory{base.function_factory()};

    ContinuousTimeHybridAutomaton system(const RealExpression& guard) const {
        ContinuousTimeHybridAutomaton result;
        result.new_mode(source, OrdinaryDifferentialAlgebraicInclusion(
                                    {dot(x) = 1, dot(y) = u}, RealVariablesBox({Dyadic(-0.125_x) <= u <= Dyadic(0.125_x)})));
        result.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(z) = 0}, {let(output) = 2 * z}));
        result.new_transition(source, impact, target, {next(z) = 7 + x}, guard >= 0, EventKind::IMPACT);
        return result;
    }

    // Exact witness u=0, x(t)=-1+(t-1), y(t)=0, absolute t in [1,3].
    // Contact with -x^2 occurs at absolute t=2 with zero guard derivative.
    HybridEnclosure exact_reach() const {
        const ExactBoxType box({
            {-1, -1},
            { 0,  0}
        });
        HybridEnclosure reach(source, RealSpace({x, y}), Enclosure(box, EnclosureConfiguration(factory)));
        reach.set_time_function(factory.create_constant(reach.parameter_domain(), Real(1)));
        const auto domain = product(box, ExactIntervalType(0, 2));
        auto flow = factory.create_zeros(2u, domain);
        flow[0u] = factory.create_coordinate(domain, 0u) + factory.create_coordinate(domain, 2u);
        flow[1u] = factory.create_coordinate(domain, 1u);
        reach.apply_full_reach_step(flow);
        return reach;
    }

    void configure(HybridContinuousTimeEvolver<ValidatedTag>& evolver, bool boxes, bool tree = false) const {
        if (boxes)
            evolver.set_integrator(FirstOrderMultifunctionIntegrator());
        else
            evolver.set_integrator(base, ZeroApproximation(), ConstantApproximation());
        evolver.configuration().set_maximum_step_size(0.03125);
        if (tree)
            evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
    }
};

bool excludes_witness(const HybridEnclosure& set, const Vector<FloatDPBounds>& point) {
    for (const auto& constraint: set.constraints()) {
        const auto value = constraint.function()(point);
        const FloatDPBounds lower(constraint.lower_bound(), dp), upper(constraint.upper_bound(), dp);
        if (value.upper_raw() < lower.lower_raw() || value.lower_raw() > upper.upper_raw())
            return true;
    }
    return false;
}
// A time-range box discards correlations between arrival and local time.
// Prove the horizon on the constrained enclosure as well: c<=u and
// c-time+T>=u imply time<=T, including on nonrectangular time domains.
bool certifies_upper_time_bound(const HybridEnclosure& set, const Real& horizon) {
    if (definitely(set.time_range().upper_bound() <= horizon))
        return true;
    for (const auto& constraint: set.constraints()) {
        const auto residual = cast_exact_interval(
            (constraint.function() - set.time_function() + horizon).range());
        const FloatDPBounds upper(constraint.upper_bound(), dp);
        if (residual.lower_bound() >= upper.upper_raw())
            return true;
    }
    return false;
}
} // namespace

class TestHybridInclusionImpacts : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    void test_exact_tangent_contact_retains_witness_and_reset_metadata() {
        ContactScenario scenario;
        const auto system = scenario.system(-sqr(scenario.x));
        const auto transitions = outer_flow_transitions(system, scenario.source);
        auto reach = scenario.exact_reach();
        apply_outer_flow_blockers(reach, transitions);
        HybridEvolutionStepResult result;
        publish_outer_jump_candidates(reach, transitions, Real(2), result);
        ARIADNE_TESTING_REQUIRE(result.successors.size() == 1u);
        const auto& successor = result.successors.front();
        const auto& reset = successor.work.enclosure;
        ARIADNE_TESTING_EXPECT_EQUAL(reset.location(), scenario.target);
        ARIADNE_TESTING_EXPECT_EQUAL(reset.state_space(), RealSpace({scenario.z}));
        ARIADNE_TESTING_EXPECT_EQUAL(reset.previous_events().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(reset.previous_events().back(), scenario.impact);
        ARIADNE_TESTING_EXPECT_EQUAL(*successor.event_kind, EventKind::IMPACT);
        ARIADNE_TESTING_PROVE(!successor.crossing.has_value());
        ARIADNE_TESTING_EXPECT_EQUAL(successor.resolution, HybridEvolutionTransitionResolution::UNRESOLVED);
        const Vector<FloatDPBounds> contact{FloatDPBounds(-1, dp), FloatDPBounds(0, dp), FloatDPBounds(1, dp)};
        ARIADNE_TESTING_PROVE(!excludes_witness(reset, contact));
        // Verify every retained constraint at the analytic tangency, rather
        // than inferring contact merely from an overlapping bounding box.
        for (const auto& constraint: reset.constraints()) {
            const auto value = constraint.function()(contact);
            const FloatDPBounds lower(constraint.lower_bound(), dp), upper(constraint.upper_bound(), dp);
            ARIADNE_TESTING_PROVE(value.lower_raw() >= lower.upper_raw());
            ARIADNE_TESTING_PROVE(value.upper_raw() <= upper.lower_raw());
        }
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(7, dp), reset.state_function()(contact)[0u]));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(14, dp), reset.auxiliary_function()(contact)[0u]));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(2, dp), reset.time_function()(contact)));
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(reset.dwell_time_range()), ExactIntervalType(0, 0)));
        ARIADNE_TESTING_PROVE(excludes_witness(reset,
                                               {FloatDPBounds(-1, dp), FloatDPBounds(0, dp), FloatDPBounds(0.5_x, dp)}));
        ARIADNE_TESTING_PROVE(excludes_witness(reset,
                                               {FloatDPBounds(-1, dp), FloatDPBounds(0, dp), FloatDPBounds(1.5_x, dp)}));
        // Boundary construction is independent of prior flow blockers.
        // A transverse guard has both signs on this raw reach, so activation
        // alone would incorrectly admit the positive-side witness too.
        const auto transverse_system = scenario.system(scenario.x);
        const auto transverse = outer_flow_transitions(transverse_system, scenario.source);
        HybridEvolutionStepResult boundary;
        publish_outer_jump_candidates(scenario.exact_reach(), transverse, Real(3), boundary);
        ARIADNE_TESTING_REQUIRE(boundary.successors.size() == 1u);
        const auto& boundary_set = boundary.successors.front().work.enclosure;
        ARIADNE_TESTING_PROVE(!excludes_witness(boundary_set, contact));
        ARIADNE_TESTING_PROVE(excludes_witness(boundary_set,
                                               {FloatDPBounds(-1, dp), FloatDPBounds(0, dp), FloatDPBounds(0.5_x, dp)}));
        ARIADNE_TESTING_PROVE(excludes_witness(boundary_set,
                                               {FloatDPBounds(-1, dp), FloatDPBounds(0, dp), FloatDPBounds(1.5_x, dp)}));
    }

    void test_closed_horizon_excludes_future_contact() {
        ContactScenario scenario;
        const auto system = scenario.system(-sqr(scenario.x));
        const auto transitions = outer_flow_transitions(system, scenario.source);
        auto reach = scenario.exact_reach();
        apply_outer_flow_blockers(reach, transitions);
        HybridEvolutionStepResult before, at_contact;
        publish_outer_jump_candidates(reach, transitions, Real(1.5_x), before);
        publish_outer_jump_candidates(reach, transitions, Real(2), at_contact);
        ARIADNE_TESTING_PROVE(before.successors.empty());
        ARIADNE_TESTING_EXPECT_EQUAL(at_contact.successors.size(), 1u);
    }

    void test_competing_blocker_excludes_contact_and_simultaneous_urgent_survives() {
        ContactScenario scenario;
        auto blocked = scenario.system(scenario.x);
        blocked.new_guard(scenario.source, DiscreteEvent("earlier"), scenario.x >= -0.25_x, EventKind::PROGRESS);
        const auto blockers = outer_flow_transitions(blocked, scenario.source);
        auto reach = scenario.exact_reach();
        apply_outer_flow_blockers(reach, blockers);
        HybridEvolutionStepResult excluded;
        publish_outer_jump_candidates(reach, blockers, Real(3), excluded);
        ARIADNE_TESTING_PROVE(excluded.successors.empty());

        auto simultaneous = scenario.system(scenario.x);
        simultaneous.new_mode(scenario.other, OrdinaryDifferentialAlgebraicEquation({dot(scenario.z) = 0}));
        const DiscreteEvent urgent("simultaneous");
        simultaneous.new_transition(scenario.source, urgent, scenario.other,
                                    {next(scenario.z) = 9}, scenario.x >= 0, EventKind::URGENT);
        const auto transitions = outer_flow_transitions(simultaneous, scenario.source);
        reach = scenario.exact_reach();
        apply_outer_flow_blockers(reach, transitions);
        HybridEvolutionStepResult result;
        publish_outer_jump_candidates(reach, transitions, Real(2), result);
        ARIADNE_TESTING_REQUIRE(result.successors.size() == 2u);
        bool contact_found = false, urgent_found = false;
        for (const auto& successor: result.successors) {
            contact_found = contact_found || successor.event == scenario.impact;
            urgent_found = urgent_found || successor.event == urgent;
            ARIADNE_TESTING_PROVE(!excludes_witness(successor.work.enclosure,
                                                    {FloatDPBounds(-1, dp), FloatDPBounds(0, dp), FloatDPBounds(1, dp)}));
            ARIADNE_TESTING_PROVE(!successor.crossing.has_value());
        }
        ARIADNE_TESTING_PROVE(contact_found && urgent_found);
    }

    void test_tangency_with_structured_and_box_methods_and_both_backends() {
        ContactScenario scenario;
        const auto system = scenario.system(scenario.y - sqr(scenario.x));
        for (bool boxes: {false, true}) {
            for (bool tree: {false, true}) {
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                scenario.configure(evolver, boxes, tree);
                const auto initial = evolver.enclosure(HybridExactBox(scenario.source,
                                                                      RealSpace({
                                                                          scenario.x, scenario.y
                }),
                                                                      ExactBoxType({{-0.0625_x, -0.0625_x}, {0, 0}})));
                const auto orbit = evolver.orbit(initial, HybridTime(0.125_x, 1), Semantics::UPPER);
                std::optional<ExactBoxType> target_states;
                std::optional<ExactIntervalType> target_times;
                for (const auto& final: orbit.final()) {
                    if (final.location() != scenario.target)
                        continue;
                    ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().back(), scenario.impact);
                    const auto box = cast_exact_box(final.state_bounding_box());
                    const auto time = cast_exact_interval(final.time_range());
                    target_states = target_states ? hull(*target_states, box) : box;
                    target_times = target_times ? hull(*target_times, time) : time;
                    ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                                          {6.875_x, 7.125_x}
                    })));
                    ARIADNE_TESTING_PROVE(time.lower_bound() >= FloatDP(0, dp));
                    ARIADNE_TESTING_PROVE(certifies_upper_time_bound(final, Real(0.125_x)));
                }
                ARIADNE_TESTING_REQUIRE(target_states.has_value() && target_times.has_value());
                // u=0 is a genuine tangent trajectory: x=0, y=0 at t=1/16.
                ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                                 {7, 7}
                }),
                                             *target_states));
                ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0.0625_x, 0.0625_x), *target_times));
                const auto terminal = evolver.orbit(initial, HybridTime(0u, 1), Semantics::UPPER);
                ARIADNE_TESTING_REQUIRE(terminal.final().size() > 0u);
                for (const auto& final: terminal.final()) {
                    ARIADNE_TESTING_EXPECT_EQUAL(final.location(), scenario.source);
                    ARIADNE_TESTING_PROVE(final.previous_events().empty());
                }
                // At t=0, g=y-x^2=-1/256 is strictly negative.
                // The structured method may retain the WHOLE initial
                // LOWER source. A box-only outer integrator still
                // cannot provide a correlated LOWER continuation.
                if (boxes) {
                    ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                                  evolver.orbit(initial, HybridTime(0u, 1), Semantics::LOWER));
                } else {
                    const auto lower_zero = evolver.orbit(
                        initial, HybridTime(0u, 1), Semantics::LOWER);
                    ARIADNE_TESTING_REQUIRE(lower_zero.final().size() == 1u);
                    const auto& lower_final = lower_zero.final().front();
                    ARIADNE_TESTING_PROVE(lower_final.location() == scenario.source);
                    ARIADNE_TESTING_PROVE(lower_final.previous_events().empty());
                    ARIADNE_TESTING_PROVE(subset(
                        ExactBoxType({
                            {-0.0625_x, -0.0625_x},
                            {        0,         0}
                    }),
                        cast_exact_box(lower_final.state_bounding_box())));
                    ARIADNE_TESTING_PROVE(subset(
                        cast_exact_interval(lower_final.time_range()),
                        ExactIntervalType(0, 0)));

                    // A tangent or input-dependent crossing is possible
                    // before t=1/8. The presently unblocked-flow-only
                    // LOWER path must not claim a legal full-source flow.
                    ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                                  evolver.orbit(initial, HybridTime(0.125_x, 1),
                                                                Semantics::LOWER));
                }
            }
        }
    }

    void test_impact_after_parameter_dependent_arrival() {
        ContactScenario scenario;
        auto system = scenario.system(scenario.y - sqr(scenario.x));
        const DiscreteLocation entry(4);
        const RealVariable s("s");
        const DiscreteEvent enter("enter inclusion");
        system.new_mode(entry, OrdinaryDifferentialAlgebraicEquation({dot(s) = 1}));
        system.new_transition(entry, enter, scenario.source,
                              {next(scenario.x) = -0.0625_x, next(scenario.y) = 0}, s >= 0, EventKind::URGENT);
        for (bool boxes: {false, true}) {
            for (bool tree: {false, true}) {
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                scenario.configure(evolver, boxes, tree);
                const auto orbit = evolver.orbit(HybridExactBox(entry, RealSpace({
                                                                           s
                }),
                                                                ExactBoxType({{-0.0625_x, -0.03125_x}})),
                                                 HybridTime(0.125_x, 2), Semantics::UPPER);
                std::optional<ExactBoxType> target_states;
                std::optional<ExactIntervalType> target_times;
                for (const auto& final: orbit.final()) {
                    if (final.location() != scenario.target)
                        continue;
                    ARIADNE_TESTING_REQUIRE(final.previous_events().size() == 2u);
                    ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events()[0u], enter);
                    ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events()[1u], scenario.impact);
                    const auto box = cast_exact_box(final.state_bounding_box());
                    const auto time = cast_exact_interval(final.time_range());
                    target_states = target_states ? hull(*target_states, box) : box;
                    target_times = target_times ? hull(*target_times, time) : time;
                    ARIADNE_TESTING_PROVE(time.lower_bound() >= FloatDP(0, dp));
                    ARIADNE_TESTING_PROVE(certifies_upper_time_bound(final, Real(0.125_x)));
                }
                ARIADNE_TESTING_REQUIRE(target_states.has_value() && target_times.has_value());
                ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                                 {7, 7}
                }),
                                             *target_states));
                // ODE entry takes [1/32,1/16], then u=0 touches after 1/16.
                ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0.09375_x, 0.125_x), *target_times));
            }
        }
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_inclusion_impacts");
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusionImpacts, test_exact_tangent_contact_retains_witness_and_reset_metadata);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusionImpacts, test_closed_horizon_excludes_future_contact);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusionImpacts, test_competing_blocker_excludes_contact_and_simultaneous_urgent_survives);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusionImpacts, test_tangency_with_structured_and_box_methods_and_both_backends);
    ARIADNE_TESTING_ADD_TEST(TestHybridInclusionImpacts, test_impact_after_parameter_dependent_arrival);
    return Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
