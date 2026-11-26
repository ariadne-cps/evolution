/***************************************************************************
 *            test_hybrid_discrete_time_step_planner.cpp
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

/** @file tests/test_hybrid_discrete_time_step_planner.cpp
 *  @brief Exhaustive structural cases for single-tick hybrid branching.
 */

#include <stdexcept>

#include "evolution/hybrid_discrete_time_step_planner.hpp"
#include "logging/test_runner.hpp"
#include "system/discrete_time_hybrid_automaton.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

class TestHybridDiscreteTimeStepPlanner : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    using Evidence = DiscreteTimePredicateEvidence;
    using Role = DiscreteTimePredicateRole;
    using Kind = HybridDiscreteTimeBranchKind;

    static DifferenceAlgebraicEquation dynamics() {
        const RealVariable x("x");
        return DifferenceAlgebraicEquation(List<PrimedRealAssignment>{next(x) = x + 1});
    }

    static HybridDiscreteTimeStepPlanner::EvidenceEvaluator constant(Evidence value) {
        return [value](DiscreteLocation, DiscreteEvent, Role, const ContinuousPredicate&) { return value; };
    }

    Void test_no_events() {
        DiscreteTimeHybridAutomaton automaton;
        const DiscreteLocation source(1);
        automaton.new_mode(source, dynamics());
        const auto plan = HybridDiscreteTimeStepPlanner::plan(automaton, source, Integer(4), constant(Evidence::UNKNOWN));
        ARIADNE_TESTING_EXPECT_EQUAL(plan.source_index, Integer(4));
        ARIADNE_TESTING_EXPECT_EQUAL(plan.successor_index, Integer(5));
        ARIADNE_TESTING_EXPECT_EQUAL(plan.branches.size(), 1u);
        ARIADNE_TESTING_PROVE(plan.branches[0u].kind == Kind::PROGRESSION);
        ARIADNE_TESTING_PROVE(plan.branches[0u].certain);
        ARIADNE_TESTING_PROVE(!plan.branches[0u].event.has_value());
        ARIADNE_TESTING_PROVE(plan.branches[0u].target == source);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      HybridDiscreteTimeStepPlanner::plan(automaton, DiscreteLocation(9), Integer(4), constant(Evidence::TRUE)));
    }

    Void test_urgent_complement() {
        DiscreteTimeHybridAutomaton automaton;
        const DiscreteLocation source(1), target(2);
        const DiscreteEvent event("urgent");
        const RealVariable x("x");
        automaton.new_mode(source, dynamics());
        automaton.new_mode(target, dynamics());
        automaton.new_urgent_transition(source, event, target, {next(x) = x}, x >= 0);

        for (const Evidence guard: {Evidence::FALSE, Evidence::TRUE, Evidence::UNKNOWN}) {
            SizeType evaluations = 0u;
            const auto eval = [&](DiscreteLocation, DiscreteEvent, Role role, const ContinuousPredicate&) {
                ARIADNE_TESTING_PROVE(role == Role::GUARD);
                ++evaluations;
                return guard;
            };
            const auto plan = HybridDiscreteTimeStepPlanner::plan(automaton, source, Integer(2), eval);
            ARIADNE_TESTING_EXPECT_EQUAL(evaluations, 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(plan.branches.size(), guard == Evidence::UNKNOWN ? 2u : 1u);
            if (guard != Evidence::FALSE) {
                ARIADNE_TESTING_PROVE(plan.branches[0u].kind == Kind::TRANSITION);
                ARIADNE_TESTING_PROVE(plan.branches[0u].event.has_value());
                ARIADNE_TESTING_PROVE(*plan.branches[0u].event == event);
                ARIADNE_TESTING_PROVE(plan.branches[0u].target == target);
                ARIADNE_TESTING_PROVE(plan.branches[0u].certain == (guard == Evidence::TRUE));
            }
            if (guard != Evidence::TRUE) {
                const auto& progress = plan.branches.back();
                ARIADNE_TESTING_PROVE(progress.kind == Kind::PROGRESSION);
                ARIADNE_TESTING_PROVE(progress.certain == (guard == Evidence::FALSE));
            }
            ARIADNE_TESTING_EXPECT_EQUAL(plan.successor_index, Integer(3));
        }
    }

    Void test_invalid_predicate_evidence() {
        DiscreteTimeHybridAutomaton automaton;
        const DiscreteLocation source(1), target(2);
        const DiscreteEvent event("urgent");
        const RealVariable x("x");
        automaton.new_mode(source, dynamics());
        automaton.new_mode(target, dynamics());
        automaton.new_urgent_transition(source, event, target, {next(x) = x}, x >= 0);

        // The evaluator contract must reject an invalid enum rather than
        // treating it as uncertain evidence.
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      HybridDiscreteTimeStepPlanner::plan(
                                          automaton, source, Integer(0),
                                          constant(static_cast<Evidence>(255))));
    }

    Void test_invalid_permissive_evidence() {
        DiscreteTimeHybridAutomaton automaton;
        const DiscreteLocation source(1), target(2);
        const DiscreteEvent event("optional");
        const RealVariable x("x");
        automaton.new_mode(source, dynamics());
        automaton.new_mode(target, dynamics());
        automaton.new_permissive_transition(source, event, target,
                                            {next(x) = x}, x >= 0, x <= 2);

        const auto invalid_guard = [](DiscreteLocation, DiscreteEvent, Role role,
                                      const ContinuousPredicate&) {
            return role == Role::GUARD ? static_cast<Evidence>(255) : Evidence::TRUE;
        };
        const auto invalid_invariant = [](DiscreteLocation, DiscreteEvent, Role role,
                                          const ContinuousPredicate&) {
            return role == Role::PROGRESSION_INVARIANT
                     ? static_cast<Evidence>(255)
                     : Evidence::TRUE;
        };
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      HybridDiscreteTimeStepPlanner::plan(automaton, source, Integer(0), invalid_guard));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      HybridDiscreteTimeStepPlanner::plan(automaton, source, Integer(0), invalid_invariant));
    }

    Void test_permissive_overlap_and_simultaneous_events() {
        DiscreteTimeHybridAutomaton automaton;
        const DiscreteLocation source(1), a(2), b(3);
        const DiscreteEvent first("a"), second("b");
        const RealVariable x("x");
        automaton.new_mode(source, dynamics());
        automaton.new_mode(a, dynamics());
        automaton.new_mode(b, dynamics());
        automaton.new_permissive_transition(source, first, a, {next(x) = x}, x >= 0, x <= 2);
        automaton.new_permissive_transition(source, second, b, {next(x) = x}, x >= 0, x <= 2);

        SizeType guard_calls = 0u, invariant_calls = 0u;
        const auto evaluated = [&](DiscreteLocation, DiscreteEvent event, Role role, const ContinuousPredicate&) {
            if (role == Role::GUARD) {
                ++guard_calls;
                return event == first ? Evidence::TRUE : Evidence::UNKNOWN;
            }
            ++invariant_calls;
            return Evidence::TRUE;
        };
        const auto plan = HybridDiscreteTimeStepPlanner::plan(automaton, source, Integer(7), evaluated);
        ARIADNE_TESTING_EXPECT_EQUAL(guard_calls, 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(invariant_calls, 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(plan.branches.size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(plan.successor_index, Integer(8));
        ARIADNE_TESTING_PROVE(plan.branches[0u].kind == Kind::TRANSITION);
        ARIADNE_TESTING_PROVE(plan.branches[1u].kind == Kind::TRANSITION);
        ARIADNE_TESTING_PROVE(plan.branches[2u].kind == Kind::PROGRESSION);
        ARIADNE_TESTING_PROVE(plan.branches[0u].certain);
        ARIADNE_TESTING_PROVE(!plan.branches[1u].certain);
        ARIADNE_TESTING_PROVE(plan.branches[2u].certain);
        ARIADNE_TESTING_PROVE(plan.branches[0u].target == a);
        ARIADNE_TESTING_PROVE(plan.branches[1u].target == b);
        ARIADNE_TESTING_PROVE(plan.branches[2u].target == source);
    }

    Void test_unknown_invariant_and_blocked_progress() {
        DiscreteTimeHybridAutomaton automaton;
        const DiscreteLocation source(1), target(2);
        const DiscreteEvent event("optional");
        const RealVariable x("x");
        automaton.new_mode(source, dynamics());
        automaton.new_mode(target, dynamics());
        automaton.new_permissive_transition(source, event, target, {next(x) = x}, x >= 0, x <= 1);
        const auto possible = [&](DiscreteLocation, DiscreteEvent, Role role, const ContinuousPredicate&) {
            return role == Role::GUARD ? Evidence::FALSE : Evidence::UNKNOWN;
        };
        const auto first_plan = HybridDiscreteTimeStepPlanner::plan(automaton, source, Integer(0), possible);
        ARIADNE_TESTING_EXPECT_EQUAL(first_plan.branches.size(), 1u);
        ARIADNE_TESTING_PROVE(first_plan.branches[0u].kind == Kind::PROGRESSION);
        ARIADNE_TESTING_PROVE(!first_plan.branches[0u].certain);

        const auto blocked = [&](DiscreteLocation, DiscreteEvent, Role role, const ContinuousPredicate&) {
            return role == Role::GUARD ? Evidence::FALSE : Evidence::FALSE;
        };
        const auto second_plan = HybridDiscreteTimeStepPlanner::plan(automaton, source, Integer(0), blocked);
        ARIADNE_TESTING_PROVE(second_plan.branches.empty());
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_discrete_time_step_planner");
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepPlanner, test_no_events);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepPlanner, test_urgent_complement);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepPlanner, test_invalid_predicate_evidence);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepPlanner, test_invalid_permissive_evidence);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepPlanner, test_permissive_overlap_and_simultaneous_events);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepPlanner, test_unknown_invariant_and_blocked_progress);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
