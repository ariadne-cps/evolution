/***************************************************************************
 *            test_hybrid_discrete_time_step_executor.cpp
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

/** @file tests/test_hybrid_discrete_time_step_executor.cpp
 *  @brief Regression of deterministic UPPER hybrid tick fan-out and provenance.
 */

#include <stdexcept>

#include "evolution/hybrid_discrete_time_step_executor.hpp"
#include "function/taylor_function.hpp"
#include "hybrid/hybrid_set.hpp"
#include "logging/test_runner.hpp"
#include "numeric/rational.hpp"
#include "system/discrete_time_hybrid_automaton.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;

class TestHybridDiscreteTimeStepExecutor : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    using Evidence = DiscreteTimePredicateEvidence;
    using Role = DiscreteTimePredicateRole;
    using Kind = HybridDiscreteTimeBranchKind;

    static EnclosureConfiguration configuration() {
        return EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8)));
    }

    static HybridEnclosure initial(DiscreteLocation source) {
        const RealVariable x("x");
        return HybridEnclosure(HybridExactBox(source, RealSpace({x}),
                               ExactBoxType({{0.0_x, 1.0_x}})), configuration());
    }

    Void test_independent_permissive_branches() {
        const RealVariable x("x");
        const DiscreteLocation source(1), destination(2);
        const DiscreteEvent event("jump");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(source, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        automaton.new_mode(destination, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x - 1}));
        automaton.new_permissive_transition(source, event, destination,
                                            {next(x) = x + 3}, x >= 0, x <= 2);

        const auto input = initial(source);
        const auto evidence = [](DiscreteLocation, DiscreteEvent, Role,
                                 const ContinuousPredicate&) { return Evidence::UNKNOWN; };
        const auto result = HybridDiscreteTimeStepExecutor::execute(
            automaton, input, Integer(3), evidence);
        ARIADNE_TESTING_EXPECT_EQUAL(result.source_index, Integer(3));
        ARIADNE_TESTING_EXPECT_EQUAL(result.successor_index, Integer(4));
        ARIADNE_TESTING_EXPECT_EQUAL(result.successors.size(), 2u);

        const auto& jump = result.successors[0u];
        const auto& progress = result.successors[1u];
        ARIADNE_TESTING_PROVE(jump.branch.kind == Kind::TRANSITION);
        ARIADNE_TESTING_PROVE(progress.branch.kind == Kind::PROGRESSION);
        ARIADNE_TESTING_PROVE(!jump.branch.certain && !progress.branch.certain);
        ARIADNE_TESTING_PROVE(jump.enclosure.location() == destination);
        ARIADNE_TESTING_PROVE(progress.enclosure.location() == source);
        ARIADNE_TESTING_PROVE(jump.branch.event.has_value());
        ARIADNE_TESTING_PROVE(*jump.branch.event == event);
        ARIADNE_TESTING_PROVE(jump.enclosure.previous_events().size() == 1u);
        ARIADNE_TESTING_PROVE(progress.enclosure.previous_events().size() == 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(jump.enclosure.time_range().lower_bound().raw().get_d(), 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(progress.enclosure.time_range().lower_bound().raw().get_d(), 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(input.time_range().lower_bound().raw().get_d(), 0.0);
        ARIADNE_TESTING_PROVE(definitely(jump.enclosure.state_bounding_box().inside(
            ExactBoxType({{2.0_x, 5.0_x}}))));
        ARIADNE_TESTING_PROVE(definitely(progress.enclosure.state_bounding_box().inside(
            ExactBoxType({{0.0_x, 3.0_x}}))));
    }

    Void test_urgent_branch_alternatives() {
        const RealVariable x("x");
        const DiscreteLocation source(1), destination(2);
        const DiscreteEvent event("urgent");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(source, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        automaton.new_mode(destination, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        automaton.new_urgent_transition(source, event, destination,
                                        {next(x) = x + 2}, x >= 0);

        const auto source_set = initial(source);
        const auto enabled = [](DiscreteLocation, DiscreteEvent, Role,
                                const ContinuousPredicate&) { return Evidence::TRUE; };
        const auto disabled = [](DiscreteLocation, DiscreteEvent, Role,
                                 const ContinuousPredicate&) { return Evidence::FALSE; };
        const auto true_result = HybridDiscreteTimeStepExecutor::execute(
            automaton, source_set, Integer(0), enabled);
        ARIADNE_TESTING_EXPECT_EQUAL(true_result.successors.size(), 1u);
        ARIADNE_TESTING_PROVE(true_result.successors[0u].branch.kind == Kind::TRANSITION);
        ARIADNE_TESTING_PROVE(true_result.successors[0u].branch.certain);
        const auto false_result = HybridDiscreteTimeStepExecutor::execute(
            automaton, source_set, Integer(0), disabled);
        ARIADNE_TESTING_EXPECT_EQUAL(false_result.successors.size(), 1u);
        ARIADNE_TESTING_PROVE(false_result.successors[0u].branch.kind == Kind::PROGRESSION);
        ARIADNE_TESTING_PROVE(false_result.successors[0u].branch.certain);
    }

    Void test_no_progress_and_no_jump() {
        const RealVariable x("x");
        const DiscreteLocation source(1), destination(2);
        const DiscreteEvent event("blocked");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(source, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        automaton.new_mode(destination, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        automaton.new_permissive_transition(source, event, destination,
                                            {next(x) = x}, x >= 0, x <= 0);
        const auto evidence = [](DiscreteLocation, DiscreteEvent, Role role,
                                 const ContinuousPredicate&) {
            return role == Role::GUARD ? Evidence::FALSE : Evidence::FALSE;
        };
        const auto result = HybridDiscreteTimeStepExecutor::execute(
            automaton, initial(source), Integer(8), evidence);
        ARIADNE_TESTING_EXPECT_EQUAL(result.source_index, Integer(8));
        ARIADNE_TESTING_EXPECT_EQUAL(result.successor_index, Integer(9));
        ARIADNE_TESTING_PROVE(result.successors.empty());
    }

    Void test_no_silent_inclusion_images() {
        const RealVariable x("x"), u("u");
        const DiscreteLocation source(1);
        DiscreteTimeHybridAutomaton automaton;
        const RealVariable a("a");
        automaton.new_mode(source, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        const DiscreteLocation target(2);
        const DiscreteEvent event("auxiliary_guard");
        automaton.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        automaton.new_urgent_transition(source, event, target, {next(x) = x}, a >= 0);
        const auto evidence = [](DiscreteLocation, DiscreteEvent, Role,
                                 const ContinuousPredicate&) { return Evidence::TRUE; };
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            HybridDiscreteTimeStepExecutor::execute(
                automaton, initial(source), Integer(0), evidence));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_discrete_time_step_executor");
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepExecutor, test_independent_permissive_branches);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepExecutor, test_urgent_branch_alternatives);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepExecutor, test_no_progress_and_no_jump);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepExecutor, test_no_silent_inclusion_images);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
