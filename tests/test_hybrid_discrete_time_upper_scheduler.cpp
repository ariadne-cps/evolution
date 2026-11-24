/***************************************************************************
 *            test_hybrid_discrete_time_upper_scheduler.cpp
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

/** @file tests/test_hybrid_discrete_time_upper_scheduler.cpp
 *  @brief Finite-horizon, branching, blocked and budget UPPER regressions.
 */

#include <stdexcept>

#include "evolution/hybrid_discrete_time_upper_scheduler.hpp"
#include "function/taylor_function.hpp"
#include "hybrid/hybrid_set.hpp"
#include "logging/test_runner.hpp"
#include "numeric/rational.hpp"
#include "system/discrete_time_hybrid_automaton.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;

class TestHybridDiscreteTimeUpperScheduler : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    using Kind = HybridDiscreteTimeBranchKind;
    using Terminal = HybridDiscreteTimeUpperTerminalKind;

    static EnclosureConfiguration configuration() {
        return EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8)));
    }

    static HybridEnclosure initial(DiscreteLocation location) {
        const RealVariable x("x");
        return HybridEnclosure(HybridExactBox(location, RealSpace({x}),
                               ExactBoxType({{0.0_x, 1.0_x}})), configuration());
    }

    Void test_horizon_and_linear_provenance() {
        const RealVariable x("x");
        const DiscreteLocation mode(1);
        DiscreteTimeHybridAutomaton system;
        system.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        const auto source = initial(mode);
        const auto zero = HybridDiscreteTimeUpperScheduler::evolve(
            system, source, Integer(0), Integer(0), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.nodes.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.terminals.size(), 1u);
        ARIADNE_TESTING_PROVE(zero.terminals[0u].kind == Terminal::HORIZON);
        ARIADNE_TESTING_PROVE(!zero.nodes[0u].parent.has_value());

        const auto full = HybridDiscreteTimeUpperScheduler::evolve(
            system, source, Integer(0), Integer(3), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(full.nodes.size(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(full.terminals.size(), 1u);
        ARIADNE_TESTING_PROVE(full.terminals[0u].kind == Terminal::HORIZON);
        ARIADNE_TESTING_EXPECT_EQUAL(full.terminals[0u].node_index, 3u);
        for (SizeType j = 1u; j < 4u; ++j) {
            ARIADNE_TESTING_EXPECT_EQUAL(full.nodes[j].index, Integer(j));
            ARIADNE_TESTING_PROVE(full.nodes[j].parent.has_value());
            ARIADNE_TESTING_EXPECT_EQUAL(*full.nodes[j].parent, j - 1u);
            ARIADNE_TESTING_PROVE(full.nodes[j].incoming_branch.has_value());
            ARIADNE_TESTING_PROVE(full.nodes[j].incoming_branch->kind == Kind::PROGRESSION);
            ARIADNE_TESTING_PROVE(full.nodes[j].incoming_branch->certain);
            ARIADNE_TESTING_EXPECT_EQUAL(
                full.nodes[j].enclosure.time_range().lower_bound().raw().get_d(), double(j));
        }
        ARIADNE_TESTING_EXPECT_EQUAL(source.time_range().lower_bound().raw().get_d(), 0.0);
    }

    Void test_branch_tree_and_horizon() {
        const RealVariable x("x");
        const DiscreteLocation mode(1), target(2);
        const DiscreteEvent event("choice");
        DiscreteTimeHybridAutomaton system;
        system.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        system.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        system.new_permissive_transition(mode, event, target, {next(x) = x + 3},
                                         x >= 0, x <= 2);
        const auto full = HybridDiscreteTimeUpperScheduler::evolve(
            system, initial(mode), Integer(0), Integer(1), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(full.nodes.size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(full.terminals.size(), 2u);
        ARIADNE_TESTING_PROVE(full.terminals[0u].kind == Terminal::HORIZON);
        ARIADNE_TESTING_PROVE(full.terminals[1u].kind == Terminal::HORIZON);
        ARIADNE_TESTING_PROVE(full.nodes[1u].parent.has_value());
        ARIADNE_TESTING_PROVE(full.nodes[2u].parent.has_value());
        ARIADNE_TESTING_EXPECT_EQUAL(*full.nodes[1u].parent, 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(*full.nodes[2u].parent, 0u);
        ARIADNE_TESTING_PROVE(full.nodes[1u].incoming_branch.has_value());
        ARIADNE_TESTING_PROVE(full.nodes[1u].incoming_branch->kind == Kind::TRANSITION);
        ARIADNE_TESTING_PROVE(full.nodes[1u].incoming_branch->event.has_value());
        ARIADNE_TESTING_PROVE(*full.nodes[1u].incoming_branch->event == event);
        ARIADNE_TESTING_PROVE(full.nodes[1u].enclosure.location() == target);
        ARIADNE_TESTING_PROVE(full.nodes[2u].incoming_branch->kind == Kind::PROGRESSION);
        ARIADNE_TESTING_PROVE(full.nodes[2u].enclosure.location() == mode);
        ARIADNE_TESTING_EXPECT_EQUAL(full.nodes[1u].enclosure.previous_events().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(full.nodes[2u].enclosure.previous_events().size(), 0u);
    }

    Void test_blocked_and_budget_contract() {
        const RealVariable x("x");
        const DiscreteLocation mode(1), target(2);
        const DiscreteEvent event("blocked");
        DiscreteTimeHybridAutomaton blocked;
        blocked.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        blocked.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        blocked.new_permissive_transition(mode, event, target, {next(x) = x},
                                          x >= 2, x <= -1);
        const auto result = HybridDiscreteTimeUpperScheduler::evolve(
            blocked, initial(mode), Integer(0), Integer(2), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.nodes.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.terminals.size(), 1u);
        ARIADNE_TESTING_PROVE(result.terminals[0u].kind == Terminal::BLOCKED);

        DiscreteTimeHybridAutomaton unblocked;
        unblocked.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        ARIADNE_TESTING_EXPECT_THROWS(std::length_error,
            HybridDiscreteTimeUpperScheduler::evolve(
                unblocked, initial(mode), Integer(0), Integer(2), 2u));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            HybridDiscreteTimeUpperScheduler::evolve(
                unblocked, initial(mode), Integer(2), Integer(1), 5u));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            HybridDiscreteTimeUpperScheduler::evolve(
                unblocked, initial(mode), Integer(0), Integer(2), 0u));
    }

    Void test_inclusion_rejected_without_sampling() {
        const RealVariable x("x"), u("u");
        const DiscreteLocation mode(1);
        DiscreteTimeHybridAutomaton system;
        system.new_mode(mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            HybridDiscreteTimeUpperScheduler::evolve(
                system, initial(mode), Integer(0), Integer(1), 10u));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_discrete_time_upper_scheduler");
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeUpperScheduler, test_horizon_and_linear_provenance);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeUpperScheduler, test_branch_tree_and_horizon);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeUpperScheduler, test_blocked_and_budget_contract);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeUpperScheduler, test_inclusion_rejected_without_sampling);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
