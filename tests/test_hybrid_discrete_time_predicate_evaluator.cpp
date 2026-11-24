/***************************************************************************
 *            test_hybrid_discrete_time_predicate_evaluator.cpp
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

/** @file tests/test_hybrid_discrete_time_predicate_evaluator.cpp
 *  @brief Sound guard evidence, correlated auxiliaries and certified tick dispatch.
 */

#include <stdexcept>

#include "../src/hybrid_guard_predicate.hpp"
#include "../src/hybrid_discrete_time_predicate_function.hpp"
#include "evolution/hybrid_discrete_time_predicate_evaluator.hpp"
#include "evolution/hybrid_discrete_time_step_executor.hpp"
#include "function/taylor_function.hpp"
#include "hybrid/hybrid_set.hpp"
#include "logging/test_runner.hpp"
#include "numeric/rational.hpp"
#include "system/discrete_time_hybrid_automaton.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;

class TestHybridDiscreteTimePredicateEvaluator : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    using Evidence = DiscreteTimePredicateEvidence;
    using Kind = HybridDiscreteTimeBranchKind;

    static EnclosureConfiguration configuration() {
        return EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8)));
    }

    static HybridEnclosure initial(DiscreteLocation mode) {
        const RealVariable x("x");
        return HybridEnclosure(HybridExactBox(mode, RealSpace({x}),
                               ExactBoxType({{0.0_x, 1.0_x}})), configuration());
    }

    Void test_certified_signs_and_boundary() {
        const RealVariable x("x");
        const DiscreteLocation mode(1);
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        const auto source = initial(mode);

        ARIADNE_TESTING_PROVE(HybridDiscreteTimePredicateEvaluator::classify(
            automaton, source, x >= -1) == Evidence::TRUE);
        ARIADNE_TESTING_PROVE(HybridDiscreteTimePredicateEvaluator::classify(
            automaton, source, x >= 2) == Evidence::FALSE);
        // x=0 meets the boundary; a strict sign cannot be certified over
        // the entire enclosure, even when the non-strict predicate holds.
        ARIADNE_TESTING_PROVE(HybridDiscreteTimePredicateEvaluator::classify(
            automaton, source, x >= 0) == Evidence::UNKNOWN);
        ARIADNE_TESTING_EXPECT_EQUAL(source.time_range().lower_bound().raw().get_d(), 0.0);
    }

    Void test_correlated_auxiliary_is_not_independent() {
        const RealVariable x("x"), a("a");
        const DiscreteLocation mode(1);
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1},
            List<RealAssignment>{let(a) = x}));
        const auto source = initial(mode);
        // Algebraic substitution proves a-x+1 = 1 for every x in the set.
        ARIADNE_TESTING_PROVE(HybridDiscreteTimePredicateEvaluator::classify(
            automaton, source, a - x + 1 >= 0) == Evidence::TRUE);
        ARIADNE_TESTING_PROVE(HybridDiscreteTimePredicateEvaluator::classify(
            automaton, source, a - x - 1 >= 0) == Evidence::FALSE);
    }

    Void test_discrete_source_guard_uses_shared_certificates() {
        const RealVariable x("x"), a("a");
        const DiscreteLocation mode(1);
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1},
            List<RealAssignment>{let(a) = x}));
        const HybridEnclosure source = initial(mode);
        const auto check = [&](const ContinuousPredicate& predicate,
                               EvolutionDetail::HybridSourceGuardStatus expected,
                               Evidence discrete_expected) {
            const auto indicator =
                EvolutionDetail::discrete_source_indicator(automaton, source, predicate);
            ARIADNE_TESTING_EXPECT_EQUAL(
                EvolutionDetail::classify_source_guard(source, indicator), expected);
            ARIADNE_TESTING_PROVE(
                HybridDiscreteTimePredicateEvaluator::classify(automaton, source, predicate)
                == discrete_expected);
        };
        check(a >= 2, EvolutionDetail::HybridSourceGuardStatus::EXCLUDED, Evidence::FALSE);
        check(a >= 1/2_q, EvolutionDetail::HybridSourceGuardStatus::UNRESOLVED, Evidence::UNKNOWN);
        check(a >= -1, EvolutionDetail::HybridSourceGuardStatus::FULL_SOURCE, Evidence::TRUE);

        const auto indicator = EvolutionDetail::discrete_source_indicator(
            automaton, source, a >= 1/2_q);
        HybridEnclosure restricted = source;
        EvolutionDetail::restrict_source_guard(restricted, DiscreteEvent("discrete"), indicator);
        ARIADNE_TESTING_EXPECT_EQUAL(restricted.number_of_constraints(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(source.number_of_constraints(), 0u);
    }

    Void test_inclusion_state_only_guard_and_urgent_complement() {
        const RealVariable x("x"), u("u");
        const DiscreteLocation mode(1), target(2);
        const DiscreteEvent event("urgent");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        automaton.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        automaton.new_urgent_transition(mode, event, target,
                                        {next(x) = x + 10}, x >= 1/2_q);
        const auto source = initial(mode);
        ARIADNE_TESTING_PROVE(HybridDiscreteTimePredicateEvaluator::classify(
            automaton, source, x >= 2) == Evidence::FALSE);
        ARIADNE_TESTING_PROVE(HybridDiscreteTimePredicateEvaluator::classify(
            automaton, source, x >= -1) == Evidence::TRUE);
        ARIADNE_TESTING_PROVE(HybridDiscreteTimePredicateEvaluator::classify(
            automaton, source, x >= 1/2_q) == Evidence::UNKNOWN);
        const auto branches = HybridDiscreteTimeStepExecutor::execute_certified(
            automaton, source, Integer(0));
        ARIADNE_TESTING_EXPECT_EQUAL(branches.successors.size(), 2u);
        const auto& jump = branches.successors[0u];
        const auto& progress = branches.successors[1u];
        ARIADNE_TESTING_PROVE(jump.branch.kind == Kind::TRANSITION);
        ARIADNE_TESTING_PROVE(progress.branch.kind == Kind::PROGRESSION);
        ARIADNE_TESTING_PROVE(jump.enclosure.location() == target);
        ARIADNE_TESTING_PROVE(progress.enclosure.location() == mode);
        ARIADNE_TESTING_EXPECT_EQUAL(jump.enclosure.number_of_constraints(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(progress.enclosure.number_of_constraints(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(jump.enclosure.continuous_set().number_of_parameters(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(progress.enclosure.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(progress.enclosure.continuous_set().variable_kinds().back(),
                                    EnclosureVariableKind::INPUT);
        ARIADNE_TESTING_EXPECT_EQUAL(jump.enclosure.previous_events().size(), 1u);
        ARIADNE_TESTING_PROVE(jump.enclosure.previous_events()[0u] == event);
        ARIADNE_TESTING_EXPECT_EQUAL(jump.enclosure.time_range().lower_bound().raw().get_d(), 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(progress.enclosure.time_range().lower_bound().raw().get_d(), 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(source.number_of_constraints(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(source.continuous_set().number_of_parameters(), 1u);
    }

    Void test_certified_upper_tick() {
        const RealVariable x("x");
        const DiscreteLocation mode(1), target(2);
        const DiscreteEvent urgent("urgent");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        automaton.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        // Guard x >= 2 is certainly false on [0,1]. No reset branch.
        automaton.new_urgent_transition(mode, urgent, target,
                                        {next(x) = x}, x >= 2);
        const auto source = initial(mode);
        const auto result = HybridDiscreteTimeStepExecutor::execute_certified(
            automaton, source, Integer(0));
        ARIADNE_TESTING_EXPECT_EQUAL(result.successor_index, Integer(1));
        ARIADNE_TESTING_EXPECT_EQUAL(result.successors.size(), 1u);
        ARIADNE_TESTING_PROVE(result.successors[0u].branch.kind == Kind::PROGRESSION);
        ARIADNE_TESTING_PROVE(result.successors[0u].branch.certain);
        ARIADNE_TESTING_PROVE(result.successors[0u].enclosure.location() == mode);
        ARIADNE_TESTING_EXPECT_EQUAL(
            result.successors[0u].enclosure.time_range().lower_bound().raw().get_d(), 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(source.time_range().lower_bound().raw().get_d(), 0.0);
    }

    Void test_boundary_keeps_both_upper_alternatives() {
        const RealVariable x("x");
        const DiscreteLocation mode(1), target(2);
        const DiscreteEvent urgent("urgent");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        automaton.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        automaton.new_urgent_transition(mode, urgent, target,
                                        {next(x) = x}, x >= 0);
        const auto result = HybridDiscreteTimeStepExecutor::execute_certified(
            automaton, initial(mode), Integer(0));
        // The interval touches zero, so neither branch can be dropped by
        // a strict-sign certificate. An overapproximation keeps both.
        ARIADNE_TESTING_EXPECT_EQUAL(result.successors.size(), 2u);
        ARIADNE_TESTING_PROVE(result.successors[0u].branch.kind == Kind::TRANSITION);
        ARIADNE_TESTING_PROVE(result.successors[1u].branch.kind == Kind::PROGRESSION);
        ARIADNE_TESTING_PROVE(!result.successors[0u].branch.certain);
        ARIADNE_TESTING_PROVE(!result.successors[1u].branch.certain);
    }

    Void test_certified_prestate_restriction_and_callback_compatibility() {
        const RealVariable x("x"), a("a");
        const DiscreteLocation mode(1), target(2);
        const DiscreteEvent event("split");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1},
            List<RealAssignment>{let(a) = x}));
        automaton.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        // A correlated auxiliary determines the split. The guard and
        // invariant cover different pre-state subsets at a nonzero boundary.
        automaton.new_permissive_transition(mode, event, target,
                                            {next(x) = x + 10},
                                            a >= 1/2_q, x <= 1/4_q);
        const HybridEnclosure source = initial(mode);
        const auto restricted = HybridDiscreteTimeStepExecutor::execute_certified(
            automaton, source, Integer(0));
        ARIADNE_TESTING_EXPECT_EQUAL(restricted.successors.size(), 2u);
        ARIADNE_TESTING_PROVE(restricted.successors[0u].branch.kind == Kind::TRANSITION);
        ARIADNE_TESTING_PROVE(restricted.successors[1u].branch.kind == Kind::PROGRESSION);
        ARIADNE_TESTING_PROVE(!restricted.successors[0u].branch.certain);
        ARIADNE_TESTING_PROVE(!restricted.successors[1u].branch.certain);
        ARIADNE_TESTING_EXPECT_EQUAL(restricted.successors[0u].enclosure.number_of_constraints(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(restricted.successors[1u].enclosure.number_of_constraints(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(source.number_of_constraints(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(source.previous_events().size(), 0u);

        // The legacy callback path still images the full enclosure:
        // evidence supplied by a caller is not a numerical cut certificate.
        const auto unknown = [](DiscreteLocation, DiscreteEvent,
                                DiscreteTimePredicateRole, const ContinuousPredicate&) {
            return Evidence::UNKNOWN;
        };
        const auto full = HybridDiscreteTimeStepExecutor::execute(
            automaton, source, Integer(0), unknown);
        ARIADNE_TESTING_EXPECT_EQUAL(full.successors.size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(full.successors[0u].enclosure.number_of_constraints(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(full.successors[1u].enclosure.number_of_constraints(), 0u);
    }

    Void test_urgent_complement_restriction() {
        const RealVariable x("x");
        const DiscreteLocation mode(1), target(2);
        const DiscreteEvent event("urgent");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        automaton.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        automaton.new_urgent_transition(mode, event, target,
                                        {next(x) = x + 10}, x >= 1/2_q);
        const auto result = HybridDiscreteTimeStepExecutor::execute_certified(
            automaton, initial(mode), Integer(0));
        ARIADNE_TESTING_EXPECT_EQUAL(result.successors.size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.successors[0u].enclosure.number_of_constraints(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.successors[1u].enclosure.number_of_constraints(), 1u);
        ARIADNE_TESTING_PROVE(result.successors[0u].branch.kind == Kind::TRANSITION);
        ARIADNE_TESTING_PROVE(result.successors[1u].branch.kind == Kind::PROGRESSION);
    }

    Void test_supported_inclusion_and_rejected_cases() {
        const RealVariable x("x"), u("u"), y("y"), a("a");
        const DiscreteLocation mode(1);
        DiscreteTimeHybridAutomaton inclusion;
        inclusion.new_mode(mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        const auto source = initial(mode);

        // The state-only predicate is now supported for a bounded-input
        // inclusion without auxiliaries. On x in [0,1], x >= 0 touches
        // the boundary and remains UNKNOWN rather than a false certificate.
        ARIADNE_TESTING_PROVE(HybridDiscreteTimePredicateEvaluator::classify(
            inclusion, source, x >= 0) == Evidence::UNKNOWN);
        ARIADNE_TESTING_PROVE(HybridDiscreteTimePredicateEvaluator::classify(
            inclusion, source, x >= -1) == Evidence::TRUE);
        ARIADNE_TESTING_PROVE(HybridDiscreteTimePredicateEvaluator::classify(
            inclusion, source, x >= 2) == Evidence::FALSE);

        // Keep an actual unsupported inclusion rejection in the same test:
        // input-dependent auxiliary outputs do not yet have a numerical
        // representation for certified hybrid discrete guard evaluation.
        DiscreteTimeHybridAutomaton with_auxiliary;
        with_auxiliary.new_mode(mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            HybridDiscreteTimePredicateEvaluator::classify(
                with_auxiliary, source, x >= 0));

        // A mismatched state space is still a contract violation.
        DiscreteTimeHybridAutomaton equation;
        equation.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(y) = y}));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            HybridDiscreteTimePredicateEvaluator::classify(equation, source, y >= 0));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_discrete_time_predicate_evaluator");
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimePredicateEvaluator, test_certified_signs_and_boundary);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimePredicateEvaluator, test_correlated_auxiliary_is_not_independent);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimePredicateEvaluator, test_discrete_source_guard_uses_shared_certificates);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimePredicateEvaluator, test_inclusion_state_only_guard_and_urgent_complement);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimePredicateEvaluator, test_certified_upper_tick);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimePredicateEvaluator, test_boundary_keeps_both_upper_alternatives);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimePredicateEvaluator, test_certified_prestate_restriction_and_callback_compatibility);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimePredicateEvaluator, test_urgent_complement_restriction);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimePredicateEvaluator, test_supported_inclusion_and_rejected_cases);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
