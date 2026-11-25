/***************************************************************************
 *            test_hybrid_discrete_time_step_image.cpp
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

/** @file tests/test_hybrid_discrete_time_step_image.cpp
 *  @brief Contract tests for validated one-tick deterministic hybrid images.
 */

#include <stdexcept>

#include "evolution/hybrid_discrete_time_step_image.hpp"
#include "function/taylor_function.hpp"
#include "hybrid/hybrid_set.hpp"
#include "numeric/rational.hpp"
#include "logging/test_runner.hpp"
#include "system/discrete_time_hybrid_automaton.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;

class TestHybridDiscreteTimeStepImage : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    static EnclosureConfiguration configuration() {
        return EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8)));
    }

    Void test_progression_and_reset_tick() {
        const RealVariable x("x");
        const DiscreteLocation source(1), destination(2);
        const DiscreteEvent event("jump");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(source, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        automaton.new_mode(destination, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x - 1}));
        automaton.new_permissive_transition(source, event, destination,
                                            {next(x) = x + 3}, x >= 0, x <= 3);

        HybridEnclosure initial(HybridExactBox(source, RealSpace({x}),
                                               ExactBoxType({{0.0_x, 1.0_x}})), configuration());
        const HybridDiscreteTimePlannedBranch progression{
            HybridDiscreteTimeBranchKind::PROGRESSION, source, std::nullopt, true};
        const HybridDiscreteTimePlannedBranch jump{
            HybridDiscreteTimeBranchKind::TRANSITION, destination, event, true};

        const auto evolved = HybridDiscreteTimeStepImage::apply(automaton, initial, progression);
        const auto reset = HybridDiscreteTimeStepImage::apply(automaton, initial, jump);
        ARIADNE_TESTING_PROVE(evolved.location() == source);
        ARIADNE_TESTING_PROVE(reset.location() == destination);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.time_range().lower_bound().raw().get_d(), 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(reset.time_range().lower_bound().raw().get_d(), 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.previous_events().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(reset.previous_events().size(), 1u);
        ARIADNE_TESTING_PROVE(reset.previous_events()[0u] == event);
        ARIADNE_TESTING_PROVE(definitely(evolved.state_bounding_box().inside(
            ExactBoxType({{0.0_x, 3.0_x}}))));
        ARIADNE_TESTING_PROVE(definitely(reset.state_bounding_box().inside(
            ExactBoxType({{2.0_x, 5.0_x}}))));
        ARIADNE_TESTING_EXPECT_EQUAL(initial.time_range().lower_bound().raw().get_d(), 0.0);
    }

    Void test_reset_target_variable_order() {
        const RealVariable x("x"), y("y");
        const DiscreteLocation source(1), destination(2);
        const DiscreteEvent event("swap");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(source, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x, next(y) = y}));
        automaton.new_mode(destination, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x, next(y) = y}));
        automaton.new_urgent_transition(source, event, destination,
                                        {next(y) = x + 10, next(x) = y + 20}, x >= 0);

        HybridEnclosure initial(HybridExactBox(source, RealSpace({x, y}),
                                               ExactBoxType({{0.0_x, 1.0_x}, {2.0_x, 3.0_x}})), configuration());
        const HybridDiscreteTimePlannedBranch branch{
            HybridDiscreteTimeBranchKind::TRANSITION, destination, event, true};
        const auto result = HybridDiscreteTimeStepImage::apply(automaton, initial, branch);
        const auto bounds = result.state_bounding_box();
        ARIADNE_TESTING_PROVE(definitely(bounds.inside(ExactBoxType({
            {21.0_x, 24.0_x}, {9.0_x, 12.0_x}}))));
        ARIADNE_TESTING_EXPECT_EQUAL(result.time_range().lower_bound().raw().get_d(), 1.0);
    }

    Void test_no_silent_inclusion_sampling() {
        const RealVariable x("x"), u("u");
        const DiscreteLocation source(1), uncertain(2);
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(source, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        automaton.new_mode(uncertain, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + u}, RealVariablesBox({-1 <= u <= 1})));
        HybridEnclosure initial(HybridExactBox(source, RealSpace({x}),
                                               ExactBoxType({{0.0_x, 1.0_x}})), configuration());

        const HybridDiscreteTimePlannedBranch progression{
            HybridDiscreteTimeBranchKind::PROGRESSION, source, std::nullopt, true};
        const HybridDiscreteTimePlannedBranch wrong{
            HybridDiscreteTimeBranchKind::PROGRESSION, uncertain, std::nullopt, true};
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            HybridDiscreteTimeStepImage::apply(automaton, initial, wrong));

        const DiscreteEvent event("to_inclusion");
        automaton.new_urgent_transition(source, event, uncertain, {next(x) = x}, x >= 0);
        const HybridDiscreteTimePlannedBranch jump{
            HybridDiscreteTimeBranchKind::TRANSITION, uncertain, event, true};
        // Reset into an inclusion is deterministic at this tick:
        // the destination input box is sampled only on progression.
        const auto arrived = HybridDiscreteTimeStepImage::apply(automaton, initial, jump);
        ARIADNE_TESTING_PROVE(arrived.location() == uncertain);
        ARIADNE_TESTING_EXPECT_EQUAL(arrived.continuous_set().number_of_parameters(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(arrived.time_range().lower_bound().raw().get_d(), 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(arrived.previous_events().size(), 1u);
        ARIADNE_TESTING_PROVE(arrived.previous_events()[0u] == event);
        const HybridDiscreteTimePlannedBranch target_progression{
            HybridDiscreteTimeBranchKind::PROGRESSION, uncertain, std::nullopt, true};
        const auto continued = HybridDiscreteTimeStepImage::apply(
            automaton, arrived, target_progression);
        ARIADNE_TESTING_EXPECT_EQUAL(continued.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(continued.continuous_set().variable_kinds().back(),
                                    EnclosureVariableKind::INPUT);
        ARIADNE_TESTING_EXPECT_EQUAL(continued.time_range().lower_bound().raw().get_d(), 2.0);
        ARIADNE_TESTING_EXPECT_NO_THROW(
            HybridDiscreteTimeStepImage::apply(automaton, initial, progression));

        HybridEnclosure from_inclusion(HybridExactBox(uncertain, RealSpace({x}),
                                               ExactBoxType({{0.0_x, 1.0_x}})), configuration());
        const HybridDiscreteTimePlannedBranch uncertain_progression{
            HybridDiscreteTimeBranchKind::PROGRESSION, uncertain, std::nullopt, true};
        const auto image = HybridDiscreteTimeStepImage::apply(
            automaton, from_inclusion, uncertain_progression);
        ARIADNE_TESTING_EXPECT_EQUAL(image.time_range().lower_bound().raw().get_d(), 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(image.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(image.continuous_set().variable_kinds().back(),
                                    EnclosureVariableKind::INPUT);
        ARIADNE_TESTING_PROVE(definitely(image.state_bounding_box().inside(
            ExactBoxType({{-2.0_x, 3.0_x}}))));
    }

    Void test_reset_into_auxiliary_inclusion_with_fresh_arrival_input() {
        const RealVariable x("x"), u("u"), a("a");
        const DiscreteLocation source(1), target(2);
        const DiscreteEvent event("to_auxiliary_inclusion");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(source, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        automaton.new_mode(target, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + a},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        automaton.new_urgent_transition(source, event, target,
                                         {next(x) = x + 1}, x >= 0);
        const HybridEnclosure initial(HybridExactBox(source, RealSpace({x}),
            ExactBoxType({{0.0_x, 0.0_x}})), configuration());
        const HybridDiscreteTimePlannedBranch jump{
            HybridDiscreteTimeBranchKind::TRANSITION, target, event, true};
        const HybridDiscreteTimePlannedBranch progress{
            HybridDiscreteTimeBranchKind::PROGRESSION, target, std::nullopt, true};
        const auto arrived = HybridDiscreteTimeStepImage::apply(automaton, initial, jump);
        ARIADNE_TESTING_PROVE(arrived.location() == target);
        ARIADNE_TESTING_EXPECT_EQUAL(arrived.auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(arrived.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(arrived.continuous_set().variable_kinds()[1u],
                                     EnclosureVariableKind::INPUT);
        ARIADNE_TESTING_EXPECT_EQUAL(arrived.time_range().lower_bound().raw().get_d(), 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(arrived.previous_events().size(), 1u);
        ARIADNE_TESTING_PROVE(arrived.previous_events()[0u] == event);
        const auto reset_bounds = arrived.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(reset_bounds.lower_bound().raw().get_d() <= 1.0);
        ARIADNE_TESTING_PROVE(reset_bounds.upper_bound().raw().get_d() >= 1.0);
        const auto& arrived_set = arrived.continuous_set();
        ARIADNE_TESTING_EXPECT_EQUAL(arrived_set.auxiliary_arguments(),
                                     EnclosureFunctionArguments::STATE_INPUTS);
        const auto arrival_inputs = arrived_set.input_parameter_function(1u);
        const auto arrival_auxiliaries = arrived_set.auxiliary_function();
        const auto arrival_states = arrived_set.state_function();
        const auto arrival_relation =
            arrival_auxiliaries[0u] - arrival_states[0u] - arrival_inputs[0u];
        const auto arrival_error = arrival_relation.range();
        ARIADNE_TESTING_PROVE(arrival_error.lower_bound().raw().get_d() >= -1e-6);
        ARIADNE_TESTING_PROVE(arrival_error.upper_bound().raw().get_d() <= 1e-6);

        const auto continued = HybridDiscreteTimeStepImage::apply(automaton, arrived, progress);
        ARIADNE_TESTING_EXPECT_EQUAL(continued.continuous_set().number_of_parameters(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(continued.auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(continued.time_range().lower_bound().raw().get_d(), 2.0);
        // x[0]=0, reset x[1]=1, a[1]=1+u[1] with u[1] in [-1,1].
        // The progression yields x[2]=2+u[1], hence the full [1,3].
        const auto bounds = continued.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= 1.0);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 3.0);
        ARIADNE_TESTING_EXPECT_EQUAL(initial.continuous_set().number_of_parameters(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(initial.auxiliary_space().dimension(), 0u);
    }

    Void test_reset_from_auxiliary_inclusion_into_auxiliary_inclusion() {
        const RealVariable x("x"), u("u"), v("v"), a("a"), b("b");
        const DiscreteLocation source(1), target(2);
        const DiscreteEvent event("reset_to_new_input");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(source, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + a},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        automaton.new_mode(target, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + b},
            List<RealAssignment>{let(b) = x + v},
            RealVariablesBox({-2 <= v <= 2})));
        automaton.new_permissive_transition(source, event, target,
                                            {next(x) = x + 3}, x >= 0, x <= 1);
        const HybridEnclosure original(HybridExactBox(source, RealSpace({x}),
            ExactBoxType({{0.0_x, 0.0_x}})), configuration());
        const HybridDiscreteTimePlannedBranch jump{
            HybridDiscreteTimeBranchKind::TRANSITION, target, event, true};
        const auto arrived = HybridDiscreteTimeStepImage::apply(automaton, original, jump);
        ARIADNE_TESTING_PROVE(arrived.location() == target);
        ARIADNE_TESTING_EXPECT_EQUAL(arrived.auxiliary_space(), RealSpace({b}));
        ARIADNE_TESTING_EXPECT_EQUAL(arrived.continuous_set().number_of_parameters(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(arrived.continuous_set().variable_kinds()[1u],
                                     EnclosureVariableKind::INPUT);
        ARIADNE_TESTING_EXPECT_EQUAL(arrived.continuous_set().variable_kinds()[2u],
                                     EnclosureVariableKind::INPUT);
        const auto reset_box = arrived.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(reset_box.lower_bound().raw().get_d() <= 3.0);
        ARIADNE_TESTING_PROVE(reset_box.upper_bound().raw().get_d() >= 3.0);
        const auto arrival_inputs = arrived.continuous_set().input_parameter_function(1u);
        const auto arrival_auxiliaries = arrived.continuous_set().auxiliary_function();
        const auto arrival_states = arrived.continuous_set().state_function();
        const auto relation =
            arrival_auxiliaries[0u] - arrival_states[0u] - arrival_inputs[0u];
        const auto error = relation.range();
        ARIADNE_TESTING_PROVE(error.lower_bound().raw().get_d() >= -1e-6);
        ARIADNE_TESTING_PROVE(error.upper_bound().raw().get_d() <= 1e-6);
        ARIADNE_TESTING_EXPECT_EQUAL(original.continuous_set().number_of_parameters(), 1u);
    }

    Void test_fresh_input_at_every_inclusion_tick() {
        const RealVariable x("x"), u("u");
        const DiscreteLocation mode(1);
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        HybridEnclosure initial(HybridExactBox(mode, RealSpace({x}),
            ExactBoxType({{0.0_x, 0.0_x}})), configuration());
        const HybridDiscreteTimePlannedBranch progress{
            HybridDiscreteTimeBranchKind::PROGRESSION, mode, std::nullopt, true};
        const auto first = HybridDiscreteTimeStepImage::apply(automaton, initial, progress);
        const auto second = HybridDiscreteTimeStepImage::apply(automaton, first, progress);
        ARIADNE_TESTING_EXPECT_EQUAL(first.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(second.continuous_set().number_of_parameters(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(second.continuous_set().variable_kinds()[1u],
                                    EnclosureVariableKind::INPUT);
        ARIADNE_TESTING_EXPECT_EQUAL(second.continuous_set().variable_kinds()[2u],
                                    EnclosureVariableKind::INPUT);
        ARIADNE_TESTING_EXPECT_EQUAL(second.time_range().lower_bound().raw().get_d(), 2.0);
        const auto extent = second.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(extent.lower_bound().raw().get_d() <= -2.0);
        ARIADNE_TESTING_PROVE(extent.upper_bound().raw().get_d() >= 2.0);
        ARIADNE_TESTING_EXPECT_EQUAL(initial.continuous_set().number_of_parameters(), 1u);
    }

    Void test_inclusion_reset_does_not_choose_new_input() {
        const RealVariable x("x"), u("u");
        const DiscreteLocation mode(1), target(2);
        const DiscreteEvent event("jump");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        automaton.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        automaton.new_permissive_transition(mode, event, target,
                                            {next(x) = x + 3}, x >= 1/2_q, x <= 1);
        HybridEnclosure source(HybridExactBox(mode, RealSpace({x}),
            ExactBoxType({{0.0_x, 1.0_x}})), configuration());
        const HybridDiscreteTimePlannedBranch jump{
            HybridDiscreteTimeBranchKind::TRANSITION, target, event, false};
        const HybridDiscreteTimePlannedBranch progress{
            HybridDiscreteTimeBranchKind::PROGRESSION, mode, std::nullopt, false};
        const auto reset = HybridDiscreteTimeStepImage::apply(automaton, source, jump);
        const auto evolved = HybridDiscreteTimeStepImage::apply(automaton, source, progress);
        ARIADNE_TESTING_EXPECT_EQUAL(reset.continuous_set().number_of_parameters(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_PROVE(reset.location() == target);
        ARIADNE_TESTING_PROVE(evolved.location() == mode);
        ARIADNE_TESTING_EXPECT_EQUAL(reset.previous_events().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.previous_events().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(reset.time_range().lower_bound().raw().get_d(), 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.time_range().lower_bound().raw().get_d(), 1.0);
    }

    Void test_event_free_input_dependent_auxiliary() {
        const RealVariable x("x"), u("u"), a("a");
        const DiscreteLocation mode(1);
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + a},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        const HybridDiscreteTimePlannedBranch progress{
            HybridDiscreteTimeBranchKind::PROGRESSION, mode, std::nullopt, true};
        const HybridEnclosure source(HybridExactBox(mode, RealSpace({x}),
            ExactBoxType({{0.0_x, 0.0_x}})), configuration());
        const auto first = HybridDiscreteTimeStepImage::apply(automaton, source, progress);
        const auto second = HybridDiscreteTimeStepImage::apply(automaton, first, progress);
        ARIADNE_TESTING_EXPECT_EQUAL(first.auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(first.continuous_set().number_of_parameters(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(second.continuous_set().number_of_parameters(), 4u);
        for (SizeType i = 1u; i != 4u; ++i)
            ARIADNE_TESTING_EXPECT_EQUAL(second.continuous_set().variable_kinds()[i],
                                         EnclosureVariableKind::INPUT);
        const auto extent = second.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(extent.lower_bound().raw().get_d() <= -3.0);
        ARIADNE_TESTING_PROVE(extent.upper_bound().raw().get_d() >= 3.0);
        const auto& set = second.continuous_set();
        ARIADNE_TESTING_EXPECT_EQUAL(set.auxiliary_arguments(),
                                     EnclosureFunctionArguments::STATE_INPUTS);
        // Keep the vector patch owners alive while forming the scalar
        // residual. Indexing a temporary function patch can leave a
        // borrowed component invalid after the full expression.
        const auto input_functions = set.input_parameter_function(1u);
        const auto auxiliary_functions = set.auxiliary_function();
        const auto state_functions = set.state_function();
        const auto current = input_functions[0u];
        const auto residual = auxiliary_functions[0u] - state_functions[0u] - current;
        const auto range = residual.range();
        ARIADNE_TESTING_PROVE(range.lower_bound().raw().get_d() >= -1e-6);
        ARIADNE_TESTING_PROVE(range.upper_bound().raw().get_d() <= 1e-6);
        ARIADNE_TESTING_EXPECT_EQUAL(source.continuous_set().number_of_parameters(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(source.auxiliary_space().dimension(), 0u);
    }

    Void test_state_only_event_from_auxiliary_inclusion() {
        const RealVariable x("x"), u("u"), a("a");
        const DiscreteLocation mode(1), target(2);
        const DiscreteEvent event("state_only");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + a},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        automaton.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        automaton.new_permissive_transition(mode, event, target,
            {next(x) = x + 3}, x >= 0, x <= 1);
        const HybridEnclosure source(HybridExactBox(mode, RealSpace({x}),
            ExactBoxType({{0.0_x, 0.0_x}})), configuration());
        const HybridDiscreteTimePlannedBranch reset{
            HybridDiscreteTimeBranchKind::TRANSITION, target, event, true};
        const HybridDiscreteTimePlannedBranch progress{
            HybridDiscreteTimeBranchKind::PROGRESSION, mode, std::nullopt, true};
        const auto jumped = HybridDiscreteTimeStepImage::apply(automaton, source, reset);
        const auto advanced = HybridDiscreteTimeStepImage::apply(automaton, source, progress);
        ARIADNE_TESTING_PROVE(jumped.location() == target);
        ARIADNE_TESTING_PROVE(advanced.location() == mode);
        ARIADNE_TESTING_EXPECT_EQUAL(jumped.auxiliary_space().dimension(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(advanced.auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(jumped.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(advanced.continuous_set().number_of_parameters(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(jumped.continuous_set().variable_kinds()[1u],
                                     EnclosureVariableKind::INPUT);
        ARIADNE_TESTING_EXPECT_EQUAL(advanced.continuous_set().variable_kinds()[2u],
                                     EnclosureVariableKind::INPUT);
        const auto reset_bounds = jumped.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(reset_bounds.lower_bound().raw().get_d() <= 3.0);
        ARIADNE_TESTING_PROVE(reset_bounds.upper_bound().raw().get_d() >= 3.0);
        const auto progress_bounds = advanced.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(progress_bounds.lower_bound().raw().get_d() <= -1.0);
        ARIADNE_TESTING_PROVE(progress_bounds.upper_bound().raw().get_d() >= 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(jumped.previous_events().size(), 1u);
        ARIADNE_TESTING_PROVE(jumped.previous_events()[0u] == event);
        ARIADNE_TESTING_EXPECT_EQUAL(jumped.time_range().lower_bound().raw().get_d(), 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(advanced.time_range().lower_bound().raw().get_d(), 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(source.continuous_set().number_of_parameters(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(source.auxiliary_space().dimension(), 0u);
    }

    Void test_input_dependent_event_reset_and_invariant_rejected() {
        const RealVariable x("x"), u("u"), a("a");
        const DiscreteLocation mode(1), target(2);
        const DiscreteEvent reset_event("aux_reset"), invariant_event("aux_invariant");
        const DifferenceAlgebraicInclusion source_system(
            List<PrimedRealAssignment>{next(x) = x + a},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1}));
        const DifferenceAlgebraicEquation target_system(
            List<PrimedRealAssignment>{next(x) = x});
        const HybridEnclosure source(HybridExactBox(mode, RealSpace({x}),
            ExactBoxType({{0.0_x, 1.0_x}})), configuration());

        DiscreteTimeHybridAutomaton reset_automaton;
        reset_automaton.new_mode(mode, source_system);
        reset_automaton.new_mode(target, target_system);
        reset_automaton.new_urgent_transition(mode, reset_event, target,
            {next(x) = a + 3}, x >= 0);
        const HybridDiscreteTimePlannedBranch progress{
            HybridDiscreteTimeBranchKind::PROGRESSION, mode, std::nullopt, false};
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            HybridDiscreteTimeStepImage::apply(reset_automaton, source, progress));

        DiscreteTimeHybridAutomaton invariant_automaton;
        invariant_automaton.new_mode(mode, source_system);
        invariant_automaton.new_mode(target, target_system);
        invariant_automaton.new_permissive_transition(mode, invariant_event, target,
            {next(x) = x}, x >= 0, a <= 1);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            HybridDiscreteTimeStepImage::apply(invariant_automaton, source, progress));
        ARIADNE_TESTING_EXPECT_EQUAL(source.auxiliary_space().dimension(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(source.continuous_set().number_of_parameters(), 1u);
    }

    Void test_input_dependent_auxiliary_event_rejected() {
        const RealVariable x("x"), u("u"), a("a");
        const DiscreteLocation mode(1), target(2);
        const DiscreteEvent event("input_dependent_guard");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        automaton.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        automaton.new_urgent_transition(mode, event, target, {next(x) = x}, a >= 0);
        HybridEnclosure initial(HybridExactBox(mode, RealSpace({x}),
            ExactBoxType({{0.0_x, 1.0_x}})), configuration());
        const HybridDiscreteTimePlannedBranch progress{
            HybridDiscreteTimeBranchKind::PROGRESSION, mode, std::nullopt, true};
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            HybridDiscreteTimeStepImage::apply(automaton, initial, progress));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_discrete_time_step_image");
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepImage, test_progression_and_reset_tick);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepImage, test_reset_target_variable_order);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepImage, test_no_silent_inclusion_sampling);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepImage, test_reset_into_auxiliary_inclusion_with_fresh_arrival_input);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepImage, test_reset_from_auxiliary_inclusion_into_auxiliary_inclusion);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepImage, test_fresh_input_at_every_inclusion_tick);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepImage, test_event_free_input_dependent_auxiliary);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepImage, test_state_only_event_from_auxiliary_inclusion);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepImage, test_input_dependent_event_reset_and_invariant_rejected);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepImage, test_input_dependent_auxiliary_event_rejected);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeStepImage, test_inclusion_reset_does_not_choose_new_input);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
