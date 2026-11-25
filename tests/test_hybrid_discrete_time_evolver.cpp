/***************************************************************************
 *            tests/test_hybrid_discrete_time_evolver.cpp
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

/** @file tests/test_hybrid_discrete_time_evolver.cpp
 *  @brief Public validated discrete hybrid evolver contracts.
 */

#include <memory>
#include <stdexcept>
#include <type_traits>

#include "evolution/hybrid_discrete_time_evolver.hpp"
#include "evolution/traits.hpp"
#include "function/taylor_function.hpp"
#include "hybrid/hybrid_set.hpp"
#include "logging/test_runner.hpp"
#include "numeric/rational.hpp"
#include "system/discrete_time_hybrid_automaton.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;
using Evolver = HybridDiscreteTimeEvolver<ValidatedTag>;
using Termination = HybridDiscreteTimeTerminationCriterion;
static_assert(std::is_base_of_v<Evolver::Interface, Evolver>);
static_assert(std::is_same_v<Evolver::ResultType, Orbit<HybridEnclosure>>);
static_assert(std::is_same_v<EvolutionTraits<ValidatedTag, DiscreteTimeHybridAutomatonInterface>::EvolverType, Evolver>);
template<class T> concept CompleteType = requires { sizeof(T); };
static_assert(!CompleteType<HybridDiscreteTimeEvolver<ApproximateTag>>);

class TestHybridDiscreteTimeEvolver : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    static EnclosureConfiguration configuration() {
        return EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8)));
    }
    static HybridEnclosure initial(DiscreteLocation location) {
        const RealVariable x("x");
        return HybridEnclosure(HybridExactBox(location, RealSpace({x}),
            ExactBoxType({{0.0_x, 1.0_x}})), configuration());
    }

    Void test_orbit_and_branches() {
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
        Evolver evolver(system);
        const auto source = initial(mode);
        const auto orbit = evolver.orbit(source, Termination(Integer(0), Integer(1)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.reach().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.intermediate().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 2u);
        const auto tree = evolver.explore(source, Termination(Integer(0), Integer(1)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(tree.nodes.size(), 3u);
        ARIADNE_TESTING_PROVE(tree.nodes[1u].incoming_branch.has_value());
        ARIADNE_TESTING_PROVE(tree.nodes[1u].incoming_branch->event.has_value());
        ARIADNE_TESTING_PROVE(*tree.nodes[1u].incoming_branch->event == event);
        ARIADNE_TESTING_EXPECT_EQUAL(
            evolver.orbit_with_blocked(source, Termination(Integer(0), Integer(1)),
                                      Semantics::UPPER).blocked.size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(source.time_range().lower_bound().raw().get_d(), 0.0);
    }

    Void test_horizon_and_blocked() {
        const RealVariable x("x");
        const DiscreteLocation mode(1), target(2);
        const DiscreteEvent event("blocked");
        DiscreteTimeHybridAutomaton linear;
        linear.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        Evolver evolver(linear);
        const auto zero = evolver.orbit(initial(mode), Termination(Integer(0), Integer(0)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.reach().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.final().size(), 1u);
        const auto orbit = evolver.orbit(initial(mode), Termination(Integer(0), Integer(3)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.reach().size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.intermediate().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().front().time_range().lower_bound().raw().get_d(), 3.0);

        DiscreteTimeHybridAutomaton blocked;
        blocked.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        blocked.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        blocked.new_permissive_transition(mode, event, target, {next(x) = x}, x >= 2, x <= -1);
        Evolver stopped(blocked);
        const auto projected = stopped.orbit_with_blocked(
            initial(mode), Termination(Integer(0), Integer(2)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(projected.orbit.final().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(projected.blocked.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            projected.blocked.front().time_range().lower_bound().raw().get_d(), 0.0);
    }

    Void test_resume_clone_and_budget() {
        const RealVariable x("x");
        const DiscreteLocation mode(1);
        DiscreteTimeHybridAutomaton system;
        system.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        Evolver evolver(system);
        evolver.configuration().set_maximum_nodes(3u);
        std::unique_ptr<Evolver::Interface> copy(evolver.clone());
        HybridEnclosure shifted = initial(mode);
        shifted.set_time_function(shifted.time_function() + StepSizeType(2, 0u));
        const auto resumed = copy->orbit(shifted, Termination(Integer(2), Integer(4)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(resumed.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(resumed.final().front().time_range().lower_bound().raw().get_d(), 4.0);
        evolver.configuration().set_maximum_nodes(2u);
        ARIADNE_TESTING_EXPECT_THROWS(std::length_error,
            evolver.orbit(shifted, Termination(Integer(2), Integer(4)), Semantics::UPPER));
        ARIADNE_TESTING_EXPECT_EQUAL(
            copy->orbit(shifted, Termination(Integer(2), Integer(4)), Semantics::UPPER).final().size(), 1u);
    }

    Void test_public_event_free_inclusion_upper() {
        const RealVariable x("x"), u("u");
        const DiscreteLocation mode(1);
        DiscreteTimeHybridAutomaton system;
        system.new_mode(mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        Evolver evolver(system);
        const auto orbit = evolver.orbit(
            initial(mode), Termination(Integer(0), Integer(2)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.reach().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            orbit.final().front().continuous_set().number_of_parameters(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            orbit.final().front().time_range().lower_bound().raw().get_d(), 2.0);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(initial(mode), Termination(Integer(0), Integer(1)), Semantics::LOWER));
    }

    Void test_public_guarded_inclusion_upper_orbit() {
        const RealVariable x("x"), u("u");
        const DiscreteLocation mode(1), target(2);
        const DiscreteEvent event("choice");
        DiscreteTimeHybridAutomaton system;
        system.new_mode(mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        system.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        system.new_urgent_transition(mode, event, target,
                                     {next(x) = x + 10}, x >= 1/2_q);
        Evolver evolver(system);
        const auto source = initial(mode);
        const auto orbit = evolver.orbit(source, Termination(Integer(0), Integer(1)),
                                        Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.reach().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 2u);
        const auto tree = evolver.explore(source, Termination(Integer(0), Integer(1)),
                                         Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(tree.nodes.size(), 3u);
        ARIADNE_TESTING_PROVE(tree.nodes[1u].incoming_branch.has_value());
        ARIADNE_TESTING_PROVE(tree.nodes[1u].incoming_branch->event.has_value());
        ARIADNE_TESTING_EXPECT_EQUAL(
            tree.nodes[2u].enclosure.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(source, Termination(Integer(0), Integer(1)), Semantics::LOWER));
    }

    Void test_public_transition_into_inclusion_and_continue_upper() {
        const RealVariable x("x"), u("u");
        const DiscreteLocation source(1), destination(2);
        const DiscreteEvent jump("to_inclusion");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(source, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        automaton.new_mode(destination, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        automaton.new_urgent_transition(source, jump, destination,
                                        {next(x) = x + 10}, x >= -1);
        Evolver evolver(automaton);
        const auto history = evolver.explore(
            initial(source), Termination(Integer(0), Integer(2)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(history.nodes.size(), 3u);
        ARIADNE_TESTING_PROVE(history.nodes[1u].enclosure.location() == destination);
        ARIADNE_TESTING_EXPECT_EQUAL(
            history.nodes[1u].enclosure.continuous_set().number_of_parameters(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            history.nodes[2u].enclosure.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            history.nodes[2u].enclosure.continuous_set().variable_kinds().back(),
            EnclosureVariableKind::INPUT);
        ARIADNE_TESTING_EXPECT_EQUAL(
            history.nodes[2u].enclosure.time_range().lower_bound().raw().get_d(), 2.0);
        ARIADNE_TESTING_PROVE(history.nodes[1u].incoming_branch.has_value());
        ARIADNE_TESTING_PROVE(history.nodes[1u].incoming_branch->event.has_value());
        ARIADNE_TESTING_PROVE(*history.nodes[1u].incoming_branch->event == jump);
        const auto orbit = evolver.orbit(
            initial(source), Termination(Integer(0), Integer(2)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 1u);
        ARIADNE_TESTING_PROVE(orbit.final().front().location() == destination);
    }

    Void test_event_free_input_auxiliary_upper_orbit() {
        const RealVariable x("x"), u("u"), a("a");
        const DiscreteLocation mode(1);
        DiscreteTimeHybridAutomaton system;
        system.new_mode(mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + a},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        Evolver evolver(system);
        const auto source = initial(mode);
        const auto zero = evolver.orbit(source,
            Termination(Integer(0), Integer(0)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.initial().auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(zero.initial().continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.reach().size(), 0u);

        const auto history = evolver.explore(source,
            Termination(Integer(0), Integer(2)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(history.nodes.size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(history.nodes[1u].enclosure.auxiliary_space(), RealSpace({a}));
        const auto& final = history.nodes[2u].enclosure;
        ARIADNE_TESTING_EXPECT_EQUAL(final.continuous_set().number_of_parameters(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(final.auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(final.time_range().lower_bound().raw().get_d(), 2.0);
        const auto extent = final.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(extent.lower_bound().raw().get_d() <= -3.0);
        ARIADNE_TESTING_PROVE(extent.upper_bound().raw().get_d() >= 3.0);
        ARIADNE_TESTING_EXPECT_EQUAL(source.continuous_set().number_of_parameters(), 1u);

        const auto continued = evolver.orbit(final,
            Termination(Integer(2), Integer(3)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(continued.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(continued.final().front().continuous_set().number_of_parameters(), 5u);
        ARIADNE_TESTING_EXPECT_EQUAL(continued.final().front().auxiliary_space(), RealSpace({a}));
    }

    Void test_state_only_event_from_auxiliary_inclusion_upper() {
        const RealVariable x("x"), u("u"), a("a");
        const DiscreteLocation source_mode(1), target_mode(2);
        const DiscreteEvent event("state_guard");
        DiscreteTimeHybridAutomaton system;
        system.new_mode(source_mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + a},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        system.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        // The event depends on x alone. Its predicate and reset must not
        // be evaluated using an unconstrained or sampled auxiliary a.
        system.new_permissive_transition(source_mode, event, target_mode,
            {next(x) = x + 3}, x >= 0, x <= 1);
        Evolver evolver(system);
        const HybridEnclosure original = initial(source_mode);
        const auto history = evolver.explore(original,
            Termination(Integer(0), Integer(1)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(history.nodes.size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(history.nodes[0u].enclosure.auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(history.nodes[0u].enclosure.continuous_set().number_of_parameters(), 2u);

        SizeType reset_count = 0u, progress_count = 0u;
        for (SizeType i = 1u; i != history.nodes.size(); ++i) {
            const auto& successor = history.nodes[i].enclosure;
            ARIADNE_TESTING_EXPECT_EQUAL(successor.time_range().lower_bound().raw().get_d(), 1.0);
            if (successor.location() == target_mode) {
                ++reset_count;
                ARIADNE_TESTING_EXPECT_EQUAL(successor.auxiliary_space().dimension(), 0u);
                ARIADNE_TESTING_EXPECT_EQUAL(successor.continuous_set().number_of_parameters(), 2u);
                ARIADNE_TESTING_EXPECT_EQUAL(successor.previous_events().size(), 1u);
                ARIADNE_TESTING_PROVE(successor.previous_events()[0u] == event);
            } else {
                ++progress_count;
                ARIADNE_TESTING_PROVE(successor.location() == source_mode);
                ARIADNE_TESTING_EXPECT_EQUAL(successor.auxiliary_space(), RealSpace({a}));
                ARIADNE_TESTING_EXPECT_EQUAL(successor.continuous_set().number_of_parameters(), 3u);
                const auto box = successor.state_bounding_box()[0u];
                ARIADNE_TESTING_PROVE(box.lower_bound().raw().get_d() <= -1.0);
                ARIADNE_TESTING_PROVE(box.upper_bound().raw().get_d() >= 1.0);
            }
        }
        ARIADNE_TESTING_EXPECT_EQUAL(reset_count, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(progress_count, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(original.continuous_set().number_of_parameters(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(original.auxiliary_space().dimension(), 0u);
        const auto orb = evolver.orbit(original,
            Termination(Integer(0), Integer(1)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(orb.final().size(), 2u);
    }

    Void test_rejections() {
        const RealVariable x("x"), u("u");
        const DiscreteLocation mode(1);
        DiscreteTimeHybridAutomaton system;
        system.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        Evolver evolver(system);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.configuration().set_maximum_nodes(0u));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            Termination(Integer(2), Integer(1)));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(initial(mode), Termination(Integer(0), Integer(0)), Semantics::LOWER));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(initial(mode), Termination(Integer(1), Integer(1)), Semantics::UPPER));
        DiscreteTimeHybridAutomaton inclusion;
        const RealVariable a("a");
        inclusion.new_mode(mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        const DiscreteLocation target(2);
        const DiscreteEvent event("auxiliary_guard");
        inclusion.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        inclusion.new_urgent_transition(mode, event, target, {next(x) = x}, a >= 0);
        Evolver no_sampling(inclusion);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            no_sampling.orbit(initial(mode), Termination(Integer(0), Integer(1)), Semantics::UPPER));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_discrete_time_evolver");
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_orbit_and_branches);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_horizon_and_blocked);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_resume_clone_and_budget);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_public_event_free_inclusion_upper);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_public_guarded_inclusion_upper_orbit);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_public_transition_into_inclusion_and_continue_upper);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_event_free_input_auxiliary_upper_orbit);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_state_only_event_from_auxiliary_inclusion_upper);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_rejections);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
