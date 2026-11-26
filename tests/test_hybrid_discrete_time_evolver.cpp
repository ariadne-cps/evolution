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
        const auto lower = evolver.orbit(
            initial(mode), Termination(Integer(0), Integer(2)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(lower.final().size() == 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.reach().size(), 2u);
        const auto& lf = lower.final().front();
        ARIADNE_TESTING_EXPECT_EQUAL(lf.continuous_set().number_of_parameters(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(lf.number_of_constraints(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(lf.previous_events().size(), 0u);
        const auto full_domain = lf.continuous_set().parameter_domain();
        for (SizeType i = 1u; i != 3u; ++i) {
            ARIADNE_TESTING_EXPECT_EQUAL(lf.continuous_set().variable_kinds()[i],
                                         EnclosureVariableKind::INPUT);
            ARIADNE_TESTING_PROVE(full_domain[i].lower_bound().raw().get_d() <= -1.0);
            ARIADNE_TESTING_PROVE(full_domain[i].upper_bound().raw().get_d() >= 1.0);
        }
        const auto bounds = lf.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= -2.0);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 3.0);
        const auto zero = evolver.orbit(
            initial(mode), Termination(Integer(0), Integer(0)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            zero.initial().continuous_set().number_of_parameters(), 1u);
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

        const auto lower = evolver.orbit(
            initial(source), Termination(Integer(0), Integer(2)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(lower.final().size() == 1u);
        const auto& lf = lower.final().front();
        ARIADNE_TESTING_PROVE(lf.location() == destination);
        ARIADNE_TESTING_EXPECT_EQUAL(lf.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(lf.number_of_constraints(), 0u);
        ARIADNE_TESTING_REQUIRE(lf.previous_events().size() == 1u);
        ARIADNE_TESTING_PROVE(lf.previous_events()[0u] == jump);
        const auto bounds = lf.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= 9.0);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 12.0);
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

        const auto lower_zero = evolver.orbit(source,
            Termination(Integer(0), Integer(0)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(lower_zero.initial().auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(
            lower_zero.initial().continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(lower_zero.initial().number_of_constraints(), 0u);
        const auto lower = evolver.orbit(source,
            Termination(Integer(0), Integer(2)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(lower.final().size() == 1u);
        const auto& lf = lower.final().front();
        ARIADNE_TESTING_EXPECT_EQUAL(lf.continuous_set().number_of_parameters(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(lf.auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(lf.number_of_constraints(), 0u);
        const auto bounds = lf.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= -3.0);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 7.0);
        const auto resumed = evolver.orbit(lf,
            Termination(Integer(2), Integer(3)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(resumed.final().size() == 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            resumed.final().front().continuous_set().number_of_parameters(), 5u);
        ARIADNE_TESTING_EXPECT_EQUAL(resumed.final().front().number_of_constraints(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(source.auxiliary_space().dimension(), 0u);
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

    Void test_public_reset_into_inclusion_with_auxiliary_and_continue() {
        const RealVariable x("x"), u("u"), a("a");
        const DiscreteLocation source_mode(1), target_mode(2);
        const DiscreteEvent jump("arrival_input");
        DiscreteTimeHybridAutomaton automaton;
        automaton.new_mode(source_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        automaton.new_mode(target_mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + a},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        // Strictly enabled on the entire initial x in [0,1], so
        // this horizon has one transition branch and one continuation.
        automaton.new_urgent_transition(source_mode, jump, target_mode,
                                        {next(x) = x + 3}, x >= -1);
        Evolver evolver(automaton);
        const HybridEnclosure source = initial(source_mode);
        const auto history = evolver.explore(source,
            Termination(Integer(0), Integer(2)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(history.nodes.size(), 3u);
        ARIADNE_TESTING_PROVE(history.nodes[1u].enclosure.location() == target_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(
            history.nodes[1u].enclosure.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            history.nodes[1u].enclosure.auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(
            history.nodes[1u].enclosure.continuous_set().variable_kinds()[1u],
            EnclosureVariableKind::INPUT);
        ARIADNE_TESTING_PROVE(history.nodes[1u].incoming_branch.has_value());
        ARIADNE_TESTING_PROVE(history.nodes[1u].incoming_branch->event.has_value());
        ARIADNE_TESTING_PROVE(*history.nodes[1u].incoming_branch->event == jump);
        const auto& final = history.nodes[2u].enclosure;
        ARIADNE_TESTING_PROVE(final.location() == target_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(final.continuous_set().number_of_parameters(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(final.auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(
            final.continuous_set().variable_kinds()[2u], EnclosureVariableKind::INPUT);
        ARIADNE_TESTING_EXPECT_EQUAL(
            final.time_range().lower_bound().raw().get_d(), 2.0);
        const auto bounds = final.state_bounding_box()[0u];
        // Source x[0] in [0,1], reset x[1]=x[0]+3 in [3,4],
        // then x[2]=2*x[1]+u[1] covers the entire [5,9].
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= 5.0);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 9.0);
        const auto result = evolver.orbit(source,
            Termination(Integer(0), Integer(2)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(result.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(source.continuous_set().number_of_parameters(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(source.auxiliary_space().dimension(), 0u);

        const auto lower = evolver.orbit(source,
            Termination(Integer(0), Integer(2)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(lower.final().size() == 1u);
        const auto& lf = lower.final().front();
        ARIADNE_TESTING_PROVE(lf.location() == target_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(lf.auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(lf.continuous_set().number_of_parameters(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(lf.number_of_constraints(), 0u);
        ARIADNE_TESTING_REQUIRE(lf.previous_events().size() == 1u);
        ARIADNE_TESTING_PROVE(lf.previous_events()[0u] == jump);
        const auto lower_bounds = lf.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(lower_bounds.lower_bound().raw().get_d() <= 5.0);
        ARIADNE_TESTING_PROVE(lower_bounds.upper_bound().raw().get_d() >= 9.0);
    }

    Void test_lower_inclusion_preserves_full_input_without_splitting() {
        const RealVariable x("x"), u("u"), a("a");
        const DiscreteLocation source_mode(1), target_mode(2);
        const DiscreteEvent event("state_only_guard");
        DiscreteTimeHybridAutomaton system;
        system.new_mode(source_mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + a},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        system.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        system.new_urgent_transition(source_mode, event, target_mode,
            {next(x) = x + a + 3}, x >= 1 / 2_q);
        Evolver evolver(system);

        // Guard-boundary sources cannot be split to publish one
        // independently claimed LOWER witness.
        const HybridEnclosure undecided_reset(HybridExactBox(
            source_mode, RealSpace({x}), ExactBoxType({{0.5_x, 1.0_x}})),
            configuration());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(undecided_reset, Termination(Integer(0), Integer(1)),
                          Semantics::LOWER));
        const HybridEnclosure undecided_progression(HybridExactBox(
            source_mode, RealSpace({x}), ExactBoxType({{0.0_x, 0.5_x}})),
            configuration());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(undecided_progression, Termination(Integer(0), Integer(1)),
                          Semantics::LOWER));
        // Caller-supplied uniformly guard-true/false whole sources
        // still preserve every input value and the unchanged map.
        const HybridEnclosure reset_source(HybridExactBox(
            source_mode, RealSpace({x}), ExactBoxType({{0.75_x, 1.0_x}})),
            configuration());
        const auto reset_orbit = evolver.orbit(
            reset_source, Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(reset_orbit.final().size() == 1u);
        const auto& reset_final = reset_orbit.final().front();
        ARIADNE_TESTING_PROVE(reset_final.location() == target_mode);
        ARIADNE_TESTING_REQUIRE(reset_final.previous_events().size() == 1u);
        ARIADNE_TESTING_PROVE(reset_final.previous_events()[0u] == event);
        ARIADNE_TESTING_EXPECT_EQUAL(reset_final.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            reset_final.continuous_set().variable_kinds()[1u], EnclosureVariableKind::INPUT);
        ARIADNE_TESTING_EXPECT_EQUAL(reset_final.number_of_constraints(), 0u);
        const auto reset_input = reset_final.continuous_set().parameter_domain()[1u];
        ARIADNE_TESTING_PROVE(reset_input.lower_bound().raw().get_d() <= -1.0);
        ARIADNE_TESTING_PROVE(reset_input.upper_bound().raw().get_d() >= 1.0);
        // With x=7/8 and u=-1 or +1, the same refined reset contains
        // both 3.75 and 5.75; no input value was selected.
        const auto reset_bounds = reset_final.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(reset_bounds.lower_bound().raw().get_d() <= 3.75);
        ARIADNE_TESTING_PROVE(reset_bounds.upper_bound().raw().get_d() >= 5.75);
        ARIADNE_TESTING_EXPECT_EQUAL(reset_source.number_of_constraints(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            reset_source.continuous_set().number_of_parameters(), 1u);

        // On the other state side, a strictly negative guard certifies
        // progression for the SAME entire input box at the next tick.
        const HybridEnclosure progression_source(HybridExactBox(
            source_mode, RealSpace({x}), ExactBoxType({{0.0_x, 0.25_x}})),
            configuration());
        const auto progression_orbit = evolver.orbit(
            progression_source, Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(progression_orbit.final().size() == 1u);
        const auto& progression_final = progression_orbit.final().front();
        ARIADNE_TESTING_PROVE(progression_final.location() == source_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(progression_final.previous_events().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            progression_final.continuous_set().number_of_parameters(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(progression_final.number_of_constraints(), 0u);
        const auto progression_input =
            progression_final.continuous_set().parameter_domain()[1u];
        ARIADNE_TESTING_PROVE(progression_input.lower_bound().raw().get_d() <= -1.0);
        ARIADNE_TESTING_PROVE(progression_input.upper_bound().raw().get_d() >= 1.0);
        const auto progression_bounds = progression_final.state_bounding_box()[0u];
        // At x=1/8 and both extreme inputs, 2*x+u is -3/4 or 5/4.
        ARIADNE_TESTING_PROVE(progression_bounds.lower_bound().raw().get_d() <= -0.75);
        ARIADNE_TESTING_PROVE(progression_bounds.upper_bound().raw().get_d() >= 1.25);
        ARIADNE_TESTING_EXPECT_EQUAL(progression_source.number_of_constraints(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            progression_source.continuous_set().number_of_parameters(), 1u);
    }

    Void test_lower_uniform_auxiliary_guard_over_full_input_box() {
        const RealVariable x("x"), u("u"), a("a");
        const DiscreteLocation source_mode(1), target_mode(2);
        const DiscreteEvent event("uniform_auxiliary_guard");
        DiscreteTimeHybridAutomaton system;
        system.new_mode(source_mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + a},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        system.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        system.new_urgent_transition(source_mode, event, target_mode,
            {next(x) = x + a + 3}, a >= -2);
        Evolver evolver(system);
        const auto source = initial(source_mode);
        const auto lower = evolver.orbit(
            source, Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(lower.final().size() == 1u);
        const auto& lf = lower.final().front();
        ARIADNE_TESTING_PROVE(lf.location() == target_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(lf.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(lf.number_of_constraints(), 0u);
        ARIADNE_TESTING_REQUIRE(lf.previous_events().size() == 1u);
        ARIADNE_TESTING_PROVE(lf.previous_events()[0u] == event);
        const auto domain = lf.continuous_set().parameter_domain();
        ARIADNE_TESTING_PROVE(domain[1u].lower_bound().raw().get_d() <= -1.0);
        ARIADNE_TESTING_PROVE(domain[1u].upper_bound().raw().get_d() >= 1.0);
        const auto bounds = lf.state_bounding_box()[0u];
        // Correlated reset x'=x+a+3 = 2*x+u+3, over [2,6].
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= 2.0);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 6.0);
        ARIADNE_TESTING_EXPECT_EQUAL(source.auxiliary_space().dimension(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(source.continuous_set().number_of_parameters(), 1u);

        DiscreteTimeHybridAutomaton ambiguous;
        ambiguous.new_mode(source_mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + a},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        ambiguous.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        ambiguous.new_urgent_transition(source_mode, event, target_mode,
            {next(x) = x + a + 3}, a >= 0);
        Evolver undecided(ambiguous);
        // A varying guard is not a licence to sample u or prefer a branch.
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            undecided.orbit(source, Termination(Integer(0), Integer(1)),
                            Semantics::LOWER));

        // The SAME guard may be optional: with a uniformly safe
        // progression invariant, LOWER ignores PERMISSIVE reset and
        // applies x'=x+a to ALL x and ALL u in [-1,1].
        DiscreteTimeHybridAutomaton optional_model;
        optional_model.new_mode(source_mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + a},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        optional_model.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        optional_model.new_permissive_transition(source_mode, event, target_mode,
            {next(x) = x + a + 3}, a >= 0, a <= 3);
        Evolver optional_evolver(optional_model);
        const auto optional_orbit = optional_evolver.orbit(
            source, Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(optional_orbit.final().size() == 1u);
        const auto& optional_final = optional_orbit.final().front();
        ARIADNE_TESTING_PROVE(optional_final.location() == source_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(optional_final.previous_events().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            optional_final.continuous_set().number_of_parameters(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(optional_final.number_of_constraints(), 0u);
        const auto optional_input = optional_final.continuous_set().parameter_domain()[1u];
        ARIADNE_TESTING_PROVE(optional_input.lower_bound().raw().get_d() <= -1.0);
        ARIADNE_TESTING_PROVE(optional_input.upper_bound().raw().get_d() >= 1.0);
        const auto optional_bounds = optional_final.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(optional_bounds.lower_bound().raw().get_d() <= -1.0);
        ARIADNE_TESTING_PROVE(optional_bounds.upper_bound().raw().get_d() >= 3.0);
    }

    Void test_lower_event_free_equation_witness() {
        const RealVariable x("x"), a("a");
        const DiscreteLocation mode(1);
        DiscreteTimeHybridAutomaton system;
        system.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1},
            List<RealAssignment>{let(a) = 2 * x}));
        Evolver evolver(system);
        const auto source = initial(mode);
        // A nonempty unconstrained ExactBox is a constructive existential
        // certificate, even though Enclosure::is_empty() may be UNKNOWN.
        ARIADNE_TESTING_EXPECT_EQUAL(source.continuous_set().number_of_constraints(), 0u);
        ARIADNE_TESTING_PROVE(
            definitely(!source.continuous_set().parameter_domain().is_empty()));

        // LOWER preserves a feasible initial witness; it is not an
        // INNER inclusion of all states and must never subdivide.
        const auto zero = evolver.orbit(
            source, Termination(Integer(0), Integer(0)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.reach().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.initial().auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(source.auxiliary_space().dimension(), 0u);

        const auto lower = evolver.orbit(
            source, Termination(Integer(0), Integer(3)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.reach().size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.intermediate().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.final().size(), 1u);
        const auto& final = lower.final().front();
        ARIADNE_TESTING_PROVE(final.location() == mode);
        ARIADNE_TESTING_EXPECT_EQUAL(final.auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(final.time_range().lower_bound().raw().get_d(), 3.0);
        const auto bounds = final.state_bounding_box()[0u];
        // Initial x[0] in [0,1] and x[n]=x[0]+n: x[0]=1/2 is
        // a concrete witness with x[3]=3.5.
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= 3.5);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 3.5);
        ARIADNE_TESTING_EXPECT_EQUAL(source.time_range().lower_bound().raw().get_d(), 0.0);
        ARIADNE_TESTING_EXPECT_EQUAL(source.auxiliary_space().dimension(), 0u);

        // Absolute integer clock is part of the contract on resume.
        auto resumed = final;
        const auto resumed_orbit = evolver.orbit(
            resumed, Termination(Integer(3), Integer(4)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(resumed_orbit.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            resumed_orbit.final().front().time_range().lower_bound().raw().get_d(), 4.0);

        evolver.configuration().set_maximum_nodes(2u);
        ARIADNE_TESTING_EXPECT_THROWS(std::length_error,
            evolver.orbit(source, Termination(Integer(0), Integer(3)), Semantics::LOWER));
        const auto short_orbit = evolver.orbit(
            source, Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(short_orbit.final().size(), 1u);
    }

    Void test_lower_constrained_point_witness() {
        const RealVariable x("x");
        const DiscreteLocation mode(1);
        DiscreteTimeHybridAutomaton system;
        system.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        Evolver evolver(system);

        // Not all x in [0,1] satisfy these constraints, but x=1/2
        // does. The validated midpoint certificate must admit LOWER
        // without dropping the original constraints or choosing a
        // different existential semantics from the ordinary evolver.
        auto feasible = initial(mode);
        Enclosure constrained = feasible.continuous_set();
        const auto coordinate = EffectiveScalarMultivariateFunction::coordinate(1u, 0u);
        constrained.new_positive_state_constraint(coordinate - Real(1 / 4_q));
        constrained.new_negative_state_constraint(coordinate - Real(3 / 4_q));
        feasible.replace_continuous_set(constrained);
        ARIADNE_TESTING_EXPECT_EQUAL(feasible.number_of_constraints(), 2u);
        const auto zero = evolver.orbit(
            feasible, Termination(Integer(0), Integer(0)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.initial().number_of_constraints(), 2u);

        const auto reached = evolver.orbit(
            feasible, Termination(Integer(0), Integer(2)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(reached.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(reached.final().front().number_of_constraints(), 2u);
        const auto bounds = reached.final().front().state_bounding_box()[0u];
        // Starting at the certified x[0]=1/2 gives x[2]=2.5.
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= 2.5);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 2.5);
        ARIADNE_TESTING_EXPECT_EQUAL(feasible.number_of_constraints(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            feasible.time_range().lower_bound().raw().get_d(), 0.0);

        // An impossible constraint is not an admissible LOWER source.
        auto impossible = initial(mode);
        Enclosure contradiction = impossible.continuous_set();
        contradiction.new_positive_state_constraint(coordinate - Real(2));
        impossible.replace_continuous_set(contradiction);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(impossible,
                          Termination(Integer(0), Integer(0)), Semantics::LOWER));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(impossible,
                          Termination(Integer(0), Integer(1)), Semantics::LOWER));
    }

    Void test_lower_forced_urgent_jump_and_disabled_progression() {
        const RealVariable x("x");
        const DiscreteLocation source_mode(1), target_mode(2);
        const DiscreteEvent event("forced_urgent");
        DiscreteTimeHybridAutomaton system;
        system.new_mode(source_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        system.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        // A strictly positive guard is uniformly true on x[0] in [0,1].
        // This urgent reset is forced, so every initial witness follows it.
        system.new_urgent_transition(source_mode, event, target_mode,
                                     {next(x) = x + 3}, x >= -1);
        Evolver evolver(system);
        const auto source = initial(source_mode);
        const auto lower = evolver.orbit(source,
            Termination(Integer(0), Integer(2)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.reach().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.intermediate().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.final().size(), 1u);
        const auto& jumped = lower.reach().front();
        ARIADNE_TESTING_PROVE(jumped.location() == target_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(jumped.previous_events().size(), 1u);
        ARIADNE_TESTING_PROVE(jumped.previous_events()[0u] == event);
        const auto& final = lower.final().front();
        ARIADNE_TESTING_PROVE(final.location() == target_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(final.time_range().lower_bound().raw().get_d(), 2.0);
        const auto bounds = final.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= 4.5);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 4.5);
        ARIADNE_TESTING_PROVE(source.location() == source_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(source.previous_events().size(), 0u);

        // When the same urgent guard is uniformly false, only the
        // original mode's deterministic progression is possible.
        DiscreteTimeHybridAutomaton disabled;
        disabled.new_mode(source_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        disabled.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        disabled.new_urgent_transition(source_mode, event, target_mode,
                                       {next(x) = x + 3}, x >= 2);
        Evolver no_jump(disabled);
        const auto progressed = no_jump.orbit(source,
            Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(progressed.final().size(), 1u);
        const auto& survivor = progressed.final().front();
        ARIADNE_TESTING_PROVE(survivor.location() == source_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(survivor.previous_events().size(), 0u);
        const auto step_bounds = survivor.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(step_bounds.lower_bound().raw().get_d() <= 1.5);
        ARIADNE_TESTING_PROVE(step_bounds.upper_bound().raw().get_d() >= 1.5);
    }

    Void test_lower_rejects_unaggregated_guard_splitting() {
        const RealVariable x("x");
        const DiscreteLocation source_mode(1), target_mode(2);
        const DiscreteEvent event("undecided");
        DiscreteTimeHybridAutomaton system;
        system.new_mode(source_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        system.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        system.new_urgent_transition(source_mode, event, target_mode,
                                     {next(x) = x + 3}, x >= 1 / 2_q);
        Evolver evolver(system);
        const auto source = initial(source_mode);

        // An ambiguous whole source cannot be replaced by an
        // independently selected child of its outer LOWER enclosure.
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(source, Termination(Integer(0), Integer(1)),
                          Semantics::LOWER));

        // The caller provides a uniformly guard-positive source,
        // with reachable reset witness x[0]=7/8. No split is made.
        const HybridEnclosure reset_source(HybridExactBox(
            source_mode, RealSpace({x}), ExactBoxType({{0.75_x, 1.0_x}})),
            configuration());
        const auto lower = evolver.orbit(
            reset_source, Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.final().size(), 1u);
        ARIADNE_TESTING_PROVE(lower.final().front().location() == target_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.final().front().previous_events().size(), 1u);
        // Verify the forced reset by both event-name and identity checks.
        const auto& reset_final = lower.final().front();
        const auto& recorded_events = reset_final.previous_events();
        ARIADNE_TESTING_REQUIRE(recorded_events.size() == 1u);
        ARIADNE_TESTING_REQUIRE(recorded_events[0u].name() == event.name());
        ARIADNE_TESTING_PROVE(lower.final().front().previous_events()[0u] == event);
        // The complete LOWER reset branch must still own a
        // one-dimensional validated reset map and a nonempty parameter
        // domain. Do not bypass its public bounding_box() operation.
        ARIADNE_TESTING_REQUIRE(reset_final.state_dimension() == 1u);
        ARIADNE_TESTING_REQUIRE(
            reset_final.continuous_set().state_function().result_size() == 1u);
        ARIADNE_TESTING_REQUIRE(
            definitely(!reset_final.continuous_set().reduced_domain().is_empty()));
        const auto reset_box = reset_final.state_bounding_box();
        ARIADNE_TESTING_REQUIRE(reset_box.dimension() == 1u);
        const auto bounds = reset_box[0u];
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= 3.875);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 3.875);
        ARIADNE_TESTING_EXPECT_EQUAL(source.previous_events().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(source.continuous_set().number_of_parameters(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(reset_source.previous_events().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(reset_source.number_of_constraints(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            evolver.orbit(source, Termination(Integer(0), Integer(0)),
                          Semantics::LOWER).final().size(), 1u);

        // The caller-supplied guard-negative source preserves
        // the event-free progression witness
        // x[0]=1/8 -> x[1]=9/8.
        const HybridEnclosure progression_source(HybridExactBox(
            source_mode, RealSpace({x}), ExactBoxType({{0.0_x, 0.25_x}})),
            configuration());
        const auto progression = evolver.orbit(
            progression_source, Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(progression.final().size() == 1u);
        const auto& progression_final = progression.final().front();
        ARIADNE_TESTING_PROVE(progression_final.location() == source_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(progression_final.previous_events().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(progression_final.number_of_constraints(), 0u);
        const auto progression_box = progression_final.state_bounding_box();
        ARIADNE_TESTING_REQUIRE(progression_box.dimension() == 1u);
        ARIADNE_TESTING_PROVE(
            progression_box[0u].lower_bound().raw().get_d() <= 1.125);
        ARIADNE_TESTING_PROVE(
            progression_box[0u].upper_bound().raw().get_d() >= 1.125);
        ARIADNE_TESTING_EXPECT_EQUAL(progression_source.number_of_constraints(), 0u);

        // Keep the caller's constraint unchanged; do not add a
        // synthetic subbox restriction to resolve the event.
        auto constrained = reset_source;
        Enclosure restricted = constrained.continuous_set();
        const auto coordinate = EffectiveScalarMultivariateFunction::coordinate(1u, 0u);
        restricted.new_positive_state_constraint(coordinate - Real(3 / 4_q));
        constrained.replace_continuous_set(restricted);
        const auto feasible = evolver.orbit(
            constrained, Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(feasible.final().size(), 1u);
        ARIADNE_TESTING_PROVE(feasible.final().front().location() == target_mode);
        // Only the caller's original x>=3/4 constraint survives.
        ARIADNE_TESTING_EXPECT_EQUAL(
            feasible.final().front().number_of_constraints(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(constrained.number_of_constraints(), 1u);
        const auto feasible_box = feasible.final().front().state_bounding_box();
        ARIADNE_TESTING_REQUIRE(feasible_box.dimension() == 1u);
        ARIADNE_TESTING_PROVE(feasible_box[0u].lower_bound().raw().get_d() <= 3.875);
        ARIADNE_TESTING_PROVE(feasible_box[0u].upper_bound().raw().get_d() >= 3.875);

        // A guard-boundary singleton cannot be treated as a
        // uniformly decided transition without inventing evidence.
        const HybridEnclosure boundary(HybridExactBox(
            source_mode, RealSpace({x}),
            ExactBoxType({{0.5_x, 0.5_x}})), configuration());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(boundary, Termination(Integer(0), Integer(1)),
                          Semantics::LOWER));
    }

    Void test_lower_preserves_unsplit_multidimensional_source() {
        const RealVariable x("x"), y("y");
        const DiscreteLocation source_mode(1), target_mode(2);
        const DiscreteEvent event("guard_on_y");
        DiscreteTimeHybridAutomaton system;
        system.new_mode(source_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1, next(y) = y}));
        system.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x, next(y) = y}));
        system.new_urgent_transition(source_mode, event, target_mode,
            {next(x) = x + 3, next(y) = y}, y >= 1 / 2_q);
        Evolver evolver(system);

        // A second coordinate with an undecided guard cannot be
        // isolated as an independent LOWER child.
        const HybridEnclosure ambiguous(HybridExactBox(
            source_mode, RealSpace({x, y}),
            ExactBoxType({{0.0_x, 1.0_x}, {0.0_x, 1.0_x}})), configuration());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(ambiguous, Termination(Integer(0), Integer(1)),
                          Semantics::LOWER));
        // Uniformly decided whole sources remain supported.
        const HybridEnclosure reset_source(HybridExactBox(
            source_mode, RealSpace({x, y}),
            ExactBoxType({{0.0_x, 1.0_x}, {0.75_x, 1.0_x}})), configuration());
        const auto reset_orbit = evolver.orbit(
            reset_source, Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(reset_orbit.final().size() == 1u);
        const auto& reset_final = reset_orbit.final().front();
        ARIADNE_TESTING_PROVE(reset_final.location() == target_mode);
        ARIADNE_TESTING_REQUIRE(reset_final.previous_events().size() == 1u);
        ARIADNE_TESTING_PROVE(reset_final.previous_events()[0u] == event);
        ARIADNE_TESTING_EXPECT_EQUAL(reset_final.number_of_constraints(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(reset_final.state_dimension(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            reset_final.time_range().lower_bound().raw().get_d(), 1.0);
        const auto reset_bounds = reset_final.state_bounding_box();
        ARIADNE_TESTING_REQUIRE(reset_bounds.dimension() == 2u);
        ARIADNE_TESTING_PROVE(reset_bounds[0u].lower_bound().raw().get_d() <= 3.5);
        ARIADNE_TESTING_PROVE(reset_bounds[0u].upper_bound().raw().get_d() >= 3.5);
        ARIADNE_TESTING_PROVE(reset_bounds[1u].lower_bound().raw().get_d() <= 0.875);
        ARIADNE_TESTING_PROVE(reset_bounds[1u].upper_bound().raw().get_d() >= 0.875);
        ARIADNE_TESTING_EXPECT_EQUAL(reset_source.number_of_constraints(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(reset_source.previous_events().size(), 0u);

        // On the other side, the only certifiable branch is progression:
        // (x,y)=(1/2,1/8) -> (3/2,1/8), with no recorded event.
        const HybridEnclosure progression_source(HybridExactBox(
            source_mode, RealSpace({x, y}),
            ExactBoxType({{0.0_x, 1.0_x}, {0.0_x, 0.25_x}})), configuration());
        const auto progression_orbit = evolver.orbit(
            progression_source, Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(progression_orbit.final().size() == 1u);
        const auto& progression_final = progression_orbit.final().front();
        ARIADNE_TESTING_PROVE(progression_final.location() == source_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(progression_final.previous_events().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(progression_final.number_of_constraints(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            progression_final.time_range().lower_bound().raw().get_d(), 1.0);
        const auto progression_bounds = progression_final.state_bounding_box();
        ARIADNE_TESTING_REQUIRE(progression_bounds.dimension() == 2u);
        ARIADNE_TESTING_PROVE(progression_bounds[0u].lower_bound().raw().get_d() <= 1.5);
        ARIADNE_TESTING_PROVE(progression_bounds[0u].upper_bound().raw().get_d() >= 1.5);
        ARIADNE_TESTING_PROVE(progression_bounds[1u].lower_bound().raw().get_d() <= 0.125);
        ARIADNE_TESTING_PROVE(progression_bounds[1u].upper_bound().raw().get_d() >= 0.125);
        ARIADNE_TESTING_EXPECT_EQUAL(progression_source.number_of_constraints(), 0u);

        // The boundary singleton is feasible but cannot force either
        // branch under strict three-valued predicate certification.
        const HybridEnclosure boundary(HybridExactBox(
            source_mode, RealSpace({x, y}),
            ExactBoxType({{0.0_x, 1.0_x}, {0.5_x, 0.5_x}})), configuration());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(boundary, Termination(Integer(0), Integer(1)),
                          Semantics::LOWER));
    }

    Void test_lower_unique_permissive_branch() {
        const RealVariable x("x");
        const DiscreteLocation source_mode(1), target_mode(2);
        const DiscreteEvent event("permissive_choice");
        DiscreteTimeHybridAutomaton system;
        system.new_mode(source_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        system.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));

        // Explicit invariant FALSE means progression is prohibited,
        // while guard TRUE forces the reset despite PERMISSIVE kind.
        system.new_permissive_transition(source_mode, event, target_mode,
            {next(x) = x + 3}, x >= -1, x <= -1);
        Evolver jump(system);
        const auto reset = jump.orbit(initial(source_mode),
            Termination(Integer(0), Integer(2)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(reset.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(reset.reach().size(), 2u);
        const auto& image = reset.reach().front();
        ARIADNE_TESTING_PROVE(image.location() == target_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(image.previous_events().size(), 1u);
        ARIADNE_TESTING_PROVE(image.previous_events()[0u] == event);
        const auto bounds = reset.final().front().state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= 4.5);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 4.5);

        // Guard FALSE with progression invariant TRUE yields only
        // the deterministic state update, with no false event.
        DiscreteTimeHybridAutomaton progress;
        progress.new_mode(source_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        progress.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        progress.new_permissive_transition(source_mode, event, target_mode,
            {next(x) = x + 3}, x >= 2, x <= 2);
        Evolver no_jump(progress);
        const auto progressed = no_jump.orbit(initial(source_mode),
            Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(progressed.final().size(), 1u);
        ARIADNE_TESTING_PROVE(progressed.final().front().location() == source_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(progressed.final().front().previous_events().size(), 0u);

        // With a PERMISSIVE reset and a certain progression invariant,
        // ignoring the optional jump is an admissible LOWER continuation.
        // Continuous hybrid LOWER already follows this exact rule.
        DiscreteTimeHybridAutomaton ambiguous;
        ambiguous.new_mode(source_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        ambiguous.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        ambiguous.new_permissive_transition(source_mode, event, target_mode,
            {next(x) = x + 3}, x >= -1, x <= 2);
        Evolver undecided(ambiguous);
        const auto optional = undecided.orbit(initial(source_mode),
            Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(optional.final().size() == 1u);
        const auto& optional_final = optional.final().front();
        ARIADNE_TESTING_PROVE(optional_final.location() == source_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(optional_final.previous_events().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(optional_final.number_of_constraints(), 0u);
        const auto optional_bounds = optional_final.state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(optional_bounds.lower_bound().raw().get_d() <= 1.0);
        ARIADNE_TESTING_PROVE(optional_bounds.upper_bound().raw().get_d() >= 2.0);

        // A varying optional guard does not prevent certain progression.
        // Even if the reset is only POSSIBLE, skipping it is permitted.
        DiscreteTimeHybridAutomaton varying;
        varying.new_mode(source_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        varying.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        varying.new_permissive_transition(source_mode, event, target_mode,
            {next(x) = x + 3}, x >= 1 / 2_q, x <= 2);
        Evolver varying_guard(varying);
        const auto optional_unknown = varying_guard.orbit(initial(source_mode),
            Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(optional_unknown.final().size() == 1u);
        ARIADNE_TESTING_PROVE(optional_unknown.final().front().location() == source_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(
            optional_unknown.final().front().previous_events().size(), 0u);
        const auto varying_bounds =
            optional_unknown.final().front().state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(varying_bounds.lower_bound().raw().get_d() <= 1.0);
        ARIADNE_TESTING_PROVE(varying_bounds.upper_bound().raw().get_d() >= 2.0);

        // When both the guard and the progression invariant are
        // unresolved, neither branch is uniformly known to be valid.
        DiscreteTimeHybridAutomaton uncertain_invariant;
        uncertain_invariant.new_mode(source_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        uncertain_invariant.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        uncertain_invariant.new_permissive_transition(source_mode, event, target_mode,
            {next(x) = x + 3}, x >= 1 / 2_q, x <= 1 / 2_q);
        Evolver blocked(uncertain_invariant);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            blocked.orbit(initial(source_mode),
                Termination(Integer(0), Integer(1)), Semantics::LOWER));
    }

    Void test_lower_unique_branch_across_competing_events() {
        const RealVariable x("x");
        const DiscreteLocation source_mode(1), target_mode(2);
        const DiscreteEvent disabled("disabled"), forced("forced");
        DiscreteTimeHybridAutomaton system;
        system.new_mode(source_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        system.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        system.new_urgent_transition(source_mode, disabled, target_mode,
            {next(x) = x + 100}, x >= 2);
        system.new_urgent_transition(source_mode, forced, target_mode,
            {next(x) = x + 3}, x >= -1);
        Evolver evolver(system);
        const auto lower = evolver.orbit(initial(source_mode),
            Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.final().size(), 1u);
        ARIADNE_TESTING_PROVE(lower.final().front().location() == target_mode);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.final().front().previous_events().size(), 1u);
        ARIADNE_TESTING_PROVE(lower.final().front().previous_events()[0u] == forced);

        // Two uniformly true urgent guards create two alternatives.
        // LOWER must not use event enumeration order as priority.
        DiscreteTimeHybridAutomaton competing;
        competing.new_mode(source_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        competing.new_mode(target_mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        competing.new_urgent_transition(source_mode, disabled, target_mode,
            {next(x) = x + 100}, x >= -1);
        competing.new_urgent_transition(source_mode, forced, target_mode,
            {next(x) = x + 3}, x >= -1);
        Evolver not_unique(competing);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            not_unique.orbit(initial(source_mode),
                Termination(Integer(0), Integer(1)), Semantics::LOWER));
    }

    Void test_lower_excludes_unproved_branches() {
        const RealVariable x("x"), u("u");
        const DiscreteLocation mode(1), target(2);
        const DiscreteEvent event("conditional_jump");
        DiscreteTimeHybridAutomaton with_event;
        with_event.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        with_event.new_mode(target, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x}));
        with_event.new_permissive_transition(mode, event, target,
            {next(x) = x + 1}, x >= 0, x <= 2);
        Evolver guarded(with_event);
        // Zero horizon requires an initial witness but no event decision.
        ARIADNE_TESTING_EXPECT_EQUAL(
            guarded.orbit(initial(mode), Termination(Integer(0), Integer(0)),
                          Semantics::LOWER).final().size(), 1u);
        // Certain progression with a PERMISSIVE alternative is a
        // legitimate full-source LOWER step, not an invented priority.
        const auto optional = guarded.orbit(
            initial(mode), Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(optional.final().size() == 1u);
        ARIADNE_TESTING_PROVE(optional.final().front().location() == mode);
        ARIADNE_TESTING_EXPECT_EQUAL(
            optional.final().front().previous_events().size(), 0u);
        const auto optional_bounds = optional.final().front().state_bounding_box()[0u];
        ARIADNE_TESTING_PROVE(optional_bounds.lower_bound().raw().get_d() <= 1.0);
        ARIADNE_TESTING_PROVE(optional_bounds.upper_bound().raw().get_d() >= 2.0);

        DiscreteTimeHybridAutomaton inclusion;
        inclusion.new_mode(mode, DifferenceAlgebraicInclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            RealVariablesBox({-1 <= u <= 1})));
        Evolver inclusion_evolver(inclusion);
        ARIADNE_TESTING_EXPECT_EQUAL(
            inclusion_evolver.orbit(initial(mode), Termination(Integer(0), Integer(0)),
                                   Semantics::LOWER).final().size(), 1u);
        const auto lower = inclusion_evolver.orbit(
            initial(mode), Termination(Integer(0), Integer(1)), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(lower.final().size() == 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            lower.final().front().continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(lower.final().front().number_of_constraints(), 0u);

        DiscreteTimeHybridAutomaton deterministic;
        deterministic.new_mode(mode, DifferenceAlgebraicEquation(
            List<PrimedRealAssignment>{next(x) = x + 1}));
        Evolver no_jump(deterministic);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            no_jump.orbit(initial(mode), Termination(Integer(1), Integer(1)),
                          Semantics::LOWER));
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
        ARIADNE_TESTING_EXPECT_EQUAL(
            evolver.orbit(initial(mode), Termination(Integer(0), Integer(0)),
                          Semantics::LOWER).final().size(), 1u);
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
        Evolver validated(inclusion);
        const auto source = initial(mode);
        const auto tree = validated.explore(
            source, Termination(Integer(0), Integer(1)), Semantics::UPPER);
        // a=x+u spans negative and positive values. The urgent guard
        // must retain BOTH the reset and the complementary progression.
        ARIADNE_TESTING_EXPECT_EQUAL(tree.nodes.size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(tree.nodes[0u].enclosure.auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(
            tree.nodes[0u].enclosure.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(tree.nodes[1u].enclosure.number_of_constraints(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(tree.nodes[2u].enclosure.number_of_constraints(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(tree.nodes[1u].enclosure.continuous_set().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(tree.nodes[2u].enclosure.continuous_set().number_of_parameters(), 3u);
        const auto orbit = validated.orbit(
            source, Termination(Integer(0), Integer(1)), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(source.auxiliary_space().dimension(), 0u);
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
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_public_reset_into_inclusion_with_auxiliary_and_continue);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_lower_inclusion_preserves_full_input_without_splitting);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_lower_uniform_auxiliary_guard_over_full_input_box);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_lower_event_free_equation_witness);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_lower_constrained_point_witness);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_lower_forced_urgent_jump_and_disabled_progression);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_lower_rejects_unaggregated_guard_splitting);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_lower_preserves_unsplit_multidimensional_source);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_lower_unique_permissive_branch);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_lower_unique_branch_across_competing_events);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_lower_excludes_unproved_branches);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeEvolver, test_rejections);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
