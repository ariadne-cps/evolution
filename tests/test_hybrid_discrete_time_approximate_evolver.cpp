/***************************************************************************
 *            test_hybrid_discrete_time_approximate_evolver.cpp
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

/** @file test_hybrid_discrete_time_approximate_evolver.cpp
 *  @brief Real branch histories for simultaneous discrete hybrid events.
 */

#include <memory>
#include <stdexcept>
#include <type_traits>

#include "logging/test_runner.hpp"
#include "system/discrete_time_hybrid_automaton.hpp"

#include "ariadne-evolution.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {
using Evolver = HybridDiscreteTimeEvolver<ApproximateTag>;
using Orbit = HybridDiscreteTimeApproximateOrbit;

static_assert(std::is_same_v<Evolver::ResultType, Orbit>);
static_assert(std::is_same_v<Evolver::Paradigm, ApproximateTag>);
static_assert(!std::is_same_v<Evolver::ResultType, HybridDiscreteTimeEvolver<ValidatedTag>::ResultType>);

class AlternatingInputSampling final : public InputSamplingPolicy {
  public:
    AlternatingInputSampling* clone() const override { return new AlternatingInputSampling(*this); }
    Vector<FloatDPApproximation> sample(const ExactBoxType&,
                                        const InputSamplingContext& ctx) const override {
        const int choice = ctx.step % 2u == 0u ? 1 : -1;
        return Vector<FloatDPApproximation>({FloatDPApproximation(choice, dp)});
    }
};

class TestHybridDiscreteTimeApproximateEvolver : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    RealVariable x{"x"};
    DiscreteLocation source{1}, first_target{2}, second_target{3};
    DiscreteEvent first{"first"}, second{"second"};

    DifferenceAlgebraicEquation increment() const {
        return DifferenceAlgebraicEquation(List<PrimedRealAssignment>{next(x) = x + 1});
    }

    DifferenceAlgebraicEquation hold() const {
        return DifferenceAlgebraicEquation(List<PrimedRealAssignment>{next(x) = x});
    }

    HybridApproximatePoint point(DiscreteLocation location, double value) const {
        return HybridApproximatePoint(location, RealSpace({x}),
                                      FloatDPApproximationPoint({ExactDouble(value)}, dp));
    }

    DiscreteTimeHybridAutomaton permissive() const {
        DiscreteTimeHybridAutomaton system;
        system.new_mode(source, increment());
        system.new_mode(first_target, hold());
        system.new_mode(second_target, hold());
        // The progression and BOTH transitions are enabled at x=1.
        system.new_permissive_transition(source, first, first_target,
                                         {next(x) = x + 10}, x <= 3 / 2_q, x <= 2);
        system.new_permissive_transition(source, second, second_target,
                                         {next(x) = x + 20}, x <= 3 / 2_q, x <= 2);
        return system;
    }

    Void test_multiple_activations_continue_as_separate_paths() {
        Evolver evolver(permissive());
        std::unique_ptr<Evolver> copy(evolver.clone());
        auto result = copy->orbit(point(source, 1.0), Integer(2));
        ARIADNE_TESTING_EXPECT_EQUAL(result.number_of_trajectories(), 3u);

        SizeType progress = 0u, first_branch = 0u, second_branch = 0u;
        for (const auto& path: result.trajectories()) {
            ARIADNE_TESTING_EXPECT_EQUAL(path.initial_index, 0u);
            ARIADNE_TESTING_EXPECT_EQUAL(path.points.size(), 3u);
            ARIADNE_TESTING_EXPECT_EQUAL(path.events.size(), 2u);
            ARIADNE_TESTING_PROVE(!path.blocked);
            ARIADNE_TESTING_EXPECT_EQUAL(path.points.front().location(), source);
            ARIADNE_TESTING_PROVE(!path.events[1u].has_value());
            if (!path.events[0u].has_value()) {
                ++progress;
                ARIADNE_TESTING_EXPECT_EQUAL(path.points[1u].location(), source);
                ARIADNE_TESTING_EXPECT_WITHIN(path.points.back().point()[0u].raw().get_d(), 3.0, 1e-14);
            } else if (*path.events[0u] == first) {
                ++first_branch;
                ARIADNE_TESTING_EXPECT_EQUAL(path.points[1u].location(), first_target);
                ARIADNE_TESTING_EXPECT_WITHIN(path.points.back().point()[0u].raw().get_d(), 11.0, 1e-14);
            } else if (*path.events[0u] == second) {
                ++second_branch;
                ARIADNE_TESTING_EXPECT_EQUAL(path.points[1u].location(), second_target);
                ARIADNE_TESTING_EXPECT_WITHIN(path.points.back().point()[0u].raw().get_d(), 21.0, 1e-14);
            } else {
                ARIADNE_TESTING_PROVE(false);
            }
        }
        ARIADNE_TESTING_EXPECT_EQUAL(progress, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(first_branch, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(second_branch, 1u);
    }

    Void test_urgent_blocks_progression_and_takes_one_tick() {
        DiscreteTimeHybridAutomaton system;
        system.new_mode(source, increment());
        system.new_mode(first_target, hold());
        system.new_mode(second_target, hold());
        system.new_urgent_transition(source, first, first_target,
                                     {next(x) = x + 10}, x >= 1 / 2_q);
        system.new_urgent_transition(source, second, second_target,
                                     {next(x) = x + 20}, x >= 1 / 2_q);
        Evolver evolver(system);
        auto paths = evolver.orbit(point(source, 1.0), Integer(1));
        ARIADNE_TESTING_EXPECT_EQUAL(paths.number_of_trajectories(), 2u);
        for (const auto& p: paths.trajectories()) {
            ARIADNE_TESTING_EXPECT_EQUAL(p.events.size(), 1u);
            ARIADNE_TESTING_PROVE(p.events[0u].has_value());
            ARIADNE_TESTING_PROVE(p.points[1u].location() != source);
        }

        auto below = evolver.orbit(point(source, -1.0), Integer(1));
        ARIADNE_TESTING_EXPECT_EQUAL(below.number_of_trajectories(), 1u);
        ARIADNE_TESTING_PROVE(!below.trajectory(0u).events[0u].has_value());
        ARIADNE_TESTING_EXPECT_EQUAL(below.trajectory(0u).points[1u].location(), source);
    }

    Void test_initial_bisection_and_explicit_path_budget() {
        Evolver evolver(permissive());
        auto initial_only = evolver.orbit(UpperBoxType({
                                              {0, 1}
        }),
                                          source, Integer(0));
        ARIADNE_TESTING_EXPECT_EQUAL(initial_only.number_of_trajectories(), 1u);
        ARIADNE_TESTING_EXPECT_WITHIN(
            initial_only.trajectory(0u).points[0u].point()[0u].raw().get_d(), 0.5, 1e-14);
        evolver.configuration().set_sampling_count(x, 3u);
        auto sampled = evolver.orbit(UpperBoxType({
                                         {0, 1}
        }),
                                     source, Integer(0));
        ARIADNE_TESTING_EXPECT_EQUAL(sampled.number_of_trajectories(), 3u);
        evolver.configuration().set_maximum_paths(2u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.orbit(UpperBoxType({
                                                           {0, 1}
        }),
                                                       source, Integer(0)));
        evolver.configuration().clear_sampling_counts();
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.orbit(point(source, 1.0), Integer(1)));
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.orbit(point(source, 1.0), Integer(-1)));
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.configuration().set_maximum_paths(0u));
    }

    Void test_inclusion_input_sampling_and_transition_from_prestate() {
        const RealVariable u("u");
        DiscreteTimeHybridAutomaton system;
        system.new_mode(source, DifferenceAlgebraicInclusion(
                                    List<PrimedRealAssignment>{next(x) = x + u},
                                    RealVariablesBox({-1 <= u <= 1})));
        system.new_mode(first_target, hold());
        system.new_permissive_transition(source, first, first_target,
                                         {next(x) = x + 10}, x <= 3 / 2_q, x <= 2);
        AlternatingInputSampling policy;
        Evolver evolver(system, policy);
        auto paths = evolver.orbit(point(source, 1.0), Integer(2));
        ARIADNE_TESTING_EXPECT_EQUAL(paths.number_of_trajectories(), 2u);
        SizeType progressed = 0u, reset = 0u;
        for (const auto& path: paths.trajectories()) {
            ARIADNE_TESTING_EXPECT_EQUAL(path.points.size(), 3u);
            if (!path.events[0u].has_value()) {
                ++progressed;
                // u[0]=+1, u[1]=-1. The value returns to the initial x.
                ARIADNE_TESTING_EXPECT_WITHIN(path.points.back().point()[0u].raw().get_d(), 1.0, 1e-14);
            } else {
                ++reset;
                // Reset acts on the source pre-state, not on x+u.
                ARIADNE_TESTING_EXPECT_WITHIN(path.points.back().point()[0u].raw().get_d(), 11.0, 1e-14);
            }
        }
        ARIADNE_TESTING_EXPECT_EQUAL(progressed, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(reset, 1u);

        Evolver missing_policy(system);
        ARIADNE_TESTING_EXPECT_EXCEPTION(missing_policy.orbit(point(source, 1.0), Integer(1)));
    }

    Void test_blocked_paths_are_not_erased() {
        DiscreteTimeHybridAutomaton system;
        system.new_mode(source, increment());
        system.new_mode(first_target, hold());
        system.new_permissive_transition(source, first, first_target,
                                         {next(x) = x + 10}, x >= 2, x <= 0);
        Evolver evolver(system);
        auto result = evolver.orbit(point(source, 1.0), Integer(3));
        ARIADNE_TESTING_EXPECT_EQUAL(result.number_of_trajectories(), 1u);
        ARIADNE_TESTING_PROVE(result.trajectory(0u).blocked);
        ARIADNE_TESTING_EXPECT_EQUAL(result.trajectory(0u).points.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.trajectory(0u).events.size(), 0u);
    }
};

} // namespace

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_discrete_time_approximate_evolver");
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeApproximateEvolver, test_multiple_activations_continue_as_separate_paths);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeApproximateEvolver, test_urgent_blocks_progression_and_takes_one_tick);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeApproximateEvolver, test_initial_bisection_and_explicit_path_budget);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeApproximateEvolver, test_inclusion_input_sampling_and_transition_from_prestate);
    ARIADNE_TESTING_ADD_TEST(TestHybridDiscreteTimeApproximateEvolver, test_blocked_paths_are_not_erased);
    return Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
