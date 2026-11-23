/***************************************************************************
 *            test_hybrid_constant_parameters.cpp
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

/** @file tests/test_hybrid_constant_parameters.cpp
 *  @brief Verify hybrid constant parameters behavior.
 */

#include <stdexcept>

#include "logging/test_runner.hpp"

#include "ariadne-evolution.hpp"
#include "ariadne-testing.hpp"
using namespace Ariadne;

class TestHybridConstantParameters : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    const RealVariable x{"x"}, p{"p"}, z{"z"}, a{"a"}, u{"u"};
    const TimeVariable t;
    const DiscreteLocation source{1}, target{2}, middle{3};

    HybridEnclosure initial(const HybridContinuousTimeEvolver<ValidatedTag>& evolver) const {
        return evolver.enclosure(HybridExactBox(source, RealSpace({
                                                            x
        }),
                                                ExactBoxType({{0, 0}})));
    }
    void configure(HybridContinuousTimeEvolver<ValidatedTag>& evolver, bool tree) const {
        evolver.set_integrator(TaylorPicardIntegrator(step_maximum_error = 1e-7));
        evolver.configuration().set_maximum_step_size(0.03125);
        if (tree)
            evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
    }
    void zero(const ValidatedScalarMultivariateFunctionPatch& residual) {
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(residual.range()), ExactIntervalType(-0x1p-24_x, 0x1p-24_x)));
    }
    void binding(const HybridEnclosure& set) {
        const auto& continuous = set.continuous_set();
        ARIADNE_TESTING_EXPECT_EQUAL(continuous.constant_parameter_space(), RealSpace({p}));
        ARIADNE_TESTING_EXPECT_EQUAL(continuous.constant_parameter_count(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(set.state_dimension(), 1u);
        ARIADNE_TESTING_PROVE(definitely(continuous.constant_parameter_function(1u).get(0u).error() <= 0));
    }
    void test_equation_flow_and_continuation() {
        for (bool timed: {false, true})
            for (bool tree: {false, true}) {
                ContinuousTimeHybridAutomaton system;
                if (timed)
                    system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = p + t}, RealVariablesBox({2 <= p <= 3}), t, {let(a) = x + p + t}));
                else
                    system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = p}, RealVariablesBox({2 <= p <= 3}), std::nullopt, {let(a) = x + p}));
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                configure(evolver, tree);
                auto current = initial(evolver);
                const auto count = current.continuous_set().number_of_parameters();
                const auto held = evolver.orbit(current, HybridTime(0u, 1), Semantics::UPPER);
                ARIADNE_TESTING_REQUIRE(held.final().size() == 1u);
                ARIADNE_TESTING_EXPECT_EQUAL(held.final()[0u].continuous_set().number_of_parameters(), count);
                for (Nat step = 1u; step != 3u; ++step) {
                    const auto orbit = evolver.orbit(current, HybridTime(Real(step) / 8, 1), Semantics::UPPER);
                    ARIADNE_TESTING_REQUIRE(orbit.final().size() == 1u);
                    current = orbit.final()[0u];
                    binding(current);
                    const auto physical = current.continuous_set().constant_parameter_function(1u).get(0u);
                    const Real duration = Real(step) / 8;
                    zero(current.state_function()[0u] - physical * duration - (timed ? duration * duration / 2 : Real(0)));
                    zero(current.auxiliary_function()[0u] - current.state_function()[0u] - physical - (timed ? current.time_function() : 0 * current.time_function()));
                }
                ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, evolver.orbit(initial(evolver), HybridTime(0u, 1), Semantics::LOWER));
                evolver.set_integrator(PreconditionedGradedTaylorSeriesIntegrator(step_maximum_error = 1e-6));
                ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, evolver.orbit(initial(evolver), HybridTime(0u, 1), Semantics::UPPER));
            }
    }
    void test_inclusion_flow_and_method_gates() {
        for (bool timed: {false, true})
            for (bool boxes: {false, true})
                for (bool tree: {false, true}) {
                    ContinuousTimeHybridAutomaton system;
                    const RealVariablesBox inputs({Real(-0.125_x) <= u <= Real(0.125_x)}), parameters({2 <= p <= 3});
                    if (timed)
                        system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion({dot(x) = p + u + t}, inputs, parameters, t));
                    else
                        system.new_mode(source, OrdinaryDifferentialAlgebraicInclusion({dot(x) = p + u}, inputs, parameters));
                    HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                    configure(evolver, tree);
                    const auto seed = initial(evolver);
                    ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, evolver.orbit(seed, HybridTime(0u, 1), Semantics::UPPER));
                    if (boxes)
                        evolver.set_integrator(FirstOrderMultifunctionIntegrator());
                    else
                        evolver.set_integrator(TaylorPicardIntegrator(step_maximum_error = 1e-7), ZeroApproximation(), ConstantApproximation());
                    const auto orbit = evolver.orbit(seed, HybridTime(0.125_x, 1), Semantics::UPPER);
                    ARIADNE_TESTING_REQUIRE(orbit.final().size() == 1u);
                    const auto& final = orbit.final()[0u];
                    binding(final);
                    const ExactIntervalType expected = timed ? ExactIntervalType(0.2421875_x, 0.3984375_x) : ExactIntervalType(0.234375_x, 0.390625_x);
                    ARIADNE_TESTING_PROVE(subset(expected, cast_exact_interval(final.state_function()[0u].range())));
                }
    }
    void test_immediate_parameter_guard_reset_and_auxiliary() {
        for (auto kind: {EventKind::URGENT, EventKind::PERMISSIVE})
            for (bool tree: {false, true}) {
                ContinuousTimeHybridAutomaton system;
                system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, RealVariablesBox({2 <= p <= 3})));
                system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, RealVariablesBox({2 <= p <= 3}), std::nullopt, {let(a) = x - 2 * p}));
                system.new_transition(source, DiscreteEvent("jump"), target, {next(x) = 2 * p}, p >= Real(2.5_x), kind);
                system.new_invariant(source, p <= 3, DiscreteEvent("bound"));
                system.check_mode(source);
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                configure(evolver, tree);
                const auto orbit = evolver.orbit(initial(evolver), HybridTime(0.03125_x, 1), Semantics::UPPER);
                bool reached = false;
                for (const auto& final: orbit.final())
                    if (final.location() == target) {
                        reached = true;
                        binding(final);
                        zero(final.state_function()[0u] - 2 * final.continuous_set().constant_parameter_function(1u).get(0u));
                        zero(final.auxiliary_function()[0u]);
                    }
                ARIADNE_TESTING_PROVE(reached);
            }
    }
    void test_timed_parameter_guard_and_reset() {
        for (auto kind: {EventKind::URGENT, EventKind::IMPACT})
            for (bool tree: {false, true}) {
                ContinuousTimeHybridAutomaton system;
                system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = p}, RealVariablesBox({2 <= p <= 3})));
                system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, RealVariablesBox({2 <= p <= 3}), std::nullopt, {let(a) = x + 2 * p + t}));
                system.new_transition(source, DiscreteEvent("timer"), target, {next(x) = x + p + t}, t >= p / 32, kind);
                system.new_invariant(source, t <= Real(0.125_x) + p, DiscreteEvent("bound"));
                HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
                configure(evolver, tree);
                const auto orbit = evolver.orbit(initial(evolver), HybridTime(0.125_x, 1), Semantics::UPPER);
                bool reached = false;
                for (const auto& final: orbit.final())
                    if (final.location() == target) {
                        reached = true;
                        binding(final);
                        const auto physical = final.continuous_set().constant_parameter_function(1u).get(0u);
                        // Reset occurs on a conservative UPPER candidate. Its clock
                        // remains correlated with state and parameters, without a first-entry claim.
                        zero(final.state_function()[0u] - physical * final.time_function() - physical - final.time_function());
                        zero(final.dwell_time_function());
                        ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().size(), 1u);
                        if (kind == EventKind::IMPACT) {
                            bool boundary = false;
                            for (const auto& constraint: final.continuous_set().constraint_models()) {
                                if (definitely(constraint.lower_bound() == ExactNumber(0) && constraint.upper_bound() == ExactNumber(0))) {
                                    const auto residual = constraint.function() - final.time_function() + physical / 32;
                                    if (definitely(subset(cast_exact_interval(residual.range()), ExactIntervalType(-0x1p-24_x, 0x1p-24_x))))
                                        boundary = true;
                                }
                            }
                            ARIADNE_TESTING_PROVE(boundary);
                        }
                        zero(final.auxiliary_function()[0u] - final.state_function()[0u] - 2 * physical - final.time_function());
                    }
                ARIADNE_TESTING_PROVE(reached);
            }
    }
    void test_target_binding_and_parameter_free_intermediate_mode() {
        for (bool tree: {false, true}) {
            ContinuousTimeHybridAutomaton introduced;
            introduced.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
            introduced.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, RealVariablesBox({2 <= p <= 3}), std::nullopt, {let(a) = p}));
            introduced.new_transition(source, DiscreteEvent("bind"), target, {next(x) = 0}, x >= -1, EventKind::URGENT);
            HybridContinuousTimeEvolver<ValidatedTag> first(introduced);
            configure(first, tree);
            const auto bound = first.orbit(initial(first), HybridTime(0.03125_x, 1), Semantics::UPPER);
            ARIADNE_TESTING_REQUIRE(bound.final().size() == 1u);
            binding(bound.final()[0u]);
            zero(bound.final()[0u].auxiliary_function()[0u] - bound.final()[0u].continuous_set().constant_parameter_function(1u).get(0u));

            ContinuousTimeHybridAutomaton carried;
            carried.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, RealVariablesBox({2 <= p <= 3})));
            carried.new_mode(middle, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, {let(a) = x}));
            carried.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, RealVariablesBox({2 <= p <= 3}), std::nullopt, {let(a) = x - p}));
            carried.new_transition(source, DiscreteEvent("out"), middle, {next(x) = p}, x >= -1, EventKind::URGENT);
            carried.new_transition(middle, DiscreteEvent("back"), target, {next(x) = x}, x >= Real(2.5_x), EventKind::PERMISSIVE);
            HybridContinuousTimeEvolver<ValidatedTag> second(carried);
            configure(second, tree);
            const auto orbit = second.orbit(initial(second), HybridTime(0.03125_x, 2), Semantics::UPPER);
            bool reached = false;
            for (const auto& final: orbit.final()) {
                binding(final);
                zero(final.state_function()[0u] - final.continuous_set().constant_parameter_function(1u).get(0u));
                if (final.location() == target) {
                    reached = true;
                    zero(final.auxiliary_function()[0u]);
                }
            }
            ARIADNE_TESTING_PROVE(reached);
        }
    }
    void test_incompatible_declarations_and_reset_rejection() {
        for (bool renamed: {false, true}) {
            ContinuousTimeHybridAutomaton system;
            system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, RealVariablesBox({2 <= p <= 3})));
            system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, renamed ? RealVariablesBox({2 <= z <= 3}) : RealVariablesBox({Real(2.5_x) <= p <= 3})));
            system.new_transition(source, DiscreteEvent("jump"), target, {next(x) = x}, x >= -1, EventKind::URGENT);
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            configure(evolver, false);
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, evolver.orbit(initial(evolver), HybridTime(0.03125_x, 1), Semantics::UPPER));
        }
        ContinuousTimeHybridAutomaton invalid;
        invalid.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, RealVariablesBox({2 <= p <= 3})));
        invalid.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, RealVariablesBox({2 <= p <= 3})));
        ARIADNE_TESTING_EXPECT_EXCEPTION(([&] {
            invalid.new_transition(source, DiscreteEvent("change_parameter"), target, {next(x) = x, next(p) = 0}, x >= -1, EventKind::URGENT);
            invalid.check_mode(source);
        }()));
    }
    void test_grouping_distinguishes_physical_tuples() {
        ContinuousTimeHybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        configure(evolver, false);
        const auto seed = initial(evolver);
        HybridEvolutionTree tree;
        auto add = [&](const RealVariable& variable) {
            auto set = seed;
            auto continuous = set.continuous_set();
            continuous.adjoin_constant_parameters(RealSpace({
                                                      variable
            }),
                                                  ExactBoxType({{2, 3}}));
            set.replace_continuous_set(continuous);
            HybridEvolutionTrajectory trajectory;
            trajectory.entries.push_back({
                {set, HybridEvolutionPhase::PROCESS_STARTING_EVENTS},
                std::nullopt,
                std::nullopt
            });
            const auto node = tree.empty() ? tree.create_root(std::move(trajectory), HybridEvolutionTrajectoryMetadata{}) : tree.add_child(tree.root().value(), std::move(trajectory), HybridEvolutionTrajectoryMetadata{});
            return HybridEvolutionCheckpointId{node, 0u};
        };
        const auto first = add(p), same = add(p), other = add(z);
        const AllCompatibleRecombinationPolicy policy;
        ARIADNE_TESTING_EXPECT_EQUAL(policy.select(tree, {first, same}).size(), 1u);
        ARIADNE_TESTING_PROVE(policy.select(tree, {first, other}).empty());
    }
};
int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_constant_parameters");
    ARIADNE_TESTING_ADD_TEST(TestHybridConstantParameters, test_equation_flow_and_continuation);
    ARIADNE_TESTING_ADD_TEST(TestHybridConstantParameters, test_inclusion_flow_and_method_gates);
    ARIADNE_TESTING_ADD_TEST(TestHybridConstantParameters, test_immediate_parameter_guard_reset_and_auxiliary);
    ARIADNE_TESTING_ADD_TEST(TestHybridConstantParameters, test_timed_parameter_guard_and_reset);
    ARIADNE_TESTING_ADD_TEST(TestHybridConstantParameters, test_target_binding_and_parameter_free_intermediate_mode);
    ARIADNE_TESTING_ADD_TEST(TestHybridConstantParameters, test_incompatible_declarations_and_reset_rejection);
    ARIADNE_TESTING_ADD_TEST(TestHybridConstantParameters, test_grouping_distinguishes_physical_tuples);
    return Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
