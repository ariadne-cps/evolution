/***************************************************************************
 *            test_differential_inclusion_evolver.cpp
 *
 *  Copyright  2008-21  Luca Geretti, Pieter Collins, Sanja Zivanovic
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

/** @file tests/test_differential_inclusion_evolver.cpp
 *  @brief Test differential inclusion evolver behavior and numerical invariants.
 */

#include <limits>
#include <sstream>

#include "../src/inclusion_evolver_state.hpp"
#include "algebra/algebra.hpp"
#include "algebra/sweeper.hpp"
#include "evolution/differential_inclusion_evolver.hpp"
#include "evolution/orbit.hpp"
#include "function/formula.hpp"
#include "function/function.hpp"
#include "function/taylor_model.hpp"
#include "geometry/box.hpp"
#include "integration/integrator.hpp"
#include "labelled/figure.hpp"
#include "log_verbosity_scope.hpp"
#include "set/function_set.hpp"
#include "test_logging.hpp"

#include "ariadne-testing.hpp"
#include "expression-set/expression_set.hpp"

using namespace Ariadne;

class TestDifferentialInclusionEvolver : public Ariadne::Testing::ContextFixture {
  public:
    void test_scheduler_preserves_rigorous_reach_when_affine_wins() {
        RealVariable x("x"), u("u");
        DifferentialInclusion system({dot(x) = u - x}, RealVariablesBox({-1 <= u <= 1}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-5);
        InclusionEvolverState state(system, {ZeroApproximation(), ConstantApproximation(), AffineApproximation()}, integrator);
        auto available = state.approximators_to_use();
        auto best = available[0];
        for (const auto& approximation: available)
            if (approximation.index() == AffineApproximation().index())
                best = approximation;
        ARIADNE_TESTING_PROVE(not best.is_reach_rigorous());
        state.update_with_best(best);
        for (const auto& approximation: available) {
            Nat expected = approximation.index() == ZeroApproximation().index() ? 1u : 0u;
            ARIADNE_TESTING_EXPECT_EQUAL(state.check_delay()[approximation], expected);
        }
        ARIADNE_TESTING_EXPECT_EQUAL(state.global_optima_count()[best], 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(state.local_optima_count()[best], 1u);
        state.reset_local_optima_count();
        ARIADNE_TESTING_EXPECT_EQUAL(state.local_optima_count()[best], 0u);
        state.next_step();
        ARIADNE_TESTING_EXPECT_EQUAL(state.step(), 1u);
        auto next = state.approximators_to_use();
        ARIADNE_TESTING_EXPECT_EQUAL(next.size(), 2u);
        for (const auto& approximation: next)
            ARIADNE_TESTING_PROVE(approximation.index() != ZeroApproximation().index());
        state.update_with_best(next[0]);
        state.next_step();
        ARIADNE_TESTING_EXPECT_EQUAL(state.step(), 2u);
        ARIADNE_TESTING_PROVE(state.approximators_to_use().size() >= 2u);
    }

    void test_scheduler_integer_boundaries() {
        const Nat maximum = std::numeric_limits<Nat>::max();
        const Nat bits = std::numeric_limits<Nat>::digits;
        ARIADNE_TESTING_EXPECT_EQUAL(EvolutionDetail::scheduled_step(0u, 0u), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(EvolutionDetail::scheduled_step(0u, bits - 1u), Nat(1u) << (bits - 1u));
        ARIADNE_TESTING_EXPECT_EQUAL(EvolutionDetail::scheduled_step(maximum - 1u, 0u), maximum);
        ARIADNE_TESTING_EXPECT_THROWS(std::overflow_error, EvolutionDetail::scheduled_step(0u, bits));
        ARIADNE_TESTING_EXPECT_THROWS(std::overflow_error, EvolutionDetail::scheduled_step(maximum, 0u));
    }
    using Ariadne::Testing::ContextFixture::ContextFixture;

    // One scalar affine inclusion: cheap even under MSVC Debug, but still
    // exercises uncertain inputs and periodic parameter reconditioning.
    void run_case(const List<InputApproximation>& approximations,
                  Nat reduction_period, ExactDouble ratio,
                  Real evolution_time, ExactDouble maximum_step, Bool combined = false) {
        RealVariable x("x"), u("u");
        DottedRealAssignments dynamics = {dot(x) = u - x};
        RealVariablesBox inputs = {0.5_dec <= u <= 1.5_dec};
        Real e = 1 / 32_q;
        RealVariablesBox initial = {-e <= x <= e};
        DifferentialInclusion inclusion(dynamics, inputs);

        ThresholdSweeperDP sweeper(DoublePrecision(), 1e-8_pr);
        TaylorPicardIntegrator integrator(
            step_maximum_error = 1e-5_pr,
            sweeper,
            lipschitz_tolerance = 0.5_x,
            minimum_temporal_order = 4,
            maximum_temporal_order = 8);
        LohnerReconditioner reconditioner(
            static_cast<Nat>(initial.variables().size()),
            static_cast<Nat>(inputs.variables().size()),
            reduction_period, ratio);

        for (auto approximation: approximations) {
            DifferentialInclusionEvolver evolver(inclusion, integrator, reconditioner);
            List<InputApproximation> selected = combined ? approximations : List<InputApproximation>{approximation};
            evolver.configuration().set_approximations(selected);
            evolver.configuration().set_maximum_step_size(maximum_step);
            auto orbit = evolver.orbit(initial, evolution_time);
            ARIADNE_TESTING_PROVE(orbit.reach().size() > 0u);
            ARIADNE_TESTING_PROVE(orbit.final().size() > 0u);
            if (combined)
                break;
        }
    }

  public:
    void test_uncertain_input() {
        // One step with each of the two input approximations.
        run_case({ConstantApproximation(), AffineApproximation()},
                 3, 3.0_x, 1 / 32_q, 1.0_x / pow(two, 5));
    }

    void test_periodic_reconditioning() {
        // Four short steps and parameter reduction every two steps.
        run_case({AffineApproximation()},
                 2, 1.25_x, 1 / 8_q, 1.0_x / pow(two, 5));
    }

    void test_scheduling_and_diagnostic_configuration() {
        LogVerbosityScope logging(100u);
        run_case({ZeroApproximation(), ConstantApproximation(), AffineApproximation(),
                  SinusoidalApproximation(), PiecewiseApproximation()},
                 2, 1.25_x, 1 / 4_q, 0.03125_x, true);
        run_case({AffineApproximation()}, 0, 3.0_x, 1 / 16_q, 0.03125_x);
        DifferentialInclusionEvolverConfiguration config;
        config.set_maximum_enclosure_radius(2.0_x);
        config.set_enable_parameter_reduction(false);
        ARIADNE_TESTING_EXPECT_EQUAL(config.maximum_enclosure_radius(), 2.0_x);
        ARIADNE_TESTING_PROVE(!config.enable_parameter_reduction());
        ARIADNE_TESTING_EXPECT_EXCEPTION(config.set_approximations({}));
        std::ostringstream text;
        text << config;
        ARIADNE_TESTING_PROVE(text.str().find("approximations=") != std::string::npos);
        LohnerReconditioner policy(1u, 1u, 2u, 1.25_x);
        std::unique_ptr<ReconditionerInterface> clone(policy.clone());
        ARIADNE_TESTING_PROVE(clone.get() != &policy);
    }

    void test_error_preserving_reconditioner_and_empty_inputs() {
        struct ErrorPreservingReconditioner : LohnerReconditioner {
            ErrorPreservingReconditioner() : LohnerReconditioner(1u, 1u, 0u, 1.0_x) {}
            Bool must_incorporate_errors(const InclusionEvolverState&) const override { return false; }
            ErrorPreservingReconditioner* clone() const override {
                return new ErrorPreservingReconditioner(*this);
            }
        } policy;
        RealVariable x("x"), u("u");
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-5);
        DifferentialInclusion system({dot(x) = u - x}, RealVariablesBox({0 <= u <= 1}));
        DifferentialInclusionEvolver evolver(system, integrator, policy);
        evolver.configuration().set_approximations({ConstantApproximation()});
        evolver.configuration().set_maximum_step_size(0.03125_x);
        const DifferentialInclusionEvolver& immutable = evolver;
        ARIADNE_TESTING_EXPECT_EQUAL(immutable.configuration().maximum_step_size(), 0.03125_x);
        auto result = evolver.orbit(RealVariablesBox({0 <= x <= 0}), Real(1 / 16_q));
        ARIADNE_TESTING_PROVE(result.final().size() > 0u);
        // The system rejects a zero-input domain before an evolver exists.
        ARIADNE_TESTING_EXPECT_EXCEPTION((DifferentialInclusion({dot(x) = -x}, RealVariablesBox{})));
    }
};

int main(int argc, const char* argv[]) {
    redirect_test_logging_to_null();
    ARIADNE_TESTING_SUITE("test_differential_inclusion_evolver");
    ARIADNE_TESTING_ADD_TEST(TestDifferentialInclusionEvolver, test_scheduler_preserves_rigorous_reach_when_affine_wins);
    ARIADNE_TESTING_ADD_TEST(TestDifferentialInclusionEvolver, test_scheduler_integer_boundaries);
    ARIADNE_TESTING_ADD_TEST(TestDifferentialInclusionEvolver, test_uncertain_input);
    ARIADNE_TESTING_ADD_TEST(TestDifferentialInclusionEvolver, test_periodic_reconditioning);
    ARIADNE_TESTING_ADD_TEST(TestDifferentialInclusionEvolver, test_scheduling_and_diagnostic_configuration);
    ARIADNE_TESTING_ADD_TEST(TestDifferentialInclusionEvolver, test_error_preserving_reconditioner_and_empty_inputs);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
