/***************************************************************************
 *            test_hybrid_continuous_time_approximate_evolver.cpp
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

/** @file tests/test_hybrid_continuous_time_approximate_evolver.cpp
 *  @brief Test hybrid evolver behavior and numerical invariants.
 */

#include <cmath>
#include <sstream>

#include "evolution/hybrid_continuous_time_approximate_evolver.hpp"
#include "evolution/hybrid_orbit.hpp"
#include "evolution/hybrid_time.hpp"
#include "log_verbosity_scope.hpp"
#include "logging/test_runner.hpp"
#include "system/hybrid_automata.hpp"

#include "ariadne-testing.hpp"
using namespace Ariadne;

// A valid progress predicate supplied through the supported automaton interface.
// ContinuousTimeHybridAutomaton's convenience invariant builder currently selects INVARIANT.
class ProgressAutomaton : public ContinuousTimeHybridAutomaton {
  public:
    ProgressAutomaton* clone() const override { return new ProgressAutomaton(*this); }
    EventKind event_kind(DiscreteLocation location, DiscreteEvent event) const override {
        auto kind = ContinuousTimeHybridAutomaton::event_kind(location, event);
        return kind == EventKind::INVARIANT ? EventKind::PROGRESS : kind;
    }
};

class TestHybridContinuousTimeApproximateEvolver : public Ariadne::Testing::ContextFixture {
  public:
    using Ariadne::Testing::ContextFixture::ContextFixture;

    RealVariable x{"x"}, y{"y"}, z{"z"};
    DiscreteLocation q{1}, target{2};
    DiscreteEvent jump{"jump"};
    using PointType = HybridContinuousTimeEvolver<ApproximateTag>::HybridApproximatePointType;
    using PointList = HybridContinuousTimeEvolver<ApproximateTag>::HybridApproximateListPointType;
    PointType point(double xv = 0.0) const {
        return PointType(q, RealSpace({x, y}), FloatDPApproximationPoint({ExactDouble(xv), 2.0_x}, dp));
    }
    Void expect_endpoint(const HybridInterpolatedCurve& curve, double time, double xv, double yv, double zv) {
        const auto& endpoint = curve.euclidean_set().end();
        ARIADNE_TESTING_PROVE(std::abs(endpoint->first.get_d() - time) < 1e-12);
        const auto& values = endpoint->second;
        ARIADNE_TESTING_PROVE(std::abs(values[curve.space()[x]].raw().get_d() - xv) < 1e-12);
        ARIADNE_TESTING_PROVE(std::abs(values[curve.space()[y]].raw().get_d() - yv) < 1e-12);
        ARIADNE_TESTING_PROVE(std::abs(values[curve.space()[z]].raw().get_d() - zv) < 1e-12);
    }
    ContinuousTimeHybridAutomaton flow_system() const {
        ContinuousTimeHybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0}, {let(z) = x + 2 * y});
        return system;
    }

  public:
    Void test_labelled_flow() {
        HybridContinuousTimeEvolver<ApproximateTag> evolver(flow_system());
        evolver.configuration().set_step_size(0.25);
        PointType reordered(q, RealSpace({y, x}), FloatDPApproximationPoint({2.0_x, 0.0_x}, dp));
        HybridRealPoint real(q, RealSpace({x, y}), RealPoint({Real(0), Real(2)}));
        auto single = evolver.orbit(point(), HybridTime(1.0_x, 3));
        auto permuted = evolver.orbit(reordered, HybridTime(1.0_x, 3));
        auto from_real = evolver.orbit(real, HybridTime(1.0_x, 3));
        expect_endpoint(single.curves()[0], 1.0, 1.0, 2.0, 5.0);
        expect_endpoint(permuted.curves()[0], 1.0, 1.0, 2.0, 5.0);
        expect_endpoint(from_real.curves()[0], 1.0, 1.0, 2.0, 5.0);
        auto list = evolver.orbit(PointList{point(), reordered}, HybridTime(1.0_x, 3));
        ARIADNE_TESTING_EXPECT_EQUAL(list.number_of_trajectories(), 2u);
        expect_endpoint(list.curves(0)[0], 1.0, 1.0, 2.0, 5.0);
        expect_endpoint(list.curves(1)[0], 1.0, 1.0, 2.0, 5.0);
        HybridTerminationCriterion no_phantom_event(Real(1), Integer(3), Set<DiscreteEvent>{DiscreteEvent()});
        auto fresh = evolver.orbit(PointList{point()}, no_phantom_event);
        expect_endpoint(fresh.curves(0u)[0], 1.0, 1.0, 2.0, 5.0);
        auto none = evolver.orbit(PointList{}, HybridTime(1.0_x, 3));
        ARIADNE_TESTING_EXPECT_EQUAL(none.number_of_trajectories(), 0u);
        auto zero = evolver.orbit(point(), HybridTime(0.0_x, 3));
        expect_endpoint(zero.curves()[0], 0.0, 0.0, 2.0, 4.0);
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.configuration().set_step_size(0.0));
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.configuration().set_step_size(-0.25));
        const auto& configured = evolver;
        ARIADNE_TESTING_EXPECT_EQUAL(configured.configuration().step_size().raw().get_d(), 0.25);
    }
    Void test_transitions_and_termination() {
        auto system = flow_system();
        system.new_mode(target, {dot(x) = 0, dot(y) = 0}, {let(z) = x + 2 * y});
        system.new_transition(q, jump, target, {next(x) = 0, next(y) = y}, x >= Real(1 / 4_q), EventKind::URGENT);
        HybridContinuousTimeEvolver<ApproximateTag> evolver(system);
        evolver.configuration().set_step_size(0.25);
        for (unsigned verbosity: {0u, 100u}) {
            LogVerbosityScope logging(verbosity);
            HybridTerminationCriterion stop(Real(1), Integer(3), Set<DiscreteEvent>{jump});
            auto single = evolver.orbit(point(), stop);
            ARIADNE_TESTING_EXPECT_EQUAL(single.size(), 2u);
            ARIADNE_TESTING_EXPECT_EQUAL(single.curves()[1].location(), target);
            expect_endpoint(single.curves()[0], 0.5, 0.5, 2.0, 4.5);
            expect_endpoint(single.curves()[1], 0.5, 0.0, 2.0, 4.0);
            auto multiple = evolver.orbit(PointList{point()}, stop);
            ARIADNE_TESTING_EXPECT_EQUAL(multiple.size(0u), 2u);
            expect_endpoint(multiple.curves(0u)[1], 0.5, 0.0, 2.0, 4.0);
            auto continued = evolver.orbit(point(), HybridTime(1.0_x, 3));
            expect_endpoint(continued.curves()[1], 1.0, 0.0, 2.0, 4.0);
            auto continued_list = evolver.orbit(PointList{point()}, HybridTime(1.0_x, 3));
            expect_endpoint(continued_list.curves(0u)[1], 1.0, 0.0, 2.0, 4.0);
        }
    }
    Void test_invariant_and_progress_stopping() {
        for (EventKind kind: {EventKind::INVARIANT, EventKind::PROGRESS}) {
            auto base = flow_system();
            DiscreteEvent block("block");
            base.new_invariant(q, x <= 1, block);
            ProgressAutomaton progress;
            progress.new_mode(q, {dot(x) = 1, dot(y) = 0}, {let(z) = x + 2 * y});
            progress.new_invariant(q, x <= 1, block);
            HybridContinuousTimeEvolver<ApproximateTag> evolver(kind == EventKind::INVARIANT
                                                                    ? static_cast<const ContinuousTimeHybridAutomatonInterface&>(base)
                                                                    : static_cast<const ContinuousTimeHybridAutomatonInterface&>(progress));
            for (unsigned verbosity: {0u, 100u}) {
                LogVerbosityScope logging(verbosity);
                auto stopped = evolver.orbit(point(2.0), HybridTime(1.0_x, 3));
                ARIADNE_TESTING_EXPECT_EQUAL(stopped.curve(0u).size(), 1u);
                expect_endpoint(stopped.curves()[0], 0.0, 2.0, 2.0, 6.0);
                auto stopped_list = evolver.orbit(PointList{point(2.0)}, HybridTime(1.0_x, 3));
                ARIADNE_TESTING_EXPECT_EQUAL(stopped_list.curve(0u, 0u).size(), 1u);
                auto valid = evolver.orbit(point(), HybridTime(0.5_x, 3));
                expect_endpoint(valid.curves()[0], 0.5, 0.5, 2.0, 4.5);
            }
        }
    }
    Void test_box_sampling() {
        auto system = flow_system();
        HybridContinuousTimeEvolver<ApproximateTag> evolver(system);
        HybridUpperBox box(q, RealSpace({
                                  x, y
        }),
                           UpperBoxType({{0, 1}, {1, 2}}));
        auto default_orbit = evolver.orbit(
            box, q, system.state_space(), HybridTime(0.0_x, 3));
        ARIADNE_TESTING_EXPECT_EQUAL(default_orbit.number_of_trajectories(), 1u);

        evolver.configuration().set_sampling_counts(SamplingCounts{
            {x, 3u},
            {y, 2u}
        });
        auto sampled = evolver.orbit(box, q, system.state_space(), HybridTime(0.0_x, 3));
        ARIADNE_TESTING_EXPECT_EQUAL(sampled.number_of_trajectories(), 6u);
        for (SizeType i = 0u; i < sampled.number_of_trajectories(); ++i) {
            const auto p = sampled.curve(static_cast<Nat>(i), 0u).begin()->second;
            ARIADNE_TESTING_PROVE(p[0].raw().get_d() >= 0.0 && p[0].raw().get_d() <= 1.0);
            ARIADNE_TESTING_PROVE(p[1].raw().get_d() >= 1.0 && p[1].raw().get_d() <= 2.0);
        }

        HybridUpperBox flat(q, RealSpace({
                                   x, y
        }),
                            UpperBoxType({{0, 0}, {1, 2}}));
        auto flat_orbit = evolver.orbit(flat, q, system.state_space(), HybridTime(0.0_x, 3));
        ARIADNE_TESTING_EXPECT_EQUAL(flat_orbit.number_of_trajectories(), 2u);
        evolver.configuration().clear_sampling_counts();
        auto singleton = evolver.orbit(
            HybridBoundedConstraintSet(q, {x == 0, y == 2}), HybridTime(0.0_x, 3));
        ARIADNE_TESTING_EXPECT_EQUAL(singleton.number_of_trajectories(), 1u);
        expect_endpoint(singleton.curves(0u)[0], 0.0, 0.0, 2.0, 4.0);
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.configuration().set_sampling_count(x, 0u));
    }
};

int main(int argc, char* argv[]) {
    ARIADNE_TESTING_SUITE("test_hybrid_continuous_time_approximate_evolver");
    ARIADNE_TESTING_ADD_TEST(TestHybridContinuousTimeApproximateEvolver, test_labelled_flow);
    ARIADNE_TESTING_ADD_TEST(TestHybridContinuousTimeApproximateEvolver, test_transitions_and_termination);
    ARIADNE_TESTING_ADD_TEST(TestHybridContinuousTimeApproximateEvolver, test_invariant_and_progress_stopping);
    ARIADNE_TESTING_ADD_TEST(TestHybridContinuousTimeApproximateEvolver, test_box_sampling);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
