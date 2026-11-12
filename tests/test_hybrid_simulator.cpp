/***************************************************************************
 *            test_hybrid_simulator.cpp
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

/** @file tests/test_hybrid_simulator.cpp
 *  @brief Test hybrid simulator behavior and numerical invariants.
 */

#include <cmath>
#include <sstream>

#include "evolution/hybrid_orbit.hpp"
#include "evolution/hybrid_simulator.hpp"
#include "evolution/hybrid_time.hpp"
#include "log_verbosity_scope.hpp"
#include "system/hybrid_automata.hpp"
#include "test_logging.hpp"

#include "ariadne-testing.hpp"
using namespace Ariadne;

// A valid progress predicate supplied through the supported automaton interface.
// HybridAutomaton's convenience invariant builder currently selects INVARIANT.
class ProgressAutomaton : public HybridAutomaton {
  public:
    ProgressAutomaton* clone() const override { return new ProgressAutomaton(*this); }
    EventKind event_kind(DiscreteLocation location, DiscreteEvent event) const override {
        auto kind = HybridAutomaton::event_kind(location, event);
        return kind == EventKind::INVARIANT ? EventKind::PROGRESS : kind;
    }
};

class TestHybridSimulator : public Ariadne::Testing::ContextFixture {
  public:
    using Ariadne::Testing::ContextFixture::ContextFixture;

    RealVariable x{"x"}, y{"y"}, z{"z"};
    DiscreteLocation q{1}, target{2};
    DiscreteEvent jump{"jump"};
    using PointType = HybridSimulator::HybridApproximatePointType;
    using PointList = HybridSimulator::HybridApproximateListPointType;
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
    HybridAutomaton flow_system() const {
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0}, {let(z) = x + 2 * y});
        return system;
    }

  public:
    Void test_labelled_flow() {
        HybridSimulator simulator(flow_system());
        simulator.configuration().set_step_size(0.25);
        PointType reordered(q, RealSpace({y, x}), FloatDPApproximationPoint({2.0_x, 0.0_x}, dp));
        HybridRealPoint real(q, RealSpace({x, y}), RealPoint({Real(0), Real(2)}));
        auto single = simulator.orbit(point(), HybridTime(1.0_x, 3));
        auto permuted = simulator.orbit(reordered, HybridTime(1.0_x, 3));
        auto from_real = simulator.orbit(real, HybridTime(1.0_x, 3));
        expect_endpoint(single.curves()[0], 1.0, 1.0, 2.0, 5.0);
        expect_endpoint(permuted.curves()[0], 1.0, 1.0, 2.0, 5.0);
        expect_endpoint(from_real.curves()[0], 1.0, 1.0, 2.0, 5.0);
        auto list = simulator.orbit(PointList{point(), reordered}, HybridTime(1.0_x, 3));
        ARIADNE_TESTING_EXPECT_EQUAL(list.number_of_trajectories(), 2u);
        expect_endpoint(list.curves(0)[0], 1.0, 1.0, 2.0, 5.0);
        expect_endpoint(list.curves(1)[0], 1.0, 1.0, 2.0, 5.0);
        HybridTerminationCriterion no_phantom_event(Real(1), Integer(3), Set<DiscreteEvent>{DiscreteEvent()});
        auto fresh = simulator.orbit(PointList{point()}, no_phantom_event);
        expect_endpoint(fresh.curves(0u)[0], 1.0, 1.0, 2.0, 5.0);
        auto none = simulator.orbit(PointList{}, HybridTime(1.0_x, 3));
        ARIADNE_TESTING_EXPECT_EQUAL(none.number_of_trajectories(), 0u);
        auto zero = simulator.orbit(point(), HybridTime(0.0_x, 3));
        expect_endpoint(zero.curves()[0], 0.0, 0.0, 2.0, 4.0);
        ARIADNE_TESTING_EXPECT_EXCEPTION(simulator.configuration().set_step_size(0.0));
        ARIADNE_TESTING_EXPECT_EXCEPTION(simulator.configuration().set_step_size(-0.25));
        const auto& configured = simulator;
        ARIADNE_TESTING_EXPECT_EQUAL(configured.configuration().step_size().raw().get_d(), 0.25);
    }
    Void test_transitions_and_termination() {
        auto system = flow_system();
        system.new_mode(target, {dot(x) = 0, dot(y) = 0}, {let(z) = x + 2 * y});
        system.new_transition(q, jump, target, {next(x) = 0, next(y) = y}, x >= Real(1 / 4_q), EventKind::URGENT);
        HybridSimulator simulator(system);
        simulator.configuration().set_step_size(0.25);
        for (unsigned verbosity: {0u, 100u}) {
            LogVerbosityScope logging(verbosity);
            HybridTerminationCriterion stop(Real(1), Integer(3), Set<DiscreteEvent>{jump});
            auto single = simulator.orbit(point(), stop);
            ARIADNE_TESTING_EXPECT_EQUAL(single.size(), 2u);
            ARIADNE_TESTING_EXPECT_EQUAL(single.curves()[1].location(), target);
            expect_endpoint(single.curves()[0], 0.5, 0.5, 2.0, 4.5);
            expect_endpoint(single.curves()[1], 0.5, 0.0, 2.0, 4.0);
            auto multiple = simulator.orbit(PointList{point()}, stop);
            ARIADNE_TESTING_EXPECT_EQUAL(multiple.size(0u), 2u);
            expect_endpoint(multiple.curves(0u)[1], 0.5, 0.0, 2.0, 4.0);
            auto continued = simulator.orbit(point(), HybridTime(1.0_x, 3));
            expect_endpoint(continued.curves()[1], 1.0, 0.0, 2.0, 4.0);
            auto continued_list = simulator.orbit(PointList{point()}, HybridTime(1.0_x, 3));
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
            HybridSimulator simulator(kind == EventKind::INVARIANT
                                          ? static_cast<const HybridAutomatonInterface&>(base)
                                          : static_cast<const HybridAutomatonInterface&>(progress));
            for (unsigned verbosity: {0u, 100u}) {
                LogVerbosityScope logging(verbosity);
                auto stopped = simulator.orbit(point(2.0), HybridTime(1.0_x, 3));
                ARIADNE_TESTING_EXPECT_EQUAL(stopped.curve(0u).size(), 1u);
                expect_endpoint(stopped.curves()[0], 0.0, 2.0, 2.0, 6.0);
                auto stopped_list = simulator.orbit(PointList{point(2.0)}, HybridTime(1.0_x, 3));
                ARIADNE_TESTING_EXPECT_EQUAL(stopped_list.curve(0u, 0u).size(), 1u);
                auto valid = simulator.orbit(point(), HybridTime(0.5_x, 3));
                expect_endpoint(valid.curves()[0], 0.5, 0.5, 2.0, 4.5);
            }
        }
    }
    Void test_box_sampling() {
        auto system = flow_system();
        HybridSimulator simulator(system);
        simulator.configuration().set_fineness(1u);
        HybridUpperBox box(q, RealSpace({
                                  x, y
        }),
                           UpperBoxType({{0, 1}, {1, 2}}));
        SizeType minces = 0u;
        for (auto type: {DiscretizationHybridType::HMince, DiscretizationHybridType::HRecombine}) {
            simulator.configuration().set_d_type(type);
            auto orbit = simulator.orbit(box, q, system.state_space(), HybridTime(0.0_x, 3));
            ARIADNE_TESTING_PROVE(orbit.number_of_trajectories() > 0u);
            if (type == DiscretizationHybridType::HMince)
                minces = orbit.number_of_trajectories();
            else
                ARIADNE_TESTING_PROVE(orbit.number_of_trajectories() <= minces);
            for (SizeType i = 0u; i < orbit.number_of_trajectories(); ++i) {
                auto p = orbit.curve(static_cast<Nat>(i), 0u).begin()->second;
                ARIADNE_TESTING_PROVE(p[0].raw().get_d() >= 0.0 && p[0].raw().get_d() <= 1.0);
                ARIADNE_TESTING_PROVE(p[1].raw().get_d() >= 1.0 && p[1].raw().get_d() <= 2.0);
            }
            std::ostringstream description;
            description << simulator.configuration();
            ARIADNE_TESTING_PROVE(description.str().find(type == DiscretizationHybridType::HMince ? "Mince" : "Recombine") != std::string::npos);
        }
        HybridUpperBox flat(q, RealSpace({
                                   x, y
        }),
                            UpperBoxType({{0, 0}, {1, 2}}));
        auto flat_orbit = simulator.orbit(flat, q, system.state_space(), HybridTime(0.0_x, 3));
        ARIADNE_TESTING_PROVE(flat_orbit.number_of_trajectories() > 0u);
        auto singleton = simulator.orbit(HybridBoundedConstraintSet(q, {x == 0, y == 2}), HybridTime(0.0_x, 3));
        ARIADNE_TESTING_EXPECT_EQUAL(singleton.number_of_trajectories(), 1u);
        expect_endpoint(singleton.curves(0u)[0], 0.0, 0.0, 2.0, 4.0);
    }
};

int main(int argc, char* argv[]) {
    redirect_test_logging_to_null();
    ARIADNE_TESTING_SUITE("test_hybrid_simulator");
    ARIADNE_TESTING_ADD_TEST(TestHybridSimulator, test_labelled_flow);
    ARIADNE_TESTING_ADD_TEST(TestHybridSimulator, test_transitions_and_termination);
    ARIADNE_TESTING_ADD_TEST(TestHybridSimulator, test_invariant_and_progress_stopping);
    ARIADNE_TESTING_ADD_TEST(TestHybridSimulator, test_box_sampling);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
