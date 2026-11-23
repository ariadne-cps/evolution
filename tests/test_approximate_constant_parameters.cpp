/***************************************************************************
 *            test_approximate_constant_parameters.cpp
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

#include <cmath>
#include <memory>
#include <stdexcept>

#include "evolution/continuous_time_approximate_evolver.hpp"
#include "evolution/hybrid_continuous_time_approximate_evolver.hpp"
#include "evolution/hybrid_orbit.hpp"
#include "logging/test_runner.hpp"
#include "system/hybrid_automata.hpp"
#include "system/odai.hpp"

#include "ariadne-testing.hpp"
using namespace Ariadne;

class TestApproximateConstantParameters : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    using Ordinary = ContinuousTimeEvolver<ApproximateTag>;
    using Hybrid = HybridContinuousTimeEvolver<ApproximateTag>;
    const RealVariable x{"x"}, y{"y"}, a{"a"}, p{"p"}, q{"q"}, unused{"unused"};
    const TimeVariable t;
    const DiscreteLocation source{1}, middle{2}, target{3};
    RealValuation selection() const { return RealValuation({
        {q, Real(-0.5_x)},
        {p,  Real(2.5_x)}
 }); }
    OrdinaryDifferentialAlgebraicEquation flow(bool timed) const {
        const RealVariablesBox parameters({2 <= p <= 3, -1 <= q <= 0});
        if (timed)
            return OrdinaryDifferentialAlgebraicEquation({dot(x) = p + q * t, dot(y) = 0}, parameters, t, {let(a) = x + 2 * p + q + t});
        return OrdinaryDifferentialAlgebraicEquation({dot(x) = p + q, dot(y) = 0}, parameters, std::nullopt, {let(a) = x + 2 * p + q});
    }
    Ordinary::ApproximatePointType ordinary_point(double xv = 1.0) const {
        return Ordinary::ApproximatePointType(RealSpace({y, x}), FloatDPApproximationPoint({3.0_x, ExactDouble(xv)}, dp));
    }
    Hybrid::HybridApproximatePointType hybrid_point(double xv = 1.0) const {
        return Hybrid::HybridApproximatePointType(source, RealSpace({y, x}), FloatDPApproximationPoint({3.0_x, ExactDouble(xv)}, dp));
    }
    void near(double actual, double expected) {
        ARIADNE_TESTING_PROVE(std::abs(actual - expected) < 1e-12);
    }
    const InterpolatedCurve& points(const LabelledInterpolatedCurve& curve) const { return curve; }
    const InterpolatedCurve& points(const HybridInterpolatedCurve& curve) const { return curve.euclidean_set(); }
    template<class Curve>
    void flow_endpoint(const Curve& curve, bool timed, double initial) {
        const auto& end = points(curve).end();
        near(end->first.get_d(), 0.25);
        const double expected = initial + (timed ? 0.609375 : 0.5);
        near(end->second[curve.space()[x]].raw().get_d(), expected);
        near(end->second[curve.space()[y]].raw().get_d(), 3);
        near(end->second[curve.space()[a]].raw().get_d(), expected + (timed ? 4.75 : 4.5));
        ARIADNE_TESTING_EXPECT_EQUAL(curve.space().dimension(), 3u);
    }
    void test_ordinary_values_are_fixed_and_named() {
        for (bool timed: {false, true}) {
            auto values = selection();
            Ordinary evolver(flow(timed), values);
            evolver.configuration().set_step_size(0.0625);
            values.set(p, Real(3)); // The caller cannot mutate the selected choices.
            near(FloatDPApproximation(evolver.constant_parameters().get(p), dp).raw().get_d(), 2.5);
            std::unique_ptr<Ordinary> clone(evolver.clone());
            const Ordinary::ApproximateListPointType initial{ordinary_point(), ordinary_point(-1)};
            auto orbit = clone->orbit(initial, Real(0.25_x));
            ARIADNE_TESTING_EXPECT_EQUAL(orbit.curves().size(), 2u);
            flow_endpoint(orbit.curve(0u), timed, 1);
            flow_endpoint(orbit.curve(1u), timed, -1);
            auto repeated = evolver.orbit(initial, Real(0.25_x));
            flow_endpoint(repeated.curve(0u), timed, 1);
            auto zero = evolver.orbit(Ordinary::ApproximateListPointType{ordinary_point()}, Real(0));
            ARIADNE_TESTING_EXPECT_EQUAL(zero.curve(0u).size(), 1u);
            near(zero.curve(0u).begin()->second[zero.curve(0u).space()[a]].raw().get_d(), 5.5);
        }
    }
    void test_ordinary_selection_validation() {
        const auto equation = flow(false);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, (void)Ordinary(equation));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, Ordinary(equation, RealValuation({
                                                                                    {p, Real(2.5_x)}
        })));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, Ordinary(equation, RealValuation({
                                                                                    {     p, Real(2.5_x)},
                                                                                    {unused,     Real(0)}
        })));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, Ordinary(equation, RealValuation({
                                                                                    {p, Real(4)},
                                                                                    {q, Real(0)}
        })));
        auto extra = selection();
        extra.set(unused, Real(1));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, Ordinary(equation, extra));
        const OrdinaryDifferentialAlgebraicEquation simple({dot(x) = p}, RealVariablesBox({2 <= p <= 3}));
        for (int value: {2, 3}) {
            Ordinary evolver(simple, RealValuation({
                                         {p, Real(value)}
            }));
            auto orbit = evolver.orbit(Ordinary::ApproximateListPointType{Ordinary::ApproximatePointType(RealSpace({x}), FloatDPApproximationPoint({0.0_x}, dp))}, Real(0.25_x));
            near(orbit.curve(0u).end()->second[0u].raw().get_d(), value * 0.25);
        }
        const OrdinaryDifferentialAlgebraicEquation no_parameters({dot(x) = 0});
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, Ordinary(no_parameters, selection()));
    }
    void test_hybrid_autonomous_and_timed_flow() {
        for (bool timed: {false, true}) {
            HybridAutomaton system;
            system.new_mode(source, flow(timed));
            Hybrid evolver(system, selection());
            evolver.configuration().set_step_size(0.0625);
            std::unique_ptr<Hybrid> clone(evolver.clone());
            const auto single = clone->orbit(hybrid_point(), HybridTime(0.25_x, 3));
            flow_endpoint(single.curves()[0u], timed, 1);
            auto multiple = evolver.orbit(Hybrid::HybridApproximateListPointType{hybrid_point(), hybrid_point(-1)}, HybridTime(0.25_x, 3));
            ARIADNE_TESTING_EXPECT_EQUAL(multiple.number_of_trajectories(), 2u);
            flow_endpoint(multiple.curves(0u)[0u], timed, 1);
            flow_endpoint(multiple.curves(1u)[0u], timed, -1);
            const auto zero = evolver.orbit(hybrid_point(), HybridTime(0u, 3));
            near(zero.curves()[0u].euclidean_set().begin()->second[zero.curves()[0u].space()[a]].raw().get_d(), 5.5);
            near(FloatDPApproximation(evolver.constant_parameters().get(p), dp).raw().get_d(), 2.5);
            Hybrid missing(system);
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, missing.orbit(hybrid_point(), HybridTime(0u, 3)));
        }
    }
    HybridAutomaton transition_system() const {
        HybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = p}, RealVariablesBox({2 <= p <= 3}), std::nullopt, {let(a) = x + p + t}));
        system.new_mode(middle, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, {let(a) = x + t}));
        system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = q + t}, RealVariablesBox({-1 <= q <= 0}), t, {let(a) = x + q + t}));
        system.new_invariant(source, x <= p, DiscreteEvent("bound"));
        system.new_invariant(target, q <= 0, DiscreteEvent("target_bound"));
        system.new_transition(source, DiscreteEvent("out"), middle, {next(x) = x + p + t}, t >= p / 32, EventKind::URGENT);
        system.new_transition(middle, DiscreteEvent("back"), target, {next(x) = x + 1}, x >= 2, EventKind::URGENT);
        return system;
    }
    void check_transition(const List<HybridInterpolatedCurve>& curves) {
        ARIADNE_TESTING_REQUIRE(curves.size() == 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(curves[1u].location(), middle);
        ARIADNE_TESTING_EXPECT_EQUAL(curves[2u].location(), target);
        const auto& entry = curves[2u].euclidean_set().begin();
        near(entry->first.get_d(), 0.09375);
        near(entry->second[curves[2u].space()[x]].raw().get_d(), 3.828125);
        const auto& end = curves[2u].euclidean_set().end();
        near(end->first.get_d(), 0.25);
        near(end->second[curves[2u].space()[x]].raw().get_d(), 3.77685546875);
        near(end->second[curves[2u].space()[a]].raw().get_d(), 3.52685546875);
        ARIADNE_TESTING_EXPECT_EQUAL(curves[2u].space().dimension(), 2u);
    }
    void test_hybrid_guards_resets_and_partial_mode_tuples() {
        Hybrid evolver(transition_system(), selection());
        evolver.configuration().set_step_size(0.03125);
        const auto start = Hybrid::HybridApproximatePointType(source, RealSpace({x}), FloatDPApproximationPoint({0.0_x}, dp));
        const auto single = evolver.orbit(start, HybridTime(0.25_x, 3));
        check_transition(single.curves());
        auto multiple = evolver.orbit(Hybrid::HybridApproximateListPointType{start, start}, HybridTime(0.25_x, 3));
        check_transition(multiple.curves(0u));
        check_transition(multiple.curves(1u));
        Hybrid missing_later(transition_system(), RealValuation({
                                                      {p, Real(2.5_x)}
        }));
        missing_later.configuration().set_step_size(0.03125);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, missing_later.orbit(start, HybridTime(0.25_x, 3)));
    }
    void test_hybrid_parameter_only_guard_and_shared_value() {
        HybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, RealVariablesBox({2 <= p <= 3})));
        system.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}, RealVariablesBox({1 <= p <= 4}), std::nullopt, {let(a) = x - p}));
        system.new_transition(source, DiscreteEvent("selected"), target, {next(x) = p}, p >= Real(2.25_x), EventKind::URGENT);
        Hybrid evolver(system, RealValuation({
                                   {p, Real(2.5_x)}
        }));
        const auto start = Hybrid::HybridApproximatePointType(source, RealSpace({x}), FloatDPApproximationPoint({0.0_x}, dp));
        const auto single = evolver.orbit(start, HybridTime(0.25_x, 3));
        ARIADNE_TESTING_REQUIRE(single.curves().size() == 2u);
        const auto& final = single.curves()[1u];
        near(final.euclidean_set().end()->second[final.space()[x]].raw().get_d(), 2.5);
        near(final.euclidean_set().end()->second[final.space()[a]].raw().get_d(), 0);
        auto multiple = evolver.orbit(Hybrid::HybridApproximateListPointType{start, start}, HybridTime(0.25_x, 3));
        for (Nat i = 0u; i != 2u; ++i) {
            const auto& end = multiple.curves(i)[1u];
            near(end.euclidean_set().end()->second[end.space()[x]].raw().get_d(), 2.5);
            near(end.euclidean_set().end()->second[end.space()[a]].raw().get_d(), 0);
        }
    }
    void test_hybrid_invariants_domains_and_first_use() {
        HybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = p}, RealVariablesBox({2 <= p <= 3}), std::nullopt, {let(a) = p}));
        system.new_invariant(source, p <= Real(2.25_x), DiscreteEvent("bound"));
        const auto start = Hybrid::HybridApproximatePointType(source, RealSpace({x}), FloatDPApproximationPoint({0.0_x}, dp));
        Hybrid stopped(system, selection());
        const auto result = stopped.orbit(start, HybridTime(0.25_x, 3));
        near(result.curves()[0u].euclidean_set().end()->first.get_d(), 0);
        Hybrid outside(system, RealValuation({
                                   {p, Real(4)}
        }));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, outside.orbit(start, HybridTime(0u, 3)));

        HybridAutomaton introduced;
        introduced.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = 0}));
        introduced.new_mode(target, OrdinaryDifferentialAlgebraicEquation({dot(x) = p}, RealVariablesBox({2 <= p <= 3}), std::nullopt, {let(a) = p}));
        introduced.new_transition(source, DiscreteEvent("bind"), target, {next(x) = x}, x >= -1, EventKind::URGENT);
        Hybrid chosen(introduced, RealValuation({
                                      {p, Real(2.5_x)}
        }));
        const auto introduced_orbit = chosen.orbit(start, HybridTime(0.25_x, 3));
        ARIADNE_TESTING_REQUIRE(introduced_orbit.curves().size() == 2u);
        near(introduced_orbit.curves()[1u].euclidean_set().end()->second[introduced_orbit.curves()[1u].space()[x]].raw().get_d(), 0.625);
        Hybrid not_chosen(introduced);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, not_chosen.orbit(start, HybridTime(0.25_x, 3)));
        const auto held = not_chosen.orbit(start, HybridTime(0u, 3));
        ARIADNE_TESTING_EXPECT_EQUAL(held.curves().size(), 1u);
        Hybrid outside_target(introduced, RealValuation({
                                              {p, Real(4)}
        }));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, outside_target.orbit(start, HybridTime(0.25_x, 3)));
        HybridAutomaton inclusion;
        const RealVariable u("u");
        inclusion.new_mode(source, OrdinaryDifferentialAlgebraicInclusion({dot(x) = p + u}, RealVariablesBox({-1 <= u <= 1}), RealVariablesBox({2 <= p <= 3})));
        Hybrid unsupported(inclusion, RealValuation({
                                          {p, Real(2.5_x)}
        }));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, unsupported.orbit(start, HybridTime(0u, 3)));
    }
};
int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_approximate_constant_parameters");
    ARIADNE_TESTING_ADD_TEST(TestApproximateConstantParameters, test_ordinary_values_are_fixed_and_named);
    ARIADNE_TESTING_ADD_TEST(TestApproximateConstantParameters, test_ordinary_selection_validation);
    ARIADNE_TESTING_ADD_TEST(TestApproximateConstantParameters, test_hybrid_autonomous_and_timed_flow);
    ARIADNE_TESTING_ADD_TEST(TestApproximateConstantParameters, test_hybrid_guards_resets_and_partial_mode_tuples);
    ARIADNE_TESTING_ADD_TEST(TestApproximateConstantParameters, test_hybrid_parameter_only_guard_and_shared_value);
    ARIADNE_TESTING_ADD_TEST(TestApproximateConstantParameters, test_hybrid_invariants_domains_and_first_use);
    return Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
