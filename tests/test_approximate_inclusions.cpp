/***************************************************************************
 *            test_approximate_inclusions.cpp
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

/** @file tests/test_approximate_inclusions.cpp
 *  @brief Verify sampled ODAI trajectories, policies and hybrid mode changes.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

#include "evolution/input_sampling.hpp"
#include "logging/test_runner.hpp"
#include "threading/thread_manager.hpp"

#include "ariadne-evolution.hpp"
#include "ariadne-testing.hpp"
using namespace Ariadne;

class AlternatingSampling : public InputSamplingPolicy {
  public:
    AlternatingSampling* clone() const override { return new AlternatingSampling(*this); }
    Vector<FloatDPApproximation> sample(const ExactBoxType& domain, const InputSamplingContext& context) const override {
        return Vector<FloatDPApproximation>(domain.dimension(), [&](SizeType i) {
            return FloatDPApproximation(context.step % 2u == 0u ? domain[i].lower_bound() : domain[i].upper_bound(), dp);
        });
    }
};
class UpperSampling : public InputSamplingPolicy {
  public:
    UpperSampling* clone() const override { return new UpperSampling(*this); }
    Vector<FloatDPApproximation> sample(const ExactBoxType& domain, const InputSamplingContext&) const override {
        return Vector<FloatDPApproximation>(domain.dimension(), [&](SizeType i) { return FloatDPApproximation(domain[i].upper_bound(), dp); });
    }
};
class InvalidSampling : public InputSamplingPolicy {
    unsigned _kind;

  public:
    explicit InvalidSampling(unsigned kind) : _kind(kind) {}
    InvalidSampling* clone() const override { return new InvalidSampling(*this); }
    Vector<FloatDPApproximation> sample(const ExactBoxType& domain, const InputSamplingContext&) const override {
        if (_kind == 0u)
            return Vector<FloatDPApproximation>(domain.dimension() + 1u, FloatDPApproximation(0, dp));
        const double value = _kind == 1u ? 100.0 : std::numeric_limits<double>::quiet_NaN();
        return Vector<FloatDPApproximation>(domain.dimension(), FloatDPApproximation(ExactDouble(value), dp));
    }
};
class TestApproximateInclusions : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    using Ordinary = ContinuousTimeEvolver<ApproximateTag>;
    using Hybrid = HybridContinuousTimeEvolver<ApproximateTag>;
    const RealVariable x{"x"}, y{"y"}, u{"u"}, v{"v"}, w{"w"}, p{"p"};
    const TimeVariable t;
    const DiscreteLocation source{1}, middle{2}, target{3};
    void near(double actual, double expected) { ARIADNE_TESTING_PROVE(std::abs(actual - expected) < 1e-12); }
    template<class A, class B>
    void same_points(const A& a, const B& b) {
        ARIADNE_TESTING_REQUIRE(a.size() == b.size());
        for (SizeType i = 0u; i != a.size(); ++i)
            ARIADNE_TESTING_EXPECT_EQUAL(a[i].raw().get_d(), b[i].raw().get_d());
    }
    RealValuation parameters() const { return RealValuation({
        {p, Real(2.5_x)}
 }); }
    OrdinaryDifferentialAlgebraicInclusion inclusion(bool timed = true) const {
        if (timed)
            return OrdinaryDifferentialAlgebraicInclusion({dot(x) = p + t + u, dot(y) = sqr(u)}, RealVariablesBox({-1 <= u <= 1}), RealVariablesBox({2 <= p <= 3}), t);
        return OrdinaryDifferentialAlgebraicInclusion({dot(x) = p + u, dot(y) = sqr(u)}, RealVariablesBox({-1 <= u <= 1}), RealVariablesBox({2 <= p <= 3}));
    }
    Ordinary::ApproximatePointType point(double initial = 0.0) const {
        return Ordinary::ApproximatePointType(RealSpace({y, x}), FloatDPApproximationPoint({0.0_x, ExactDouble(initial)}, dp));
    }
    Hybrid::HybridApproximatePointType hybrid_point() const {
        return Hybrid::HybridApproximatePointType(source, RealSpace({y, x}), FloatDPApproximationPoint({0.0_x, 0.0_x}, dp));
    }
    void test_policy_boxes() {
        const ExactBoxType domain({
            {-2, 3},
            { 7, 8},
            { 1, 1}
        });
        const InputSamplingContext context{4u, 2u, FloatDPApproximation(0.5_x, dp), FloatDPApproximation(0.125_x, dp)};
        UniformInputSampling uniform(42u);
        const auto values = uniform.sample(domain, context), repeated = uniform.sample(domain, context);
        same_points(values, repeated);
        ARIADNE_TESTING_EXPECT_EQUAL(values.size(), 3u);
        for (SizeType i = 0u; i != 3u; ++i) {
            const double value = values[i].raw().get_d();
            ARIADNE_TESTING_PROVE(value >= domain[i].lower_bound().get_d() && value <= domain[i].upper_bound().get_d());
        }
        struct RestoreRounding {
            FloatDP::RoundingModeType mode = FloatDP::get_rounding_mode();
            ~RestoreRounding() { FloatDP::set_rounding_mode(mode); }
        } restore_rounding;
        for (SizeType step = 0u; step != 32u; ++step) {
            auto request = context;
            request.step = step;
            FloatDP::set_rounding_mode(FloatDP::ROUND_UPWARD);
            const auto sampled_up = uniform.sample(domain, request);
            const auto upward_mode = FloatDP::get_rounding_mode();
            FloatDP::set_rounding_mode(FloatDP::ROUND_DOWNWARD);
            const auto sampled_down = uniform.sample(domain, request);
            const auto downward_mode = FloatDP::get_rounding_mode();
            same_points(sampled_up, sampled_down);
            ARIADNE_TESTING_EXPECT_EQUAL(upward_mode, static_cast<int>(FloatDP::ROUND_UPWARD));
            ARIADNE_TESTING_EXPECT_EQUAL(downward_mode, static_cast<int>(FloatDP::ROUND_DOWNWARD));
        }
        const auto center = MidpointInputSampling().sample(domain, context);
        near(center[0u].raw().get_d(), 0.5);
        near(center[1u].raw().get_d(), 7.5);
        near(center[2u].raw().get_d(), 1);
    }
    void test_ordinary_piecewise_inputs_and_external_time() {
        for (bool timed: {false, true}) {
            Ordinary evolver(inclusion(timed), AlternatingSampling(), parameters());
            evolver.configuration().set_step_size(0.0625);
            std::unique_ptr<Ordinary::Interface> clone(evolver.clone());
            ARIADNE_TESTING_EXPECT_EQUAL(clone->system().dimension(), 2u);
            ARIADNE_TESTING_EXPECT_EQUAL(clone->system().auxiliary_space().dimension(), 0u);
            auto orbit = clone->orbit(Ordinary::ApproximateListPointType{point(), point(1)}, Real(0.25_x));
            for (SizeType i = 0u; i != 2u; ++i) {
                const auto& curve = orbit.curve(i);
                near(curve.end()->first.get_d(), 0.25);
                near(curve.end()->second[curve.space()[x]].raw().get_d(), static_cast<double>(i) + (timed ? 0.65625 : 0.625));
                near(curve.end()->second[curve.space()[y]].raw().get_d(), 0.25);
                ARIADNE_TESTING_EXPECT_EQUAL(curve.space().dimension(), 2u);
            }
            auto held = evolver.orbit(Ordinary::ApproximateListPointType{point()}, Real(0));
            ARIADNE_TESTING_EXPECT_EQUAL(held.curve(0u).size(), 1u);
        }
    }
    void test_provider_coordinate_normalization_and_midpoint() {
        const OrdinaryDifferentialAlgebraicInclusion scaled({dot(x) = 3 * u}, RealVariablesBox({2 <= u <= 4}));
        Ordinary upper(scaled, UpperSampling()), center(scaled, MidpointInputSampling());
        const Ordinary::ApproximateListPointType initial{Ordinary::ApproximatePointType(RealSpace({x}), FloatDPApproximationPoint({0.0_x}, dp))};
        auto a = upper.orbit(initial, Real(0.25_x)), b = center.orbit(initial, Real(0.25_x));
        near(a.curve(0u).end()->second[0u].raw().get_d(), 3.0);
        near(b.curve(0u).end()->second[0u].raw().get_d(), 2.25);
    }
    void test_uniform_reproducibility_across_scheduling() {
        struct RestoreConcurrency {
            std::size_t value = ThreadManager::instance().concurrency();
            ~RestoreConcurrency() { ThreadManager::instance().set_concurrency(value); }
        } restore;
        Ordinary evolver(inclusion(), UniformInputSampling(42u), parameters());
        evolver.configuration().set_step_size(0.0625);
        const Ordinary::ApproximateListPointType initial{point(), point(1), point(2)};
        ThreadManager::instance().set_concurrency(1u);
        auto sequential = evolver.orbit(initial, Real(0.25_x));
        ThreadManager::instance().set_concurrency(std::min<std::size_t>(2u, ThreadManager::instance().maximum_concurrency()));
        auto parallel = evolver.orbit(initial, Real(0.25_x));
        std::unique_ptr<Ordinary> clone(evolver.clone());
        auto cloned = clone->orbit(initial, Real(0.25_x));
        Ordinary other(inclusion(), UniformInputSampling(43u), parameters());
        other.configuration().set_step_size(0.0625);
        auto different = other.orbit(initial, Real(0.25_x));
        bool changed = false;
        for (SizeType i = 0u; i != initial.size(); ++i) {
            const auto& reference = sequential.curve(i);
            ARIADNE_TESTING_REQUIRE(reference.size() == parallel.curve(i).size());
            same_points(reference.end()->second, parallel.curve(i).end()->second);
            auto expected = reference.begin(), actual = parallel.curve(i).begin();
            for (; expected != reference.end(); ++expected, ++actual) {
                ARIADNE_TESTING_EXPECT_EQUAL(expected->first, actual->first);
                same_points(expected->second, actual->second);
            }
            same_points(reference.end()->second, cloned.curve(i).end()->second);
            const double value = reference.end()->second[reference.space()[x]].raw().get_d();
            ARIADNE_TESTING_PROVE(value >= static_cast<double>(i) + 0.40625 && value <= static_cast<double>(i) + 0.90625);
            changed = changed || value != different.curve(i).end()->second[different.curve(i).space()[x]].raw().get_d();
        }
        ARIADNE_TESTING_PROVE(changed);
    }
    void test_hybrid_flow_single_and_multiple() {
        for (bool timed: {false, true}) {
            HybridAutomaton system;
            system.new_mode(source, inclusion(timed));
            Hybrid evolver(system, AlternatingSampling(), parameters());
            evolver.configuration().set_step_size(0.0625);
            const auto single = evolver.orbit(hybrid_point(), HybridTime(0.25_x, 3));
            auto multiple = evolver.orbit(Hybrid::HybridApproximateListPointType{hybrid_point(), hybrid_point()}, HybridTime(0.25_x, 3));
            const auto& curve = single.curves()[0u];
            near(curve.euclidean_set().end()->second[curve.space()[x]].raw().get_d(), timed ? 0.65625 : 0.625);
            near(curve.euclidean_set().end()->second[curve.space()[y]].raw().get_d(), 0.25);
            for (Nat i = 0u; i != 2u; ++i)
                same_points(curve.euclidean_set().end()->second, multiple.curves(i)[0u].euclidean_set().end()->second);
            Hybrid uniform(system, UniformInputSampling(42u), parameters());
            uniform.configuration().set_step_size(0.0625);
            const auto one = uniform.orbit(hybrid_point(), HybridTime(0.25_x, 3));
            auto batch = uniform.orbit(Hybrid::HybridApproximateListPointType{hybrid_point(), hybrid_point()}, HybridTime(0.25_x, 3));
            auto repeated = uniform.orbit(Hybrid::HybridApproximateListPointType{hybrid_point(), hybrid_point()}, HybridTime(0.25_x, 3));
            same_points(one.curves()[0u].euclidean_set().end()->second, batch.curves(0u)[0u].euclidean_set().end()->second);
            for (Nat i = 0u; i != 2u; ++i)
                same_points(batch.curves(i)[0u].euclidean_set().end()->second, repeated.curves(i)[0u].euclidean_set().end()->second);
            ARIADNE_TESTING_PROVE(batch.curves(0u)[0u].euclidean_set().end()->second[0u].raw().get_d() != batch.curves(1u)[0u].euclidean_set().end()->second[0u].raw().get_d());
        }
    }
    HybridAutomaton transitions() const {
        HybridAutomaton system;
        system.new_mode(source, OrdinaryDifferentialAlgebraicEquation({dot(x) = 1}, RealVariablesBox({2 <= p <= 3})));
        system.new_mode(middle, OrdinaryDifferentialAlgebraicInclusion({dot(x) = u}, RealVariablesBox({-1 <= u <= 1}), RealVariablesBox({2 <= p <= 3})));
        system.new_mode(target, OrdinaryDifferentialAlgebraicInclusion({dot(x) = v + sqr(w) + t}, RealVariablesBox({-1 <= v <= 1, -2 <= w <= 2})));
        system.new_transition(source, DiscreteEvent("enter"), middle, {next(x) = p + t}, t >= Real(0.0625_x), EventKind::URGENT);
        system.new_transition(middle, DiscreteEvent("change"), target, {next(x) = x + p + t}, t >= Real(0.09375_x), EventKind::URGENT);
        system.new_invariant(middle, p <= 3, DiscreteEvent("parameter_bound"));
        return system;
    }
    void check_transitions(const List<HybridInterpolatedCurve>& curves) {
        ARIADNE_TESTING_REQUIRE(curves.size() == 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(curves[1u].location(), middle);
        near(curves[1u].euclidean_set().begin()->first.get_d(), 0.09375);
        near(curves[1u].euclidean_set().end()->first.get_d(), 0.125);
        near(curves[1u].euclidean_set().end()->second[0u].raw().get_d(), 2.625);
        ARIADNE_TESTING_EXPECT_EQUAL(curves[2u].location(), target);
        near(curves[2u].euclidean_set().begin()->second[0u].raw().get_d(), 5.25);
        near(curves[2u].euclidean_set().end()->first.get_d(), 0.25);
        near(curves[2u].euclidean_set().end()->second[0u].raw().get_d(), 5.7734375);
        ARIADNE_TESTING_EXPECT_EQUAL(curves[2u].space().dimension(), 1u);
    }
    void test_hybrid_transitions_clocks_and_input_dimensions() {
        const auto system = transitions();
        Hybrid evolver(system, AlternatingSampling(), parameters());
        evolver.configuration().set_step_size(0.03125);
        const auto start = Hybrid::HybridApproximatePointType(source, RealSpace({x}), FloatDPApproximationPoint({0.0_x}, dp));
        const auto single = evolver.orbit(start, HybridTime(0.25_x, 3));
        check_transitions(single.curves());
        auto multiple = evolver.orbit(Hybrid::HybridApproximateListPointType{start, start}, HybridTime(0.25_x, 3));
        check_transitions(multiple.curves(0u));
        check_transitions(multiple.curves(1u));
        Hybrid no_policy(system, parameters());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, no_policy.orbit(start, HybridTime(0.25_x, 3)));
        const auto held = no_policy.orbit(start, HybridTime(0u, 3));
        ARIADNE_TESTING_EXPECT_EQUAL(held.curves().size(), 1u);
    }
    void test_sampler_output_and_representation_gates() {
        const Ordinary::ApproximateListPointType initial{point()};
        HybridAutomaton system;
        system.new_mode(source, inclusion());
        for (unsigned kind: {0u, 1u, 2u}) {
            Ordinary evolver(inclusion(), InvalidSampling(kind), parameters());
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, evolver.orbit(initial, Real(0.125_x)));
            auto held = evolver.orbit(initial, Real(0));
            ARIADNE_TESTING_EXPECT_EQUAL(held.curve(0u).size(), 1u);
            Hybrid hybrid(system, InvalidSampling(kind), parameters());
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, hybrid.orbit(hybrid_point(), HybridTime(0.125_x, 3)));
        }
        Hybrid missing(system, parameters());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, missing.orbit(hybrid_point(), HybridTime(0u, 3)));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, Ordinary(inclusion(), MidpointInputSampling()));
        const ExactBoxType domain({
            {-1, 1}
        });
        const auto uncertain = ValidatedIntervalTaylorFunctionModel<FloatDP>::constant(domain, FloatDPUpperInterval(-1, 1, dp), ThresholdSweeper<FloatDP>(dp, 1e-10));
        const OrdinaryDifferentialAlgebraicInclusion interval(RealSpace({x}), ValidatedVectorIntervalTaylorFunctionModel<FloatDP>({uncertain}));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, Ordinary(interval, MidpointInputSampling()));
        HybridAutomaton interval_system;
        interval_system.new_mode(source, interval);
        Hybrid interval_hybrid(interval_system, MidpointInputSampling());
        const auto interval_start = Hybrid::HybridApproximatePointType(source, RealSpace({x}), FloatDPApproximationPoint({0.0_x}, dp));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, interval_hybrid.orbit(interval_start, HybridTime(0u, 3)));
    }
};
int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_approximate_inclusions");
    ARIADNE_TESTING_ADD_TEST(TestApproximateInclusions, test_policy_boxes);
    ARIADNE_TESTING_ADD_TEST(TestApproximateInclusions, test_ordinary_piecewise_inputs_and_external_time);
    ARIADNE_TESTING_ADD_TEST(TestApproximateInclusions, test_provider_coordinate_normalization_and_midpoint);
    ARIADNE_TESTING_ADD_TEST(TestApproximateInclusions, test_uniform_reproducibility_across_scheduling);
    ARIADNE_TESTING_ADD_TEST(TestApproximateInclusions, test_hybrid_flow_single_and_multiple);
    ARIADNE_TESTING_ADD_TEST(TestApproximateInclusions, test_hybrid_transitions_clocks_and_input_dimensions);
    ARIADNE_TESTING_ADD_TEST(TestApproximateInclusions, test_sampler_output_and_representation_gates);
    return Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
