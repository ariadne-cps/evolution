/***************************************************************************
 *            test_discrete_time_approximate_evolver.cpp
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

/** @file tests/test_discrete_time_approximate_evolver.cpp
 *  @brief Check deterministic sampled tick trajectories, state labels and box sampling.
 */

#include <cmath>
#include <memory>
#include <type_traits>

#include "ariadne-evolution.hpp"
#include "ariadne-testing.hpp"
#include "logging/test_runner.hpp"
#include "system/difference_algebraic_inclusion.hpp"

using namespace Ariadne;

namespace {

using Evolver = DiscreteTimeEvolver<ApproximateTag>;
template<class E>
concept HasSetSemantics = requires(const E& e, const typename E::InitialType& pts,
                                   const typename E::TerminationType& n) {
    e.orbit(pts, n, Semantics::UPPER);
};
static_assert(!HasSetSemantics<Evolver>);
static_assert(std::is_same_v<typename EvolutionTraits<ApproximateTag, DifferenceAlgebraicEquation>::EvolverType,
                             Evolver>);

// Deliberately trajectory- and tick-dependent input choices. These are
// deterministic numerical samples, not UPPER/LOWER interpretations.
class AlternatingInputSampling final : public InputSamplingPolicy {
  public:
    AlternatingInputSampling* clone() const override {
        return new AlternatingInputSampling(*this);
    }
    Vector<FloatDPApproximation> sample(const ExactBoxType&, const InputSamplingContext& context) const override {
        const int sign = (context.step % 2u == 0u) ? 1 : -1;
        const int value = context.trajectory == 0u ? sign : -sign;
        return Vector<FloatDPApproximation>({FloatDPApproximation(value, dp)});
    }
};

class OutsideInputSampling final : public InputSamplingPolicy {
  public:
    OutsideInputSampling* clone() const override { return new OutsideInputSampling(*this); }
    Vector<FloatDPApproximation> sample(const ExactBoxType&, const InputSamplingContext&) const override {
        return Vector<FloatDPApproximation>({FloatDPApproximation(2, dp)});
    }
};

class WrongDimensionInputSampling final : public InputSamplingPolicy {
  public:
    WrongDimensionInputSampling* clone() const override { return new WrongDimensionInputSampling(*this); }
    Vector<FloatDPApproximation> sample(const ExactBoxType&, const InputSamplingContext&) const override {
        return Vector<FloatDPApproximation>(0u, FloatDPApproximation(0, dp));
    }
};

class TestDiscreteTimeApproximateEvolver : public Testing::ContextFixture {
  public:
    using Testing::ContextFixture::ContextFixture;

    Void test_explicit_points_and_polymorphism() {
        const RealVariable x("x"), y("y"), z("z"), missing("missing");
        const DifferenceAlgebraicEquation map(
            {next(x) = x + 1, next(y) = y + x},
            {let(z) = x + 2 * y});
        Evolver evolver(map);
        std::unique_ptr<Evolver::Interface> sampled(evolver.clone());
        ARIADNE_TESTING_EXPECT_EQUAL(sampled->system().state_space(), RealSpace({x, y}));

        const Evolver::ApproximatePointType canonical(
            RealSpace({x, y}), FloatDPApproximationPoint({0.0_x, 2.0_x}, dp));
        const Evolver::ApproximatePointType reordered(
            RealSpace({y, x}), FloatDPApproximationPoint({2.0_x, 0.0_x}, dp));
        const Evolver::InitialType starts{canonical, reordered};

        auto result = sampled->orbit(starts, Integer(3));
        ARIADNE_TESTING_EXPECT_EQUAL(result.curves().size(), 2u);
        for (const auto& curve: result.curves()) {
            ARIADNE_TESTING_EXPECT_EQUAL(curve.size(), 4u);
            ARIADNE_TESTING_EXPECT_WITHIN(curve.end()->first.get_d(), 3.0, 1e-14);
            const auto& state = curve.end()->second;
            ARIADNE_TESTING_EXPECT_WITHIN(state[curve.space()[x]].raw().get_d(), 3.0, 1e-14);
            ARIADNE_TESTING_EXPECT_WITHIN(state[curve.space()[y]].raw().get_d(), 5.0, 1e-14);
            ARIADNE_TESTING_EXPECT_WITHIN(state[curve.space()[z]].raw().get_d(), 13.0, 1e-14);
        }

        auto zero = evolver.orbit(starts, Integer(0));
        ARIADNE_TESTING_EXPECT_EQUAL(zero.curves().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.curve(0u).size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(evolver.orbit(Evolver::InitialType{}, Integer(2)).curves().size(), 0u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.orbit(starts, Integer(-1)));

        const Evolver::ApproximatePointType incompatible(
            RealSpace({missing, y}), FloatDPApproximationPoint({0.0_x, 2.0_x}, dp));
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.orbit(Evolver::InitialType{incompatible}, Integer(1)));
    }

    Void test_paving_independent_sampling() {
        const RealVariable x("x"), y("y"), missing("missing");
        const DifferenceAlgebraicEquation map({next(x) = x, next(y) = y});
        Evolver evolver(map);
        const UpperBoxType box({{0, 1}, {2, 4}});
        auto defaults = evolver.orbit(box, Integer(0));
        ARIADNE_TESTING_EXPECT_EQUAL(defaults.curves().size(), 1u);
        ARIADNE_TESTING_EXPECT_WITHIN(defaults.curve(0u).begin()->second[0u].raw().get_d(), 0.5, 1e-14);
        ARIADNE_TESTING_EXPECT_WITHIN(defaults.curve(0u).begin()->second[1u].raw().get_d(), 3.0, 1e-14);

        evolver.configuration().set_sampling_counts(SamplingCounts{{x, 3u}, {y, 2u}});
        auto samples = evolver.orbit(box, Integer(1));
        ARIADNE_TESTING_EXPECT_EQUAL(samples.curves().size(), 6u);
        for (const auto& curve: samples.curves()) {
            ARIADNE_TESTING_EXPECT_EQUAL(curve.size(), 2u);
            const auto& point = curve.begin()->second;
            ARIADNE_TESTING_PROVE(point[0u].raw().get_d() >= 0.0);
            ARIADNE_TESTING_PROVE(point[0u].raw().get_d() <= 1.0);
            ARIADNE_TESTING_PROVE(point[1u].raw().get_d() >= 2.0);
            ARIADNE_TESTING_PROVE(point[1u].raw().get_d() <= 4.0);
        }

        const UpperBoxType flat({{0, 0}, {2, 4}});
        auto flat_samples = evolver.orbit(flat, Integer(0));
        ARIADNE_TESTING_EXPECT_EQUAL(flat_samples.curves().size(), 2u);

        evolver.configuration().clear_sampling_counts();
        ARIADNE_TESTING_EXPECT_EQUAL(evolver.orbit(flat, Integer(0)).curves().size(), 1u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.configuration().set_sampling_count(x, 0u));
        evolver.configuration().set_sampling_counts(SamplingCounts{{missing, 2u}});
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.orbit(box, Integer(0)));
        evolver.configuration().clear_sampling_counts();
        ARIADNE_TESTING_EXPECT_EXCEPTION(evolver.orbit(UpperBoxType({{0, 1}}), Integer(0)));
    }

    Void test_bounded_input_inclusion_tick_sampling_and_auxiliaries() {
        const RealVariable x("x"), a("a"), u("u");
        const DifferenceAlgebraicInclusion inclusion(
            {next(x) = x + u}, {let(a) = x + u}, RealVariablesBox({-1 <= u <= 1}));
        AlternatingInputSampling policy;
        Evolver evolver(inclusion, policy);
        std::unique_ptr<Evolver::Interface> copy(evolver.clone());
        ARIADNE_TESTING_EXPECT_EQUAL(copy->system().state_space(), RealSpace({x}));
        ARIADNE_TESTING_PROVE(copy->system().has_inclusion());
        const Evolver::InitialType initial{
            Evolver::ApproximatePointType(RealSpace({x}), FloatDPApproximationPoint({0.0_x}, dp)),
            Evolver::ApproximatePointType(RealSpace({x}), FloatDPApproximationPoint({0.0_x}, dp))
        };
        auto orbit = copy->orbit(initial, Integer(2));
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.curves().size(), 2u);
        for (SizeType j = 0u; j < 2u; ++j) {
            const auto& curve = orbit.curve(j);
            ARIADNE_TESTING_EXPECT_EQUAL(curve.size(), 3u);
            const double sign = j == 0u ? 1.0 : -1.0;
            ARIADNE_TESTING_EXPECT_WITHIN(curve.begin()->second[curve.space()[a]].raw().get_d(), sign, 1e-14);
            auto middle = curve.begin();
            ++middle;
            ARIADNE_TESTING_EXPECT_WITHIN(middle->first.get_d(), 1.0, 1e-14);
            ARIADNE_TESTING_EXPECT_WITHIN(middle->second[curve.space()[x]].raw().get_d(), sign, 1e-14);
            ARIADNE_TESTING_EXPECT_WITHIN(middle->second[curve.space()[a]].raw().get_d(), 0.0, 1e-14);
            const auto& last = curve.end()->second;
            ARIADNE_TESTING_EXPECT_WITHIN(last[curve.space()[x]].raw().get_d(), 0.0, 1e-14);
            ARIADNE_TESTING_EXPECT_WITHIN(last[curve.space()[a]].raw().get_d(), sign, 1e-14);
        }

        evolver.configuration().set_sampling_count(x, 3u);
        auto sampled = evolver.orbit(UpperBoxType({{0, 1}}), Integer(1));
        ARIADNE_TESTING_EXPECT_EQUAL(sampled.curves().size(), 3u);

        MidpointInputSampling midpoint;
        Evolver centred(inclusion, midpoint);
        auto centred_orbit = centred.orbit(Evolver::InitialType{initial.front()}, Integer(2));
        const auto& end = centred_orbit.curve(0u).end()->second;
        ARIADNE_TESTING_EXPECT_WITHIN(end[centred_orbit.curve(0u).space()[x]].raw().get_d(), 0.0, 1e-14);
        ARIADNE_TESTING_EXPECT_WITHIN(end[centred_orbit.curve(0u).space()[a]].raw().get_d(), 0.0, 1e-14);

        OutsideInputSampling outside;
        WrongDimensionInputSampling wrong_size;
        Evolver out_of_box(inclusion, outside);
        Evolver invalid_dimension(inclusion, wrong_size);
        ARIADNE_TESTING_EXPECT_EXCEPTION(out_of_box.orbit(initial, Integer(1)));
        ARIADNE_TESTING_EXPECT_EXCEPTION(invalid_dimension.orbit(initial, Integer(1)));
        // No samples are required when no initial trajectories exist.
        ARIADNE_TESTING_EXPECT_EQUAL(evolver.orbit(Evolver::InitialType{}, Integer(2)).curves().size(), 0u);
    }

};

} // namespace

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_discrete_time_approximate_evolver");
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeApproximateEvolver, test_explicit_points_and_polymorphism);
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeApproximateEvolver, test_paving_independent_sampling);
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeApproximateEvolver, test_bounded_input_inclusion_tick_sampling_and_auxiliaries);
    return Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
