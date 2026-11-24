/***************************************************************************
 *            tests/test_simulator_sampling_contract.cpp
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

/** @file tests/test_simulator_sampling_contract.cpp
 *  @brief Verify sampling cardinality, ordering, and domain membership.
 */

#include "simulator_sampling.hpp"

#include "logging/test_runner.hpp"
#include "symbolic/space.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {
class TestSimulatorSamplingContract : public Testing::ContextFixture {
  public:
    using Testing::ContextFixture::ContextFixture;

    Void test_continuous_counts_and_membership() {
        const RealVariable x("x"), y("y"), missing("missing");
        const RealSpace space({x, y});
        const UpperBoxType domain({{0, 1}, {-2, 2}});

        const auto defaults = EvolutionDetail::sample_vector_field_initial_points(
            domain, space, SamplingCounts{});
        ARIADNE_TESTING_EXPECT_EQUAL(defaults.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(defaults.at(0)[0].raw().get_d(), 0.5);

        const auto samples = EvolutionDetail::sample_vector_field_initial_points(
            domain, space, SamplingCounts{{x, 3u}, {y, 2u}});
        ARIADNE_TESTING_EXPECT_EQUAL(samples.size(), 6u);
        for (SizeType i = 0u; i < samples.size(); ++i) {
            const auto& point = samples.at(i);
            ARIADNE_TESTING_PROVE(point.space() == space);
            const double xv = point[0].raw().get_d();
            const double yv = point[1].raw().get_d();
            ARIADNE_TESTING_PROVE(xv >= 0.0 && xv <= 1.0);
            ARIADNE_TESTING_PROVE(yv >= -2.0 && yv <= 2.0);
            if (i != 0u)
                ARIADNE_TESTING_PROVE(xv >= samples.at(i - 1u)[0].raw().get_d());
        }
        ARIADNE_TESTING_EXPECT_EXCEPTION(
            EvolutionDetail::sample_vector_field_initial_points(
                domain, space, SamplingCounts{{missing, 2u}}));
        ARIADNE_TESTING_EXPECT_EXCEPTION(
            EvolutionDetail::sample_vector_field_initial_points(
                domain, space, SamplingCounts{{x, 0u}}));
    }

    Void test_hybrid_counts_and_flat_coordinate() {
        const RealVariable x("x"), y("y");
        const RealSpace space({x, y});
        const DiscreteLocation loc(1);
        const UpperBoxType domain({{0, 1}, {1, 2}});
        const auto samples = EvolutionDetail::sample_hybrid_initial_points(
            domain, loc, space, SamplingCounts{{x, 3u}, {y, 2u}});
        ARIADNE_TESTING_EXPECT_EQUAL(samples.size(), 6u);
        for (SizeType i = 0u; i < samples.size(); ++i) {
            const auto& point = samples[i];
            ARIADNE_TESTING_PROVE(point.location() == loc);
            ARIADNE_TESTING_PROVE(point.space() == space);
            const double xv = point.point()[0].raw().get_d();
            const double yv = point.point()[1].raw().get_d();
            ARIADNE_TESTING_PROVE(xv >= 0.0 && xv <= 1.0);
            ARIADNE_TESTING_PROVE(yv >= 1.0 && yv <= 2.0);
        }
        const UpperBoxType flat({{0, 0}, {1, 2}});
        const auto flat_points = EvolutionDetail::sample_hybrid_initial_points(
            flat, loc, space, SamplingCounts{{x, 5u}, {y, 2u}});
        ARIADNE_TESTING_EXPECT_EQUAL(flat_points.size(), 2u);
        for (SizeType i = 0u; i < flat_points.size(); ++i)
            ARIADNE_TESTING_EXPECT_EQUAL(flat_points[i].point()[0].raw().get_d(), 0.0);
        const UpperBoxType singleton({{0, 0}, {2, 2}});
        const auto one = EvolutionDetail::sample_hybrid_initial_points(
            singleton, loc, space, SamplingCounts{{x, 5u}, {y, 2u}});
        ARIADNE_TESTING_EXPECT_EQUAL(one.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(one[0].point()[1].raw().get_d(), 2.0);
    }
};
} // namespace

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_simulator_sampling_contract");
    ARIADNE_TESTING_ADD_TEST(TestSimulatorSamplingContract, test_continuous_counts_and_membership);
    ARIADNE_TESTING_ADD_TEST(TestSimulatorSamplingContract, test_hybrid_counts_and_flat_coordinate);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
