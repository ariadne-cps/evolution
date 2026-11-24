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

    Void test_continuous_resolution_and_subspace() {
        const RealVariable x("x"), y("y"), missing("missing");
        const RealSpace space({x, y});
        const UpperBoxType domain({{0, 1}, {-2, 2}});
        for (auto type: {DiscretisationType::None,
                         DiscretisationType::Mince,
                         DiscretisationType::Recombine}) {
            const auto samples = EvolutionDetail::sample_vector_field_initial_points(
                domain, space, subspace{}, type, 2u);
            const SizeType expected = type == DiscretisationType::Mince ? 16u :
                type == DiscretisationType::Recombine ? 4u : 1u;
            ARIADNE_TESTING_EXPECT_EQUAL(samples.size(), expected);
            for (const auto& sample: samples) {
                ARIADNE_TESTING_PROVE(sample.space() == space);
                const double xv = sample[0].raw().get_d();
                const double yv = sample[1].raw().get_d();
                ARIADNE_TESTING_PROVE(xv >= 0.0 && xv <= 1.0);
                ARIADNE_TESTING_PROVE(yv >= -2.0 && yv <= 2.0);
            }
        }

        const auto selected = EvolutionDetail::sample_vector_field_initial_points(
            domain, space, subspace{{x, 3u}, {y, 2u}},
            DiscretisationType::Mince, 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(selected.size(), 6u);
        for (SizeType i = 1u; i < selected.size(); ++i) {
            const double previous = selected.at(i - 1u)[0].raw().get_d();
            const double current = selected.at(i)[0].raw().get_d();
            ARIADNE_TESTING_PROVE(current >= previous);
        }
        ARIADNE_TESTING_EXPECT_EXCEPTION(
            EvolutionDetail::sample_vector_field_initial_points(
                domain, space, subspace{{missing, 2u}},
                DiscretisationType::Mince, 0u));
    }

    Void test_hybrid_resolution_and_flat_coordinate() {
        const RealVariable x("x"), y("y");
        const RealSpace space({x, y});
        const DiscreteLocation loc(1);
        const UpperBoxType domain({{0, 1}, {1, 2}});
        const auto fine = EvolutionDetail::sample_hybrid_initial_points(
            domain, loc, space, DiscretizationHybridType::HMince, 2u);
        const auto coarse = EvolutionDetail::sample_hybrid_initial_points(
            domain, loc, space, DiscretizationHybridType::HRecombine, 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(fine.size(), 16u);
        ARIADNE_TESTING_EXPECT_EQUAL(coarse.size(), 4u);
        for (const auto& point: fine) {
            ARIADNE_TESTING_PROVE(point.location() == loc);
            ARIADNE_TESTING_PROVE(point.space() == space);
            const double xv = point.point()[0].raw().get_d();
            const double yv = point.point()[1].raw().get_d();
            ARIADNE_TESTING_PROVE(xv >= 0.0 && xv <= 1.0);
            ARIADNE_TESTING_PROVE(yv >= 1.0 && yv <= 2.0);
        }
        const UpperBoxType flat({{0, 0}, {1, 2}});
        const auto flat_points = EvolutionDetail::sample_hybrid_initial_points(
            flat, loc, space, DiscretizationHybridType::HMince, 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(flat_points.size(), 4u);
        for (const auto& point: flat_points)
            ARIADNE_TESTING_EXPECT_EQUAL(point.point()[0].raw().get_d(), 0.0);
        const UpperBoxType singleton({{0, 0}, {2, 2}});
        const auto one = EvolutionDetail::sample_hybrid_initial_points(
            singleton, loc, space, DiscretizationHybridType::HMince, 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(one.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(one[0].point()[1].raw().get_d(), 2.0);
    }
};
} // namespace

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_simulator_sampling_contract");
    ARIADNE_TESTING_ADD_TEST(TestSimulatorSamplingContract, test_continuous_resolution_and_subspace);
    ARIADNE_TESTING_ADD_TEST(TestSimulatorSamplingContract, test_hybrid_resolution_and_flat_coordinate);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
