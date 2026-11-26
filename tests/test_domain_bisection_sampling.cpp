/***************************************************************************
 *            tests/test_domain_bisection_sampling.cpp
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

/** @file tests/test_domain_bisection_sampling.cpp
 *  @brief Verify deterministic bisection without paving.
 */

#include "evolution/domain_bisection_sampling.hpp"
#include "logging/test_runner.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {
class TestDomainBisectionSampling : public Testing::ContextFixture {
  public:
    using Testing::ContextFixture::ContextFixture;

    Void test_two_axis_order_and_coordinates() {
        const UpperBoxType domain({
            {0, 1},
            {0, 1}
        });
        const auto pieces = EvolutionDetail::bisect_domain(domain, {2u, 2u});
        ARIADNE_TESTING_EXPECT_EQUAL(pieces.size(), 4u);
        const double expected[4][2] = {
            {0.25, 0.25},
            {0.25, 0.75},
            {0.75, 0.25},
            {0.75, 0.75}
        };
        for (SizeType i = 0u; i < pieces.size(); ++i) {
            const auto point = pieces[i].midpoint();
            ARIADNE_TESTING_EXPECT_EQUAL(point[0].raw().get_d(), expected[i][0]);
            ARIADNE_TESTING_EXPECT_EQUAL(point[1].raw().get_d(), expected[i][1]);
        }
    }

    Void test_unequal_counts_and_flat_axis() {
        const UpperBoxType domain({
            {0, 1},
            {2, 2}
        });
        const auto pieces = EvolutionDetail::bisect_domain(domain, {3u, 1u});
        ARIADNE_TESTING_EXPECT_EQUAL(pieces.size(), 3u);
        const double expected[3] = {0.25, 0.625, 0.875};
        for (SizeType i = 0u; i < pieces.size(); ++i) {
            const auto point = pieces[i].midpoint();
            ARIADNE_TESTING_EXPECT_EQUAL(point[0].raw().get_d(), expected[i]);
            ARIADNE_TESTING_EXPECT_EQUAL(point[1].raw().get_d(), 2.0);
        }
        ARIADNE_TESTING_EXPECT_EXCEPTION(
            EvolutionDetail::bisect_domain(domain, {2u}));
        ARIADNE_TESTING_EXPECT_EXCEPTION(
            EvolutionDetail::bisect_domain(domain, {0u, 1u}));
    }
};
} // namespace

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_domain_bisection_sampling");
    ARIADNE_TESTING_ADD_TEST(TestDomainBisectionSampling, test_two_axis_order_and_coordinates);
    ARIADNE_TESTING_ADD_TEST(TestDomainBisectionSampling, test_unequal_counts_and_flat_axis);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
