/***************************************************************************
 *            test_discrete_time_inclusion_evolver.cpp
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

/** @file tests/test_discrete_time_inclusion_evolver.cpp
 *  @brief Witness-independent UPPER propagation of independent input boxes.
 */

#include <memory>
#include <stdexcept>

#include "evolution/discrete_time_inclusion_evolver.hpp"
#include "logging/test_runner.hpp"
#include "system/difference_algebraic_inclusion.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;

class TestDiscreteTimeInclusionEvolver : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_full_box_each_tick_and_copy() {
        const RealVariable x("x"), u("u");
        const DifferenceAlgebraicInclusion inclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            RealVariablesBox({-1 <= u <= 1}));
        DiscreteTimeInclusionEvolver evolver(inclusion);
        const auto initial = evolver.enclosure(ExactBoxType({{0.0_x, 0.0_x}}));

        const auto zero = evolver.orbit(initial, Integer(0), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.reach().size(), 0u);

        const auto result = evolver.orbit(initial, Integer(2), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(result.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.reach().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.intermediate().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.final()[0u].number_of_parameters(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.final()[0u].variable_kinds()[1u],
                                     EnclosureVariableKind::INPUT);
        ARIADNE_TESTING_EXPECT_EQUAL(result.final()[0u].variable_kinds()[2u],
                                     EnclosureVariableKind::INPUT);
        const auto extent = result.final()[0u].euclidean_set().bounding_box()[0u];
        // The independent input choices (+1,+1) and (-1,-1) must both be covered.
        ARIADNE_TESTING_PROVE(extent.lower_bound().raw().get_d() <= -2.0);
        ARIADNE_TESTING_PROVE(extent.upper_bound().raw().get_d() >= 2.0);
        ARIADNE_TESTING_EXPECT_EQUAL(initial.number_of_parameters(), 1u);

        std::unique_ptr<DiscreteTimeInclusionEvolver> copy(evolver.clone());
        copy->configuration().set_enable_subdivisions(true);
        ARIADNE_TESTING_PROVE(!evolver.configuration().enable_subdivisions());
        ARIADNE_TESTING_EXPECT_EQUAL(copy->system().state_space(), RealSpace({x}));
    }

    Void test_rejection_is_explicit_even_at_zero_horizon() {
        const RealVariable x("x"), y("y"), u("u"), a("a");
        const DifferenceAlgebraicInclusion inclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            RealVariablesBox({-1 <= u <= 1}));
        DiscreteTimeInclusionEvolver evolver(inclusion);
        const auto initial = evolver.enclosure(ExactBoxType({{0.0_x, 1.0_x}}));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(initial, Integer(0), Semantics::LOWER));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(initial, Integer(1), Semantics::LOWER));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(initial, Integer(-1), Semantics::UPPER));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.enclosure(ExactBoxType({{0.0_x, 1.0_x}, {0.0_x, 1.0_x}})));
        const LabelledEnclosure wrong(ExactBoxType({{0.0_x, 1.0_x}}),
                                      RealSpace({y}),
                                      EnclosureConfiguration(evolver.function_factory()));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            evolver.orbit(wrong, Integer(0), Semantics::UPPER));

        const DifferenceAlgebraicInclusion auxiliaries(
            List<PrimedRealAssignment>{next(x) = x + u},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1}));
        DiscreteTimeInclusionEvolver unsupported(auxiliaries);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            unsupported.orbit(unsupported.enclosure(ExactBoxType({{0.0_x, 1.0_x}})),
                              Integer(0), Semantics::UPPER));
    }
};

int main(int argc, const char* argv[]) {
    ARIADNE_TESTING_SUITE("test_discrete_time_inclusion_evolver");
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeInclusionEvolver, test_full_box_each_tick_and_copy);
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeInclusionEvolver, test_rejection_is_explicit_even_at_zero_horizon);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
