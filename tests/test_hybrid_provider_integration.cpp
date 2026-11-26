/***************************************************************************
 *            test_hybrid_provider_integration.cpp
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

/** @file tests/test_hybrid_provider_integration.cpp
 *  @brief Exercise symbolic initial sets with the integrated Hybrid provider.
 */

#include <memory>

#include "evolution/hybrid_evolver.hpp"
#include "hybrid/hybrid_expression_set.hpp"
#include "logging/test_runner.hpp"
#include "system/hybrid_automaton.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

class TestHybridProviderIntegration : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    void test_symbolic_initial_sets_and_location_spaces() {
        const RealVariable x("x"), y("y");
        const DiscreteLocation first(1), second(2);
        HybridAutomaton system;
        system.new_mode(first, {dot(x) = 1});
        system.new_mode(second, {dot(y) = 0, dot(x) = 0});
        GeneralHybridEvolver evolver(system);

        const HybridBoxSet box(first, {0 <= x <= 0});
        const HybridBoundedConstraintSet constrained(first, {0 <= x <= 0});
        const auto from_box = evolver.enclosure(
            HybridExactBox(first, RealSpace({x}), ExactBoxType({{0, 0}})));
        const auto from_constraint = evolver.enclosure(constrained);
        ARIADNE_TESTING_PROVE(from_constraint.state_space() == RealSpace({x}));
        ARIADNE_TESTING_EXPECT_EQUAL(from_constraint.location(), first);
        ARIADNE_TESTING_EXPECT_EQUAL(from_constraint.state_bounding_box(), from_box.state_bounding_box());
        std::unique_ptr<HybridEnclosure> copy(from_constraint.clone());
        ARIADNE_TESTING_EXPECT_EQUAL(copy->state_bounding_box(), from_constraint.state_bounding_box());

        const HybridBoundedConstraintSet other(second, {1 <= x <= 1, 2 <= y <= 2});
        const auto other_enclosure = evolver.enclosure(other);
        const RealSpace ordered = system.continuous_state_space(second);
        ARIADNE_TESTING_PROVE(other_enclosure.state_space() == ordered);
        ARIADNE_TESTING_EXPECT_EQUAL(other_enclosure.state_bounding_box().dimension(), 2u);
        const ExactBoxType expected = ordered[0] == x
            ? ExactBoxType({{0.999_x, 1.001_x}, {1.999_x, 2.001_x}})
            : ExactBoxType({{1.999_x, 2.001_x}, {0.999_x, 1.001_x}});
        ARIADNE_TESTING_PROVE(definitely(other_enclosure.state_bounding_box().inside(
            expected)));

        evolver.configuration().set_maximum_step_size(0.125_x);
        const auto orbit = evolver.orbit(box, HybridTime(0.125_x, 1), Semantics::UPPER);
        ARIADNE_TESTING_PROVE(orbit.final().size() > 0u);
        for (const auto& final_set: orbit.final()) {
            ARIADNE_TESTING_EXPECT_EQUAL(final_set.location(), first);
            ARIADNE_TESTING_PROVE(final_set.state_space() == RealSpace({x}));
            ARIADNE_TESTING_PROVE(definitely(final_set.state_bounding_box().inside(
                ExactBoxType({{0.124_x, 0.126_x}}))));
        }

        const auto stopped = evolver.orbit(constrained, HybridTime(0.0_x, 1), Semantics::UPPER);
        ARIADNE_TESTING_PROVE(stopped.final().size() > 0u);
        for (const auto& final_set: stopped.final()) {
            ARIADNE_TESTING_EXPECT_EQUAL(final_set.location(), first);
            ARIADNE_TESTING_EXPECT_EQUAL(final_set.state_bounding_box(), from_constraint.state_bounding_box());
        }
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_provider_integration");
    ARIADNE_TESTING_ADD_TEST(TestHybridProviderIntegration, test_symbolic_initial_sets_and_location_spaces);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
