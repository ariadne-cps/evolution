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

    Void test_constrained_multidimensional_inputs_preserve_history() {
        const RealVariable x("x"), u("u"), v("v");
        const DifferenceAlgebraicInclusion inclusion(
            List<PrimedRealAssignment>{next(x) = x + u + v},
            RealVariablesBox({-1 <= u <= 1, 0 <= v <= 2}));
        DiscreteTimeInclusionEvolver evolver(inclusion);
        const auto identity = EffectiveVectorMultivariateFunction::identity(1u);
        const auto coordinate = identity[0u];
        auto initial = evolver.enclosure(ExactBoxType({{-1.0_x, 1.0_x}}));
        initial.new_positive_state_constraint(coordinate - Real(1 / 4_q));
        initial.new_negative_state_constraint(coordinate - Real(1 / 2_q));

        const auto result = evolver.orbit(initial, Integer(2), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(result.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.final()[0u].number_of_parameters(), 5u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.final()[0u].euclidean_set().number_of_constraints(), 2u);
        for (SizeType i = 1u; i != 5u; ++i)
            ARIADNE_TESTING_EXPECT_EQUAL(result.final()[0u].variable_kinds()[i],
                                         EnclosureVariableKind::INPUT);

        // A single feasible source x0=3/8 and distinct extreme input
        // sequences are represented by the same UPPER enclosure.
        // Both tick choices are free; no first-tick input is frozen.
        const auto bounds = result.final()[0u].euclidean_set().bounding_box()[0u];
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= -1.625);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 6.375);
        ARIADNE_TESTING_EXPECT_EQUAL(initial.number_of_parameters(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(initial.euclidean_set().number_of_constraints(), 2u);

        // The public validated interface must retain its own UPPER orbit.
        const DiscreteTimeInclusionEvolver::Interface& generic = evolver;
        const auto interface_result = generic.orbit(initial, Integer(1), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(interface_result.final().size(), 1u);
    }

    Void test_upper_subdivision_keeps_all_successors() {
        const RealVariable x("x"), u("u");
        const DifferenceAlgebraicInclusion inclusion(
            List<PrimedRealAssignment>{next(x) = x + u},
            RealVariablesBox({-1 <= u <= 1}));
        DiscreteTimeInclusionEvolver evolver(inclusion);
        evolver.configuration().set_enable_subdivisions(true);
        evolver.configuration().set_maximum_enclosure_radius(0.75);
        const auto initial = evolver.enclosure(ExactBoxType({{-1.0_x, 1.0_x}}));
        const auto result = evolver.orbit(initial, Integer(1), Semantics::UPPER);
        ARIADNE_TESTING_PROVE(result.final().size() > 1u);
        bool covers_negative = false;
        bool covers_positive = false;
        for (const auto& leaf : result.final()) {
            const auto extent = leaf.euclidean_set().bounding_box()[0u];
            covers_negative = covers_negative ||
                              extent.lower_bound().raw().get_d() <= -2.0;
            covers_positive = covers_positive ||
                              extent.upper_bound().raw().get_d() >= 2.0;
        }
        ARIADNE_TESTING_PROVE(covers_negative);
        ARIADNE_TESTING_PROVE(covers_positive);
        ARIADNE_TESTING_EXPECT_EQUAL(initial.number_of_parameters(), 1u);
    }

    Void test_input_dependent_auxiliary_and_fresh_tick_inputs() {
        const RealVariable x("x"), u("u"), a("a");
        const DifferenceAlgebraicInclusion inclusion(
            List<PrimedRealAssignment>{next(x) = x + a},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1}));
        DiscreteTimeInclusionEvolver evolver(inclusion);
        const auto initial = evolver.enclosure(ExactBoxType({{0.0_x, 0.0_x}}));

        // u[0] is a free coordinate already at horizon zero: a[0]=u[0].
        const auto zero = evolver.orbit(initial, Integer(0), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.initial().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.initial().auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(zero.initial().auxiliary_arguments(),
                                     EnclosureFunctionArguments::STATE_INPUTS);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.initial().variable_kinds()[1u],
                                     EnclosureVariableKind::INPUT);
        const auto zero_aux = zero.initial().auxiliary_function()[0u].range();
        ARIADNE_TESTING_PROVE(zero_aux.lower_bound().raw().get_d() <= -1.0);
        ARIADNE_TESTING_PROVE(zero_aux.upper_bound().raw().get_d() >= 1.0);

        // x[1]=u[0], x[2]=2*u[0]+u[1] and a[2]=x[2]+u[2].
        // The three inputs are universally quantified independently.
        const auto evolved = evolver.orbit(initial, Integer(2), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.final()[0u].number_of_parameters(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(evolved.final()[0u].auxiliary_space(), RealSpace({a}));
        for (SizeType i = 1u; i != 4u; ++i)
            ARIADNE_TESTING_EXPECT_EQUAL(evolved.final()[0u].variable_kinds()[i],
                                         EnclosureVariableKind::INPUT);
        const auto bounds = evolved.final()[0u].euclidean_set().bounding_box()[0u];
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= -3.0);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 3.0);
        const auto aux = evolved.final()[0u].auxiliary_function();
        const auto current_input = evolved.final()[0u].configuration()
            .function_factory().create_coordinate(evolved.final()[0u].parameter_domain(), 3u);
        const auto relation = aux[0u] - evolved.final()[0u].state_function()[0u] - current_input;
        const auto error = relation.range();
        ARIADNE_TESTING_PROVE(error.lower_bound().raw().get_d() >= -1e-6);
        ARIADNE_TESTING_PROVE(error.upper_bound().raw().get_d() <= 1e-6);
        ARIADNE_TESTING_EXPECT_EQUAL(initial.number_of_parameters(), 1u);

        // Reuse the last full current-input box without silently resampling
        // a[2]; progression consumes u[2] and allocates a new u[3].
        const auto resumed = evolver.orbit(evolved.final()[0u], Integer(1),
                                            Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(resumed.final().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(resumed.final()[0u].number_of_parameters(), 5u);
        ARIADNE_TESTING_EXPECT_EQUAL(resumed.final()[0u].auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(resumed.final()[0u].variable_kinds()[4u],
                                     EnclosureVariableKind::INPUT);
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
            List<PrimedRealAssignment>{next(x) = x + a},
            List<RealAssignment>{let(a) = x + u},
            RealVariablesBox({-1 <= u <= 1}));
        DiscreteTimeInclusionEvolver supported(auxiliaries);
        const auto zero = supported.orbit(
            supported.enclosure(ExactBoxType({{0.0_x, 0.0_x}})),
            Integer(0), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(zero.initial().auxiliary_space(), RealSpace({a}));
        ARIADNE_TESTING_EXPECT_EQUAL(zero.initial().number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            supported.orbit(supported.enclosure(ExactBoxType({{0.0_x, 0.0_x}})),
                            Integer(0), Semantics::LOWER));
    }
};

int main(int argc, const char* argv[]) {
    ARIADNE_TESTING_SUITE("test_discrete_time_inclusion_evolver");
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeInclusionEvolver, test_full_box_each_tick_and_copy);
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeInclusionEvolver, test_constrained_multidimensional_inputs_preserve_history);
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeInclusionEvolver, test_upper_subdivision_keeps_all_successors);
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeInclusionEvolver, test_input_dependent_auxiliary_and_fresh_tick_inputs);
    ARIADNE_TESTING_ADD_TEST(TestDiscreteTimeInclusionEvolver, test_rejection_is_explicit_even_at_zero_horizon);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
