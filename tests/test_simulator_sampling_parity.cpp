/***************************************************************************
 *            test_simulator_sampling_parity.cpp
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

/** @file tests/test_simulator_sampling_parity.cpp
 *  @brief Compare simulator point sampling with its pre-extraction paving path.
 */

#include "algebra/algebra.hpp"
#include "evolution/hybrid_simulator.hpp"
#include "evolution/vector_field_simulator.hpp"
#include "geometry/box.hpp"
#include "logging/test_runner.hpp"
#include "paving/binary_tree_paving.hpp"
#include "simulator_paving_sampling.hpp"
#include "symbolic/space.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

class TestSimulatorSamplingParity : public Testing::ContextFixture {
  public:
    using Testing::ContextFixture::ContextFixture;

    Void check_continuous_sequence(const UpperBoxType& box, const RealSpace& space,
                                   const subspace& selected, DiscretisationType type,
                                   Nat subdivisions) {
        // Independent legacy reference: construct the widened grid, form its
        // outer approximation, and enumerate the cells in backend order.
        const FloatDP epsilon(0.01_q, ROUND_TO_NEAREST, dp);
        const UpperBoxType padded = widen(box, epsilon);
        const Point<FloatDP> origin(padded.lower_bounds());
        Vector<FloatDP> lengths(padded.widths().size(), FloatDP(0, dp));
        for (SizeType axis = 0; axis < space.dimension(); ++axis) {
            const Nat divisor = selected.has_key(space[axis]) ? selected.get(space[axis]) : 1u;
            lengths[axis] = (padded[axis].width().raw() / divisor).value();
        }
        BinaryTreePaving reference(Grid(origin, lengths));
        reference.adjoin_outer_approximation(box, subdivisions);
        if (type == DiscretisationType::Mince)
            reference.mince(subdivisions);
        else if (type == DiscretisationType::Recombine)
            reference.recombine();

        const auto actual = EvolutionDetail::sample_vector_field_initial_points(
            box, space, selected, type, subdivisions);
        ARIADNE_TESTING_EXPECT_EQUAL(actual.size(), reference.size());
        SizeType index = 0u;
        for (const auto& cell: reference) {
            const UpperBoxType expected_box = cell.box();
            const auto expected_point = expected_box.centre();
            ARIADNE_TESTING_PROVE(actual.at(index).space() == space);
            for (SizeType axis = 0; axis < space.dimension(); ++axis) {
                ARIADNE_TESTING_EXPECT_EQUAL(
                    actual.at(index)[axis].raw().get_d(),
                    expected_point[axis].raw().get_d());
            }
            ++index;
        }
        ARIADNE_TESTING_EXPECT_EQUAL(index, actual.size());
    }

    Void test_continuous_subspace_and_order() {
        const RealVariable x("x"), y("y"), missing("missing");
        const RealSpace space({x, y});
        const UpperBoxType box({{0, 1}, {-2, 2}});
        for (const auto& selected: {subspace{}, subspace{{x, 2u}},
                                   subspace{{x, 2u}, {y, 3u}}}) {
            for (auto type: {DiscretisationType::None,
                             DiscretisationType::Mince,
                             DiscretisationType::Recombine}) {
                for (Nat subdivisions: {0u, 1u, 2u}) {
                    check_continuous_sequence(box, space, selected, type,
                                              selected.empty() ? subdivisions : 0u);
                }
            }
        }
        ARIADNE_TESTING_EXPECT_EXCEPTION(
            EvolutionDetail::sample_vector_field_initial_points(
                box, space, subspace{{missing, 2u}},
                DiscretisationType::Mince, 0u));
    }

    Void check_hybrid_sequence(const UpperBoxType& box, const UpperBoxType& original,
                               DiscretizationHybridType type, Nat fineness) {
        const RealVariable x("x"), y("y");
        const RealSpace space({x, y});
        const DiscreteLocation location(1);

        // Independent legacy reference: the hybrid path uses box centre as
        // grid origin, then clips each cell before taking its midpoint.
        const Point<FloatDP> origin(box.centre());
        Vector<FloatDP> lengths(box.widths().size(), FloatDP(0, dp));
        for (SizeType axis = 0; axis < lengths.size(); ++axis)
            lengths[axis] = box.widths().at(axis).raw();
        BinaryTreePaving reference(Grid(origin, lengths));
        reference.adjoin_outer_approximation(box, fineness);
        if (type == DiscretizationHybridType::HMince)
            reference.mince(fineness);
        else
            reference.recombine();

        const auto actual = EvolutionDetail::sample_hybrid_initial_points(
            box, original, location, space, type, fineness);
        ARIADNE_TESTING_EXPECT_EQUAL(actual.size(), reference.size());
        SizeType index = 0u;
        for (const auto& cell: reference) {
            const auto clipped = intersection(UpperBoxType(cell.box()), original);
            const auto expected_point = clipped.midpoint();
            ARIADNE_TESTING_PROVE(actual.at(index).location() == location);
            ARIADNE_TESTING_PROVE(actual.at(index).space() == space);
            for (SizeType axis = 0; axis < space.dimension(); ++axis) {
                ARIADNE_TESTING_EXPECT_EQUAL(
                    actual.at(index).point()[axis].raw().get_d(),
                    expected_point[axis].raw().get_d());
            }
            ++index;
        }
        ARIADNE_TESTING_EXPECT_EQUAL(index, actual.size());
    }

    Void test_hybrid_clipping_and_order() {
        const UpperBoxType box({{0, 2}, {-1, 1}});
        const UpperBoxType sampling_box({{0.1_dec, 1.9_dec},
                                         {-0.8_dec, 0.7_dec}});
        for (auto type: {DiscretizationHybridType::HMince,
                         DiscretizationHybridType::HRecombine}) {
            for (Nat fineness: {0u, 1u, 2u}) {
                check_hybrid_sequence(box, sampling_box, type, fineness);
            }
        }

        // The public hybrid simulator widens a partially flat initial box
        // before calling the sampler. The original box remains the clipping
        // domain, including its zero-width coordinate.
        const UpperBoxType partially_flat({{0, 0}, {1, 2}});
        const UpperBoxType widened = widen(partially_flat,
                                           FloatDPUpperBound(0.0001_q, dp));
        check_hybrid_sequence(widened, partially_flat,
                              DiscretizationHybridType::HMince, 1u);
        check_hybrid_sequence(widened, partially_flat,
                              DiscretizationHybridType::HRecombine, 1u);
    }
};

} // namespace

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_simulator_sampling_parity");
    ARIADNE_TESTING_ADD_TEST(TestSimulatorSamplingParity, test_continuous_subspace_and_order);
    ARIADNE_TESTING_ADD_TEST(TestSimulatorSamplingParity, test_hybrid_clipping_and_order);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
