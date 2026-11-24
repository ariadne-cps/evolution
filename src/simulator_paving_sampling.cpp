/***************************************************************************
 *            simulator_paving_sampling.cpp
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

/** @file src/simulator_paving_sampling.cpp
 *  @brief Legacy tree-based discretisation of continuous and hybrid initial boxes.
 */

#include "simulator_paving_sampling.hpp"

#include "algebra/algebra.hpp"
#include "paving/binary_tree_paving.hpp"
#include "symbolic/space.hpp"

namespace Ariadne {

namespace {

// These implementations are moved verbatim from the two simulator source files.
// Retain the same Grid origin, cell iteration, and midpoint semantics until a
// grid-independent replacement is proved equivalent.
BinaryTreePaving create_paving(UpperBoxType box, Space<Real> spc, subspace sspc) {
    for (const auto& entry: sspc) {
        ARIADNE_PRECONDITION(spc.contains(entry.first));
    }
    Nat divider = 1u;
    FloatDP eps_box(0.01_q, Ariadne::ROUND_TO_NEAREST, dp);

    // Widen Initial set
    box = widen(box, eps_box);

    // Point<FloatDP> tmpPointCenter = Point(box.centre());
    Point<FloatDP> tmpPointCenter = Point(box.lower_bounds());
    Vector<FloatDP> origin(tmpPointCenter.dimension(), FloatDP(0, dp));
    Vector<FloatDP> lengths(box.widths().size(), FloatDP(0, dp));

    // Widening a finite coordinate by eps_box>0 gives it positive width.
    // In zero dimensions, lower_bounds() and midpoint() are both empty points.
    for (SizeType i = 0; i < tmpPointCenter.dimension(); i++) {
        divider = 1;
        if (sspc.has_key(spc[i])) {
            Nat val = sspc.get(spc[i]);
            divider = val;
        }
        lengths[i] = (box[i].width().raw() / divider).value();
        origin[i] = tmpPointCenter[i];
    }
    Grid grid(origin, lengths);
    BinaryTreePaving gridPaving(grid);

    return gridPaving;
}

List<LabelledPoint<FloatDPApproximation>> create_point_list(BinaryTreePaving& paving, RealSpace spc, DiscretisationType discretisation_type, Nat mince_dimension) {
    if (discretisation_type == DiscretisationType::Mince)
        paving.mince(mince_dimension);
    else if (discretisation_type == DiscretisationType::Recombine)
        paving.recombine();

    BinaryTreePaving::ConstIterator iter = paving.begin();
    List<LabelledPoint<FloatDPApproximation>> result(paving.size(), LabelledPoint(spc, Point<FloatDPApproximation>()));
    SizeType k(0);
    for (; iter != paving.end(); ++iter) {
        UpperBoxType cell = iter->box();
        auto midpoint = cell.centre();
        LabelledPoint<FloatDPApproximation> pt(spc, midpoint);
        result.at(k) = LabelledPoint<FloatDPApproximation>(spc, pt);
        k++;
    }

    return result;
}

Vector<HybridPoint<FloatDPApproximation>> create_hybrid_point_list(BinaryTreePaving& grid, DiscreteLocation loc, RealSpace spc, DiscretizationHybridType _dtype, Nat fineness, const UpperBoxType& sampling_box) {
    if (_dtype == DiscretizationHybridType::HMince) {
        grid.mince(fineness);
    } else {
        grid.recombine();
    }

    BinaryTreePaving::ConstIterator iter = grid.begin();
    Vector<HybridPoint<FloatDPApproximation>> approxPointList(grid.size(), HybridPoint<FloatDPApproximation>());
    SizeType k(0);
    for (; iter != grid.end(); ++iter) {
        // Outer approximation may include cells crossing the requested boundary.
        // Sample their intersection so every initial point belongs to the box.
        auto sampled = intersection(UpperBoxType(iter->box()), sampling_box);
        auto midpoint = sampled.midpoint();
        HybridPoint<FloatDPApproximation> pt(loc, spc, midpoint);
        approxPointList.at(k) = pt;
        k++;
    }

    return approxPointList;
}

BinaryTreePaving create_hybrid_grid(UpperBoxType box) {
    Point<FloatDP> tmpPointCenter = Point(box.centre());
    Vector<FloatDP> origin(tmpPointCenter.dimension(), FloatDP(0, dp));
    Vector<FloatDP> lengths(box.widths().size(), FloatDP(0, dp));
    for (SizeType i = 0; i < tmpPointCenter.dimension(); i++) {
        origin[i] = tmpPointCenter[i];
        lengths[i] = box.widths().at(i).raw();
    }
    Grid grid(origin, lengths);
    BinaryTreePaving gridPaving(grid);

    return gridPaving;
}

} // namespace

namespace EvolutionDetail {

VectorFieldSimulator::ApproximateListPointType sample_vector_field_initial_points(
    const UpperBoxType& initial_box, const RealSpace& space,
    const subspace& selected_subspace, DiscretisationType discretisation,
    Nat subdivisions) {
    BinaryTreePaving paving = create_paving(initial_box, space, selected_subspace);
    paving.adjoin_outer_approximation(initial_box, subdivisions);
    return create_point_list(paving, space, discretisation, subdivisions);
}

HybridSimulator::HybridApproximateListPointType sample_hybrid_initial_points(
    const UpperBoxType& box, const UpperBoxType& sampling_box,
    DiscreteLocation location, const RealSpace& space,
    DiscretizationHybridType discretisation, Nat fineness) {
    BinaryTreePaving paving = create_hybrid_grid(box);
    paving.adjoin_outer_approximation(box, fineness);
    return create_hybrid_point_list(paving, location, space,
                                    discretisation, fineness, sampling_box);
}

} // namespace EvolutionDetail
} // namespace Ariadne
