/***************************************************************************
 *            src/simulator_sampling.cpp
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

/** @file src/simulator_sampling.cpp
 *  @brief Box-domain sampling for continuous and hybrid evolution.
 */

#include "simulator_sampling.hpp"

#include <vector>

#include "domain_bisection_sampling.hpp"
#include "symbolic/space.hpp"
#include "utility/macros.hpp"

namespace Ariadne::EvolutionDetail {

VectorFieldSimulator::ApproximateListPointType sample_vector_field_initial_points(
    const UpperBoxType& initial_box, const RealSpace& space,
    const subspace& selected_subspace, DiscretisationType discretisation,
    Nat subdivisions) {
    ARIADNE_PRECONDITION(initial_box.dimension() == space.dimension());
    for (const auto& entry: selected_subspace) {
        ARIADNE_PRECONDITION(space.contains(entry.first));
        ARIADNE_PRECONDITION(entry.second > 0u);
    }

    const Nat default_count =
        discretisation == DiscretisationType::None ? 1u :
        bisection_count(subdivisions,
                        discretisation == DiscretisationType::Recombine);
    std::vector<Nat> counts(space.dimension(), 1u);
    const auto widths = initial_box.widths();
    for (SizeType axis = 0u; axis < space.dimension(); ++axis) {
        if (selected_subspace.has_key(space[axis]))
            counts[axis] = selected_subspace.get(space[axis]);
        else if (widths[axis].get_d() > 0.0)
            counts[axis] = default_count;
    }

    const auto domains = bisect_domain(initial_box, counts);
    VectorFieldSimulator::ApproximateListPointType result(
        domains.size(), LabelledPoint<FloatDPApproximation>(
            space, Point<FloatDPApproximation>()));
    for (SizeType i = 0u; i < domains.size(); ++i)
        result.at(i) = LabelledPoint<FloatDPApproximation>(
            space, domains[i].centre());
    return result;
}

HybridSimulator::HybridApproximateListPointType sample_hybrid_initial_points(
    const UpperBoxType& domain, DiscreteLocation location,
    const RealSpace& space, DiscretizationHybridType discretisation,
    Nat fineness) {
    ARIADNE_PRECONDITION(domain.dimension() == space.dimension());
    const Nat count = bisection_count(
        fineness, discretisation == DiscretizationHybridType::HRecombine);
    const auto widths = domain.widths();
    std::vector<Nat> counts(space.dimension(), 1u);
    for (SizeType axis = 0u; axis < space.dimension(); ++axis)
        if (widths[axis].get_d() > 0.0)
            counts[axis] = count;

    const auto domains = bisect_domain(domain, counts);
    HybridSimulator::HybridApproximateListPointType result(
        domains.size(), HybridSimulator::HybridApproximatePointType());
    for (SizeType i = 0u; i < domains.size(); ++i)
        result.at(i) = HybridSimulator::HybridApproximatePointType(
            location, space, domains[i].midpoint());
    return result;
}

} // namespace Ariadne::EvolutionDetail
