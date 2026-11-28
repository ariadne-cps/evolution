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

namespace {

/** @brief Map named counts into domain-axis order.
 * Flat axes are sampled once to avoid duplicate points.
 */
std::vector<Nat> ordered_sampling_counts(
    const UpperBoxType& domain, const RealSpace& space,
    const SamplingCounts& selected) {
    ARIADNE_PRECONDITION(domain.dimension() == space.dimension());
    for (const auto& entry: selected) {
        ARIADNE_PRECONDITION(space.contains(entry.first));
        ARIADNE_PRECONDITION(entry.second > 0u);
    }
    std::vector<Nat> counts(space.dimension(), 1u);
    const auto widths = domain.widths();
    for (SizeType axis = 0u; axis < space.dimension(); ++axis)
        if (widths[axis].get_d() > 0.0 && selected.has_key(space[axis]))
            counts[axis] = selected.get(space[axis]);
    return counts;
}

} // namespace

VectorFieldSimulator::ApproximateListPointType sample_vector_field_initial_points(
    const UpperBoxType& domain, const RealSpace& space,
    const SamplingCounts& counts) {
    const auto domains = bisect_domain(domain, ordered_sampling_counts(domain, space, counts));
    VectorFieldSimulator::ApproximateListPointType result(
        domains.size(), LabelledPoint<FloatDPApproximation>(
                            space, Point<FloatDPApproximation>()));
    for (SizeType i = 0u; i < domains.size(); ++i)
        result.at(i) = LabelledPoint<FloatDPApproximation>(space, domains[i].centre());
    return result;
}

HybridSimulator::HybridApproximateListPointType sample_hybrid_initial_points(
    const UpperBoxType& domain, DiscreteLocation location,
    const RealSpace& space, const SamplingCounts& counts) {
    const auto domains = bisect_domain(domain, ordered_sampling_counts(domain, space, counts));
    HybridSimulator::HybridApproximateListPointType result(
        domains.size(), HybridSimulator::HybridApproximatePointType());
    for (SizeType i = 0u; i < domains.size(); ++i)
        result.at(i) = HybridSimulator::HybridApproximatePointType(
            location, space, domains[i].midpoint());
    return result;
}

} // namespace Ariadne::EvolutionDetail
