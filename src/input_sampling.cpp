/***************************************************************************
 *            input_sampling.cpp
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

/** @file src/input_sampling.cpp
 *  @brief Implement deterministic piecewise-constant input sampling.
 */

#include "evolution/input_sampling.hpp"

#include <cmath>
#include <stdexcept>

#include "numeric/rounding.hpp"
namespace Ariadne {
namespace {
std::uint64_t mix(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30u)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27u)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31u);
}
Vector<FloatDPApproximation> sample_box(const ExactBoxType& domain, const InputSamplingContext& context,
                                        std::uint64_t seed, bool midpoint) {
    RoundingSession<FloatDP> rounding_session;
    rounding_session.set(FloatDP::ROUND_TO_NEAREST);
    Vector<FloatDPApproximation> result(domain.dimension(), FloatDPApproximation(0, dp));
    for (SizeType i = 0u; i != domain.dimension(); ++i) {
        const double lower = domain[i].lower_bound().get_d(), upper = domain[i].upper_bound().get_d();
        if (!std::isfinite(lower) || !std::isfinite(upper) || lower > upper)
            throw std::invalid_argument("Input sampling requires a finite nonempty numerical input box");
        const auto bits = mix(seed ^ mix(static_cast<std::uint64_t>(context.trajectory)) ^
                              mix(static_cast<std::uint64_t>(context.step) + 0x517cc1b727220a95ULL) ^
                              mix(static_cast<std::uint64_t>(i) + 0x6eed0e9da4d94a4fULL));
        const double unit = midpoint ? 0.5 : static_cast<double>(bits >> 11u) * 0x1p-53;
        result[i] = FloatDPApproximation(ExactDouble(std::lerp(lower, upper, unit)), dp);
    }
    return result;
}
} // namespace
Vector<FloatDPApproximation> MidpointInputSampling::sample(const ExactBoxType& domain, const InputSamplingContext& context) const {
    return sample_box(domain, context, 0u, true);
}
Vector<FloatDPApproximation> UniformInputSampling::sample(const ExactBoxType& domain, const InputSamplingContext& context) const {
    return sample_box(domain, context, _seed, false);
}
} // namespace Ariadne
