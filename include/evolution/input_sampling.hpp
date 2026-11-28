/***************************************************************************
 *            evolution/input_sampling.hpp
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

/** @file evolution/input_sampling.hpp
 *  @brief Explicit reproducible input policies for approximate inclusion evolution.
 */

#ifndef ARIADNE_INPUT_SAMPLING_HPP
#define ARIADNE_INPUT_SAMPLING_HPP
#include <cstdint>

#include "algebra/vector.hpp"
#include "geometry/box.hpp"
#include "numeric/floatdp.hpp"
namespace Ariadne {
/** @brief One accepted continuous step; indices are independent of scheduling. */
struct InputSamplingContext {
    SizeType trajectory;
    SizeType step;
    FloatDPApproximation time;
    FloatDPApproximation duration;
};
/** @brief Select numerical provider inputs once per RK4 step.
 *  Implementations must be const, thread-safe and deterministic for a context.
 *  The domain is the structured provider's normalized box, not a physical
 *  input valuation. Returned values must belong to this box.
 */
class InputSamplingPolicy {
  public:
    virtual ~InputSamplingPolicy() = default;
    virtual InputSamplingPolicy* clone() const = 0;
    virtual Vector<FloatDPApproximation> sample(const ExactBoxType&, const InputSamplingContext&) const = 0;
};
/** @brief Explicitly hold the center of each numerical input interval. */
class MidpointInputSampling final : public InputSamplingPolicy {
  public:
    MidpointInputSampling* clone() const override { return new MidpointInputSampling(*this); }
    Vector<FloatDPApproximation> sample(const ExactBoxType&, const InputSamplingContext&) const override;
};
/** @brief Counter-based uniform samples, held constant during each step. */
class UniformInputSampling final : public InputSamplingPolicy {
    std::uint64_t _seed;

  public:
    explicit UniformInputSampling(std::uint64_t seed = 0u) : _seed(seed) {}
    std::uint64_t seed() const { return _seed; }
    UniformInputSampling* clone() const override { return new UniformInputSampling(*this); }
    Vector<FloatDPApproximation> sample(const ExactBoxType&, const InputSamplingContext&) const override;
};
} // namespace Ariadne
#endif
