/***************************************************************************
 *            discrete_time_approximate_evolver.hpp
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

/** @file evolution/discrete_time_approximate_evolver.hpp
 *  @brief Approximate point trajectories of difference equations and bounded-input inclusions.
 */

#ifndef ARIADNE_DISCRETE_TIME_APPROXIMATE_EVOLVER_HPP
#define ARIADNE_DISCRETE_TIME_APPROXIMATE_EVOLVER_HPP

#include <memory>
#include <utility>

#include "evolution/configuration_interface.hpp"
#include "evolution/discrete_dynamics.hpp"
#include "evolution/evolver.decl.hpp"
#include "evolution/evolver_interface.hpp"
#include "evolution/input_sampling.hpp"
#include "evolution/orbit.hpp"
#include "evolution/sampling_counts.hpp"
#include "geometry/box.hpp"
#include "labelled/point.hpp"
#include "utility/macros.hpp"

namespace Ariadne {

/** @brief Initial-box subdivision counts for approximate discrete evolution.
 *  @details One midpoint is sampled by default; selected nonflat axes may
 *  produce more points. No paving or validated-set semantics are involved.
 */
template<>
class DiscreteTimeEvolverConfiguration<ApproximateTag> : public ConfigurationInterface {
  public:
    using Paradigm = ApproximateTag;
    DiscreteTimeEvolverConfiguration() = default;

    Void set_sampling_counts(SamplingCounts counts) {
        for (const auto& entry: counts)
            ARIADNE_PRECONDITION(entry.second > 0u);
        _sampling_counts = std::move(counts);
    }
    Void set_sampling_count(RealVariable variable, Nat count) {
        ARIADNE_PRECONDITION(count > 0u);
        ARIADNE_PRECONDITION(!_sampling_counts.has_key(variable));
        _sampling_counts.insert(variable, count);
    }
    Void clear_sampling_counts() { _sampling_counts.remove_keys(_sampling_counts.keys()); }
    const SamplingCounts& sampling_counts() const { return _sampling_counts; }

    OutputStream& _write(OutputStream& os) const override;

  private:
    SamplingCounts _sampling_counts;
};

/** @brief Approximate integer-tick trajectories of a deterministic difference map.
 *  @details Explicit initial points are never resampled. For a box, sampling
 *  uses the paving-independent midpoint policy shared with the continuous
 *  approximate evolvers. Inclusion inputs require a separate explicit policy
 *  and are sampled anew at each integer index, independently for each path.
 *  This class has no UPPER or LOWER set semantics.
 */
template<>
class DiscreteTimeEvolver<ApproximateTag>
    : public EvolverInterface<ApproximateTag, DiscreteDynamics,
                              List<LabelledPoint<FloatDPApproximation>>, Integer> {
  public:
    using Paradigm = ApproximateTag;
    using SystemType = DiscreteDynamics;
    using TerminationType = Integer;
    using TimeType = Integer;
    using ApproximatePointType = LabelledPoint<FloatDPApproximation>;
    using InitialType = List<ApproximatePointType>;
    using Interface = EvolverInterface<ApproximateTag, SystemType, InitialType, TerminationType>;
    using OrbitType = Orbit<InitialType>;
    using ResultType = OrbitType;
    using ConfigurationType = DiscreteTimeEvolverConfiguration<ApproximateTag>;

    explicit DiscreteTimeEvolver(const DifferenceAlgebraicEquation& system);
    /** @brief A bounded-input inclusion requires an explicit sampling policy. */
    DiscreteTimeEvolver(const DifferenceAlgebraicInclusion& system, const InputSamplingPolicy& policy);
    DiscreteTimeEvolver(const DiscreteTimeEvolver& other);

    DiscreteTimeEvolver* clone() const override { return new DiscreteTimeEvolver(*this); }
    const SystemType& system() const override { return *_system; }
    ConfigurationType& configuration() { return *_configuration; }
    const ConfigurationType& configuration() const { return *_configuration; }

    /** @brief Iterate canonicalized, labelled initial points for the given ticks. */
    OrbitType orbit(const InitialType&, const TerminationType&) const override;
    /** @brief Sample a box once using the named initial sampling counts. */
    OrbitType orbit(const UpperBoxType&, const TerminationType&) const;
    OrbitType orbit(const ExactBoxType&, const TerminationType&) const;

  private:
    /** @brief One sample per trajectory and source tick, including auxiliaries. */
    std::shared_ptr<const InputSamplingPolicy> _input_sampling;
    std::shared_ptr<const SystemType> _system;
    std::shared_ptr<ConfigurationType> _configuration;
};

} // namespace Ariadne

#endif // ARIADNE_DISCRETE_TIME_APPROXIMATE_EVOLVER_HPP
