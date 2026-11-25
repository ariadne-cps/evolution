/***************************************************************************
 *            discrete_time_inclusion_evolver.hpp
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

/** @file evolution/discrete_time_inclusion_evolver.hpp
 *  @brief Validated UPPER evolution of ordinary bounded-input difference inclusions.
 */

#ifndef ARIADNE_DISCRETE_TIME_INCLUSION_EVOLVER_HPP
#define ARIADNE_DISCRETE_TIME_INCLUSION_EVOLVER_HPP

#include <memory>

#include "evolution/discrete_time_evolver.hpp"
#include "evolution/orbit.hpp"
#include "system/difference_algebraic_inclusion.hpp"

namespace Ariadne {

/** @brief A validated ordinary discrete inclusion evolver.
 *  @details Each progression universally encloses every input from the
 *  declared bounded box, adding a fresh independent INPUT coordinate per
 *  tick. Existing parameter correlations and constraints are retained.
 *  Only UPPER is supported. Explicit auxiliary outputs a[n]=g(x[n],u[n])
 *  bind to the current full input box at each tick; an update consumes the
 *  current u[n] and installs a fresh independent u[n+1] for new outputs.
 *  No input coordinates are sampled or frozen across different ticks.
 *
 *  This class is separate from the equation-only DiscreteTimeEvolver<ValidatedTag>
 *  to preserve that class's existing strongly typed SystemType contract.
 */
class DiscreteTimeInclusionEvolver
    : public EvolverInterface<ValidatedTag, DifferenceAlgebraicInclusion, LabelledEnclosure, Integer> {
  public:
    using Paradigm = ValidatedTag;
    using SystemType = DifferenceAlgebraicInclusion;
    using TimeType = Integer;
    using TerminationType = Integer;
    using EnclosureType = LabelledEnclosure;
    using ConfigurationType = DiscreteTimeEvolverConfiguration<ValidatedTag>;
    using FunctionFactoryType = ValidatedFunctionPatchFactory;
    using OrbitType = Orbit<EnclosureType>;
    using Interface = EvolverInterface<ValidatedTag, SystemType, EnclosureType, TerminationType>;

    explicit DiscreteTimeInclusionEvolver(const SystemType& system);
    DiscreteTimeInclusionEvolver(const DiscreteTimeInclusionEvolver& other);
    DiscreteTimeInclusionEvolver* clone() const override {
        return new DiscreteTimeInclusionEvolver(*this);
    }

    const SystemType& system() const override { return *_system; }
    ConfigurationType& configuration() { return *_configuration; }
    const ConfigurationType& configuration() const { return *_configuration; }

    FunctionFactoryType function_factory() const;
    EnclosureType enclosure(const ExactBoxType&) const;

    OrbitType orbit(const EnclosureType&, const TerminationType&, Semantics) const override;
    OrbitType orbit(const RealExpressionBoundedConstraintSet&, const TerminationType&, Semantics) const;

  private:
    std::shared_ptr<const SystemType> _system;
    std::shared_ptr<ConfigurationType> _configuration;
};

} // namespace Ariadne
#endif // ARIADNE_DISCRETE_TIME_INCLUSION_EVOLVER_HPP
