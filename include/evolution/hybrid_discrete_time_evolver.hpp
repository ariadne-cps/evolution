/***************************************************************************
 *            evolution/hybrid_discrete_time_evolver.hpp
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

/** @file evolution/hybrid_discrete_time_evolver.hpp
 *  @brief Validated UPPER evolution of deterministic discrete hybrid automata.
 */

#ifndef ARIADNE_HYBRID_DISCRETE_TIME_EVOLVER_HPP
#define ARIADNE_HYBRID_DISCRETE_TIME_EVOLVER_HPP

#include <memory>
#include "evolution/configuration_interface.hpp"
#include "evolution/evolver.decl.hpp"
#include "evolution/evolver_interface.hpp"
#include "evolution/hybrid_discrete_time_upper_scheduler.hpp"
#include "system/discrete_time_hybrid_automaton_interface.hpp"

namespace Ariadne {

//! \brief Exact absolute integer indices of an initial enclosure and horizon.
//! \details Avoid inferring an integer index from a validated real-time range.
class HybridDiscreteTimeTerminationCriterion {
  public:
    HybridDiscreteTimeTerminationCriterion(Integer initial_index, Integer termination_index);
    const Integer& initial_index() const { return _initial_index; }
    const Integer& termination_index() const { return _termination_index; }
  private:
    Integer _initial_index;
    Integer _termination_index;
};

//! \brief Hard node budget; exhausting it aborts an incomplete UPPER exploration.
template<>
class HybridDiscreteTimeEvolverConfiguration<ValidatedTag> : public ConfigurationInterface {
  public:
    using Paradigm = ValidatedTag;
    HybridDiscreteTimeEvolverConfiguration();
    SizeType maximum_nodes() const { return _maximum_nodes; }
    Void set_maximum_nodes(SizeType value);
    OutputStream& _write(OutputStream& os) const override;
  private:
    SizeType _maximum_nodes;
};

//! \brief Public UPPER subset of validated hybrid discrete-time evolution.
//! \details All possible certified branches consume one tick; uncertain
//! predicates retain alternatives. No branch priority, time-zero jump,
//! implicit merging is imposed. Bounded-input inclusion progressions and
//! transitions to equation or inclusion modes are supported with or
//! without algebraic auxiliary outputs. Guard, invariant and reset
//! expressions may refer to these outputs via the correlated enclosure.
//! A transition into an auxiliary-bearing inclusion target installs a
//! fresh independent arrival input only after applying the source reset. Full input
//! boxes are attached at the current index; each progression introduces
//! a new independent input box.
//! LOWER is unsupported, even if a branch precondition is certain.
//! A blocked leaf is not an orbit.final() enclosure: retrieve it via
//! orbit_with_blocked() or inspect provenance through explore().
template<>
class HybridDiscreteTimeEvolver<ValidatedTag>
    : public EvolverInterface<ValidatedTag, DiscreteTimeHybridAutomatonInterface,
                              HybridEnclosure, HybridDiscreteTimeTerminationCriterion> {
  public:
    using Paradigm = ValidatedTag;
    using Interface = EvolverInterface<ValidatedTag, DiscreteTimeHybridAutomatonInterface,
                                      HybridEnclosure, HybridDiscreteTimeTerminationCriterion>;
    using SystemType = DiscreteTimeHybridAutomatonInterface;
    using ConfigurationType = HybridDiscreteTimeEvolverConfiguration<ValidatedTag>;
    using TimeType = Integer;
    using TerminationType = HybridDiscreteTimeTerminationCriterion;
    using EnclosureType = HybridEnclosure;
    using InitialType = EnclosureType;
    using OrbitType = Orbit<EnclosureType>;
    using ResultType = OrbitType;

    explicit HybridDiscreteTimeEvolver(const SystemType& system);
    HybridDiscreteTimeEvolver(const HybridDiscreteTimeEvolver& other);
    HybridDiscreteTimeEvolver<ValidatedTag>* clone() const override;
    const SystemType& system() const override { return *_system; }
    ConfigurationType& configuration() { return *_configuration; }
    const ConfigurationType& configuration() const { return *_configuration; }

    OrbitType orbit(const EnclosureType& initial, const TerminationType& termination,
                    Semantics semantics) const override;
    //! \brief Preserve blocked leaves separately from horizon-only orbit.final().
    HybridDiscreteTimeUpperOrbitProjection orbit_with_blocked(
        const EnclosureType& initial, const TerminationType& termination,
        Semantics semantics) const;
    //! \brief Full UPPER branch tree with index, event and parent provenance.
    HybridDiscreteTimeUpperExploration explore(
        const EnclosureType& initial, const TerminationType& termination,
        Semantics semantics) const;
  private:
    std::shared_ptr<const SystemType> _system;
    std::shared_ptr<ConfigurationType> _configuration;
};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_DISCRETE_TIME_EVOLVER_HPP
