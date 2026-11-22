/***************************************************************************
 *            evolver_interface.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/** @file evolution/evolver_interface.hpp
 *  @brief Interfaces and semantics of time evolution.
 */

#ifndef ARIADNE_EVOLVER_INTERFACE_HPP
#define ARIADNE_EVOLVER_INTERFACE_HPP

#include <cstdint>

#include "evolution/evolver.decl.hpp"

namespace Ariadne {

template<class ES>
class ListSet;
template<class ES>
class Orbit;

//! \ingroup DynamicsModule
//! \brief The semantics used to determine the trajectories of the system.
//! \see EvolverInterface
enum class Semantics : std::uint8_t {
    LOWER, //!< Under-approximation with trajectories terminating at spacial discontinuities.
    UPPER  //!< Over-approximations with all possibilities included as spacial discontinuities.
};

/** @brief Paradigm-specific evolution contract. Unsupported paradigms are incomplete. */
template<class P, class SYS, class ST, class TRM>
class EvolverInterface;

/** @brief Validated evolution returns enclosures under explicit set semantics. */
template<class SYS, class ES, class TRM>
class EvolverInterface<ValidatedTag, SYS, ES, TRM> {
  public:
    using Paradigm = ValidatedTag;
    using SystemType = SYS;
    using TimeType = typename SYS::TimeType;
    using TerminationType = TRM;
    using InitialType = ES;
    using EnclosureType = ES;
    using EnclosureListType = ListSet<ES>;
    using OrbitType = Orbit<ES>;
    using ResultType = OrbitType;
    virtual ~EvolverInterface() = default;
    virtual EvolverInterface* clone() const = 0;
    virtual const SystemType& system() const = 0;
    virtual OrbitType orbit(const InitialType&, const TerminationType&, Semantics) const = 0;
};

/** @brief Approximate evolution returns sampled trajectories without set semantics. */
template<class SYS, class PT, class TRM>
class EvolverInterface<ApproximateTag, SYS, PT, TRM> {
  public:
    using Paradigm = ApproximateTag;
    using SystemType = SYS;
    using TimeType = typename SYS::TimeType;
    using TerminationType = TRM;
    using InitialType = PT;
    using OrbitType = Orbit<PT>;
    using ResultType = OrbitType;
    virtual ~EvolverInterface() = default;
    virtual EvolverInterface* clone() const = 0;
    virtual const SystemType& system() const = 0;
    virtual OrbitType orbit(const InitialType&, const TerminationType&) const = 0;
};

/** @brief Factory retaining the paradigm-specific input and result contract. */
template<class P, class SYS, class ST, class TRM>
class EvolverFactoryInterface {
  public:
    virtual ~EvolverFactoryInterface() = default;
    virtual EvolverFactoryInterface* clone() const = 0;
    virtual EvolverInterface<P, SYS, ST, TRM>* create(const SYS&) const = 0;
};

} // namespace Ariadne
#endif // ARIADNE_EVOLVER_INTERFACE_HPP
