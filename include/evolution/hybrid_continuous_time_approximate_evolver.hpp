/***************************************************************************
 *            hybrid_continuous_time_approximate_evolver.hpp
 *
 *  Copyright  2009-20  Luca Geretti, Mirko Albanese
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

/** @file evolution/hybrid_continuous_time_approximate_evolver.hpp
 *  @brief Simulation algorithms for hybrid dynamical systems.
 */

#ifndef ARIADNE_HYBRID_CONTINUOUS_TIME_APPROXIMATE_EVOLVER_HPP
#define ARIADNE_HYBRID_CONTINUOUS_TIME_APPROXIMATE_EVOLVER_HPP

#include "evolution/configuration_interface.hpp"
#include "evolution/evolver.decl.hpp"
#include "evolution/hybrid_continuous_time_evolver_interface.hpp"
#include "hybrid/hybrid_paving.hpp"
#include "hybrid/hybrid_set.hpp"
#include "logging/logging.hpp"
#include "numeric/numeric.hpp"
#include "symbolic/valuation.hpp"

namespace Ariadne {

/** @brief DiscretizationHybridType represents simulation algorithms for hybrid dynamical systems. */
enum class DiscretizationHybridType {
    HMince,
    HRecombine
};

/** @brief Apply the requested operation. */
inline OutputStream& operator<<(OutputStream& os, const DiscretizationHybridType& dtype) {
    if (dtype == DiscretizationHybridType::HMince) {
        os << "Mince";
    } else {
        os << "Recombine";
    }

    return os;
}

/** @brief HybridTerminationCriterion represents simulation algorithms for hybrid dynamical systems. */
class HybridTerminationCriterion;
/** @brief HybridAutomatonInterface represents simulation algorithms for hybrid dynamical systems. */
class HybridAutomatonInterface;
/** @brief DiscreteEvent represents simulation algorithms for hybrid dynamical systems. */
class DiscreteEvent;
/** @brief DiscreteLocation represents simulation algorithms for hybrid dynamical systems. */
class DiscreteLocation;

/** @brief HybridContinuousTimeEvolverConfiguration<ApproximateTag> represents simulation algorithms for hybrid dynamical systems. */
template<>
class HybridContinuousTimeEvolverConfiguration<ApproximateTag>;

template<class T>
class Orbit;

/*! \brief A class for computing the evolution of a hybrid system.
 */
template<>
class HybridContinuousTimeEvolver<ApproximateTag> : public HybridContinuousTimeEvolverInterface<ApproximateTag> {
  public:
    using Paradigm = ApproximateTag;
    using Interface = HybridContinuousTimeEvolverInterface<ApproximateTag>;
    typedef HybridPoint<FloatDPApproximation> HybridApproximatePointType;
    typedef Vector<HybridApproximatePointType> HybridApproximateListPointType;
    typedef Point<FloatDPApproximation> ApproximatePointType;
    typedef Vector<ApproximatePointType> ApproximatePointListType;
    typedef HybridContinuousTimeEvolverConfiguration<ApproximateTag> ConfigurationType;
    typedef HybridAutomatonInterface SystemType;
    typedef Orbit<HybridApproximatePointType> OrbitType;
    typedef HybridTerminationCriterion TerminationType;
    typedef Orbit<HybridApproximateListPointType> OrbitListType;

  private:
    SharedPointer<SystemType> _sys_ptr;
    SharedPointer<ConfigurationType> _configuration;
    RealValuation _constant_parameters;

  public:
    //! \brief Default constructor.
    HybridContinuousTimeEvolver(const SystemType& system);
    //! \brief Select named parameter values retained through all mode changes.
    HybridContinuousTimeEvolver(const SystemType& system, const RealValuation& constant_parameters);
    const RealValuation& constant_parameters() const { return _constant_parameters; }
    HybridContinuousTimeEvolver* clone() const override { return new HybridContinuousTimeEvolver(*this); }
    const SystemType& system() const override { return *_sys_ptr; }

    //!@{
    //! \name Configuration for the class.
    //! \brief A reference to the configuration controlling the evolution.
    ConfigurationType& configuration() { return *this->_configuration; }
    const ConfigurationType& configuration() const { return *this->_configuration; }

    //!@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics.
    Orbit<HybridApproximatePointType> orbit(const HybridApproximatePointType& initial_point, const TerminationType& termination) const override;
    /** @brief Compute or update evolution information. */
    Orbit<HybridApproximatePointType> orbit(const HybridRealPoint& initial_point, const TerminationType& termination) const;
    //! \brief Currently simulates from the midpoint of the set in the first location only.
    // Orbit<HybridApproximatePointType> orbit(const HybridBoundedConstraintSet& initial_set, const TerminationType& termination) const;

    OrbitListType orbit(const HybridApproximateListPointType& initial_list, const TerminationType& termination) const; // TO DO
    OrbitListType orbit(const HybridBoundedConstraintSet& initial_set, const TerminationType& termination) const;      // TO DO
    /** @brief Compute or update evolution information. */
    OrbitListType orbit(const HybridUpperBox& initial_box, DiscreteLocation loc, HybridSpace spc, const TerminationType& termination) const;

  private:
    /** @brief Perform the specified evolution operation. */
    Map<DiscreteEvent, EffectiveScalarMultivariateFunction> _guard_functions(const DiscreteLocation& location) const;
    /** @brief Perform the specified evolution operation. */
    Bool _satisfies_invariants(const DiscreteLocation& location, const Point<FloatDPApproximation>& point, const FloatDPApproximation& time) const;
};

/** @brief HybridContinuousTimeEvolverConfiguration<ApproximateTag> represents simulation algorithms for hybrid dynamical systems. */
template<>
class HybridContinuousTimeEvolverConfiguration<ApproximateTag> : public ConfigurationInterface {
  public:
    using Paradigm = ApproximateTag;
    //! \brief Default constructor gives reasonable values.
    HybridContinuousTimeEvolverConfiguration();

    virtual ~HybridContinuousTimeEvolverConfiguration() = default;

  private:
    //! \brief The step size for integration.
    //! Decreasing this value increases the accuracy of the computation.
    FloatDPApproximation _step_size;
    Nat _fineness;
    DiscretizationHybridType _discretization_type;

  public:
    const FloatDPApproximation& step_size() const { return _step_size; }
    Void set_step_size(ApproximateDouble value) {
        ARIADNE_PRECONDITION(decide(value > 0));
        _step_size = FloatDPApproximation(value, double_precision);
    }
    const Nat& fineness() const { return _fineness; }
    Void set_fineness(Nat value) { _fineness = value; }
    Void set_d_type(DiscretizationHybridType type) { _discretization_type = type; }
    const DiscretizationHybridType& d_type() const { return _discretization_type; }

    /** @brief Perform the specified evolution operation. */
    virtual OutputStream& _write(OutputStream& os) const;
};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_CONTINUOUS_TIME_APPROXIMATE_EVOLVER_HPP
