/***************************************************************************
 *            continuous_time_evolver.hpp
 *
 *  Copyright  2007-20  Alberto Casagrande, Pieter Collins
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

/** @file evolution/continuous_time_evolver.hpp
 *  @brief Validated evolvers for vector fields.
 */

#ifndef ARIADNE_CONTINUOUS_TIME_EVOLVER_HPP
#define ARIADNE_CONTINUOUS_TIME_EVOLVER_HPP

#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "evolution/configuration_interface.hpp"
#include "evolution/enclosure.hpp"
#include "evolution/evolver_interface.hpp"
#include "evolution/orbit.hpp"
#include "function/function_interface.hpp"
#include "integration/integrator_interface.hpp"
#include "logging/logging.hpp"
#include "system/odae.hpp"
#include "threading/workload.hpp"
#include "utility/tuple.hpp"

namespace Ariadne {

/** @brief OrdinaryDifferentialAlgebraicEquation: validated evolvers for vector fields. */
class OrdinaryDifferentialAlgebraicEquation;
template<class ES>
class Orbit;

using Mutex = std::mutex;
template<class T>
using LockGuard = std::lock_guard<T>;

/** @brief ContinuousTimeEvolverConfiguration: validated evolvers for vector fields. */
class ContinuousTimeEvolverConfiguration;
/** @brief PreconditionedTaylorSeriesState: validated evolvers for vector fields. */
class PreconditionedTaylorSeriesState;

/** @brief RealExpressionBoundedConstraintSet: validated evolvers for vector fields. */
class RealExpressionBoundedConstraintSet;

//! \brief A class for computing the evolution of a vector_field system.
//!
//! The actual evolution steps are performed by the Integrator class.
class ContinuousTimeEvolver
    : public EvolverInterface<OrdinaryDifferentialAlgebraicEquation, LabelledEnclosure, typename OrdinaryDifferentialAlgebraicEquation::TimeType> {
  public:
    typedef EvolverInterface<OrdinaryDifferentialAlgebraicEquation, LabelledEnclosure, typename OrdinaryDifferentialAlgebraicEquation::TimeType> Interface;
    typedef ContinuousTimeEvolverConfiguration ConfigurationType;
    typedef OrdinaryDifferentialAlgebraicEquation SystemType;
    typedef IntegratorInterface IntegratorType;
    typedef typename OrdinaryDifferentialAlgebraicEquation::TimeType TimeType;
    typedef Dyadic TimeStepType;
    typedef TimeType TerminationType;
    typedef LabelledEnclosure EnclosureType;
    typedef ListSet<EnclosureType> EnclosureListType;
    /** @brief TimedEnclosureType: validated evolvers for vector fields. */
    struct TimedEnclosureType {
        TimeStepType first;
        EnclosureType second;
        SharedPointer<PreconditionedTaylorSeriesState> preconditioned_state;

        TimedEnclosureType(const TimeStepType& time, const EnclosureType& enclosure)
            : first(time), second(enclosure), preconditioned_state() {}

        TimedEnclosureType(
            const TimeStepType& time,
            const EnclosureType& enclosure,
            SharedPointer<PreconditionedTaylorSeriesState> state)
            : first(time), second(enclosure),
              preconditioned_state(std::move(state)) {}
    };
    typedef Orbit<EnclosureType> OrbitType;
    typedef ValidatedFunctionPatchFactory FunctionFactoryType;

  private:
    //! \brief Synchronised wrapping of orbit to allow concurrent adjoining
    struct SynchronisedOrbit : public OrbitType {
        SynchronisedOrbit(const EnclosureType& initial) : OrbitType(initial) {}
        Void adjoin_reach(const EnclosureType& set) override {
            LockGuard<Mutex> lock(_mux);
            OrbitType::adjoin_reach(set);
        }
        Void adjoin_intermediate(const EnclosureType& set) override {
            LockGuard<Mutex> lock(_mux);
            OrbitType::adjoin_intermediate(set);
        }
        Void adjoin_final(const EnclosureType& set) override {
            LockGuard<Mutex> lock(_mux);
            OrbitType::adjoin_final(set);
        }
        SizeType reach_size() {
            LockGuard<Mutex> lock(_mux);
            return OrbitType::reach().size();
        }

      private:
        Mutex _mux;
    };
    typedef DynamicWorkload<TimedEnclosureType, const TimeType&, Semantics, SharedPointer<SynchronisedOrbit>> WorkloadType;

  public:
    //! \brief Construct from an ODAE without trajectory-constant parameters.
    //! \details Autonomous and explicit-time ODAEs share the same evolution loop.
    //! \throws std::invalid_argument if constant parameters are declared, if
    //! explicit-time auxiliary outputs require an unavailable time-aware
    //! enclosure mapping, or if the selected preconditioned integrator lacks
    //! explicit-time support.
    ContinuousTimeEvolver(const SystemType&, const IntegratorInterface&);

    //! \brief Make a dynamically-allocated copy.
    ContinuousTimeEvolver* clone() const { return new ContinuousTimeEvolver(*this); }

    //! \brief Get the internal system.
    virtual const SystemType& system() const { return *_system; }

    //! \brief Make an enclosure from a computed box set.
    EnclosureType enclosure(const ExactBoxType&) const;

    //!@{
    //! \name Configuration for the class.
    //! \brief A reference to the configuration controlling the evolution.
    ConfigurationType& configuration() { return *this->_configuration; }
    const ConfigurationType& configuration() const { return *this->_configuration; }

    //! \brief The class which constructs functions for the enclosures.
    const FunctionFactoryType& function_factory() const;

    //!@}

    //!@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics.
    OrbitType orbit(const EnclosureType& initial_set, const TimeType& time, Semantics semantics) const;
    /** @brief Perform the specified evolution operation. */
    OrbitType orbit(const RealVariablesBox& initial_set, const TimeType& time, Semantics semantics) const;
    /** @brief Perform the specified evolution operation. */
    OrbitType orbit(const RealExpressionBoundedConstraintSet& initial_set, const TimeType& time, Semantics semantics) const;

    //!@}

  protected:
    Void _process_timed_enclosure(WorkloadType::Access& workload, const TimedEnclosureType& current_timed_set,
                                  const TimeType& maximum_time, Semantics semantics, SharedPointer<SynchronisedOrbit> result) const;

    Void _process_timed_enclosure_step(WorkloadType::Access& workload, const TimedEnclosureType& current_timed_set,
                                       const TimeType& maximum_time, Semantics semantics, SharedPointer<SynchronisedOrbit> result) const;

    /** @brief Perform the specified evolution operation. */
    Void _append_initial_set(WorkloadType& workload, const TimeStepType& initial_time, const EnclosureType& current_set) const;

  private:
    SharedPointer<SystemType> _system;
    SharedPointer<IntegratorType> _integrator;
    SharedPointer<ConfigurationType> _configuration;
};

//! \brief Configuration for a ContinuousTimeEvolver, essentially for controlling the accuracy of continuous evolution methods.
class ContinuousTimeEvolverConfiguration : public ConfigurationInterface {
  public:
    typedef ExactDouble RealType;
    typedef ApproximateDouble ApproximateRealType;

    //! \brief Default constructor gives reasonable values.
    ContinuousTimeEvolverConfiguration();

    ~ContinuousTimeEvolverConfiguration() override = default;

  private:
    //! \brief The maximum allowable step size for integration.
    //! Decreasing this value increases the accuracy of the computation.
    RealType _maximum_step_size;

    //! \brief The maximum allowable radius of a basic set during integration.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    RealType _maximum_enclosure_radius;

    //! \brief The maximum allowable approximation error in the parameter-to-space mapping of an enclosure set.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    RealType _maximum_spacial_error;

    //! \brief Enable reconditioning of basic sets.
    Bool _enable_reconditioning;

    //! \brief Enable subdivisions of basic sets along evolution for upper semantics.
    //! \details Subdivisions are always allowed on the initial set.
    Bool _enable_subdivisions;

  public:
    const RealType& maximum_step_size() const { return _maximum_step_size; }
    Void set_maximum_step_size(const ApproximateRealType value) { _maximum_step_size = cast_exact(value); }

    const RealType& maximum_enclosure_radius() const { return _maximum_enclosure_radius; }
    Void set_maximum_enclosure_radius(const ApproximateRealType value) { _maximum_enclosure_radius = cast_exact(value); }

    const RealType& maximum_spacial_error() const { return _maximum_spacial_error; }
    Void set_maximum_spacial_error(const ApproximateRealType value) { _maximum_spacial_error = cast_exact(value); }

    const Bool& enable_reconditioning() const { return _enable_reconditioning; }
    Void set_enable_reconditioning(const Bool value) { _enable_reconditioning = value; }

    const Bool& enable_subdivisions() const { return _enable_subdivisions; }
    Void set_enable_subdivisions(const Bool value) { _enable_subdivisions = value; }

  public:
    /** @brief Perform the specified evolution operation. */
    OutputStream& _write(OutputStream& os) const override;
};

} // namespace Ariadne

#endif // ARIADNE_CONTINUOUS_TIME_EVOLVER_HPP
