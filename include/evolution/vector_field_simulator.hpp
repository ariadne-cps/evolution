/***************************************************************************
 *            vector_field_simulator.hpp
 *
 *  Copyright  2009-21  Luca Geretti, Mirko Albanese
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

/** @file evolution/vector_field_simulator.hpp
 *  @brief Numerical simulation of continuous vector-field systems.
 */

#ifndef ARIADNE_VECTOR_FIELD_SIMULATOR_HPP
#define ARIADNE_VECTOR_FIELD_SIMULATOR_HPP

#include <utility>

#include "evolution/configuration_interface.hpp"
#include "evolution/orbit.hpp"
#include "evolution/sampling_counts.hpp"
#include "geometry/point.hpp"
#include "logging/logging.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/floatdp.hpp"
#include "system/vector_field.hpp"
#include "threading/using.hpp"
#include "threading/workload.hpp"

namespace Ariadne {

/** @brief VectorField: numerical simulation of continuous vector-field systems. */
class VectorField;
/** @brief VectorFieldSimulatorConfiguration: numerical simulation of continuous vector-field systems. */
class VectorFieldSimulatorConfiguration;

/** @brief RealExpressionBoundedConstraintSet: numerical simulation of continuous vector-field systems. */
class RealExpressionBoundedConstraintSet;

template<class T>
class Orbit;

/*! \brief A class for computing the simulated evolution of a continuous system.
 */
class VectorFieldSimulator {
  public:
    typedef FloatDPApproximation ApproximateTimeType;
    typedef LabelledPoint<FloatDPApproximation> ApproximatePointType;
    typedef List<ApproximatePointType> ApproximateListPointType;
    typedef LabelledPoint<Real> RealPointType;
    typedef List<RealPointType> RealListPointType;
    typedef RealVariablesBox RealBoxType;
    typedef Real TerminationType;
    typedef VectorField SystemType;
    typedef VectorFieldSimulatorConfiguration ConfigurationType;
    typedef Orbit<ApproximateListPointType> OrbitType;

    //! \brief Synchronised wrapping of orbit to allow concurrent adjoining
    struct SynchronisedOrbit : public OrbitType {
        SynchronisedOrbit(const ApproximateListPointType& initial_points) : OrbitType(initial_points) {}
        /** @brief Perform the specified evolution operation. */
        void insert(const FloatDP& t, const ApproximatePointType& pt, const SizeType& curve_number) {
            std::lock_guard<std::mutex> lock(_mux);
            OrbitType::insert(t, pt, curve_number);
        }

      private:
        std::mutex _mux;
    };
    typedef StaticWorkload<Pair<SizeType, ApproximatePointType>, const TerminationType&, SharedPointer<SynchronisedOrbit>> WorkloadType;

  public:
    //! \brief Default constructor.
    VectorFieldSimulator(const SystemType& system);

    //!@{
    //! \name Configuration for the class.
    //! \brief A reference to the configuration controlling the evolution.
    ConfigurationType& configuration() { return *this->_configuration; }
    const ConfigurationType& configuration() const { return *this->_configuration; }

    //!@{
    //! \name Simulation using points.
    //! \brief Compute an approximation to the orbit set.
    OrbitType orbit(const RealBoxType& initial_box, const TerminationType& termination) const;
    /** @brief Perform the specified evolution operation. */
    OrbitType orbit(const RealExpressionBoundedConstraintSet& initial_set, const TerminationType& termination) const;
    /** @brief Perform the specified evolution operation. */
    OrbitType orbit(const ApproximateListPointType& initial_list, const TerminationType& termination) const;
    /** @brief Perform the specified evolution operation. */
    OrbitType orbit(UpperBoxType& initial_box, const TerminationType& termination) const;

  private:
    /** @brief Perform the specified evolution operation. */
    void _simulate_from_point(Pair<SizeType, ApproximatePointType> indexed_initial, const TerminationType& termination, SharedPointer<SynchronisedOrbit> orbit) const;

  private:
    SharedPointer<SystemType> _system;
    SharedPointer<ConfigurationType> _configuration;
};

//! \brief Configuration for a VectorFieldSimulator, essentially to control accuracy of evolution.
class VectorFieldSimulatorConfiguration : public ConfigurationInterface {
  public:
    //! \brief Default constructor gives reasonable values.
    VectorFieldSimulatorConfiguration();

    virtual ~VectorFieldSimulatorConfiguration() = default;

  private:
    //! \brief The fixed integration step size to use.
    FloatDPApproximation _step_size;
    SamplingCounts _sampling_counts;

  public:
    Void set_step_size(double h) {
        ARIADNE_PRECONDITION(h > 0);
        _step_size = h;
    }
    const FloatDPApproximation& step_size() const { return _step_size; }

    /** @brief Configure the number of sampling subintervals per named axis.
     * Axes not present in the map contribute one sample.
     */
    Void set_sampling_counts(SamplingCounts counts) {
        for (const auto& entry: counts)
            ARIADNE_PRECONDITION(entry.second > 0u);
        _sampling_counts = std::move(counts);
    }

    /** @brief Set a positive count for a previously unconfigured axis. */
    Void set_sampling_count(RealVariable variable, Nat count) {
        ARIADNE_PRECONDITION(count > 0u);
        ARIADNE_PRECONDITION(!_sampling_counts.has_key(variable));
        _sampling_counts.insert(variable, count);
    }

    /** @brief Reset to one sample per axis. */
    Void clear_sampling_counts() { _sampling_counts.remove_keys(_sampling_counts.keys()); }

    const SamplingCounts& sampling_counts() const { return _sampling_counts; }

    /** @brief Perform the specified evolution operation. */
    virtual OutputStream& _write(OutputStream& os) const;
};

} // namespace Ariadne

#endif // ARIADNE_VECTOR_FIELD_SIMULATOR_HPP
