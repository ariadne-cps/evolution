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

#include "evolution/configuration_interface.hpp"
#include "evolution/orbit.hpp"
#include "geometry/point.hpp"
#include "logging/logging.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/floatdp.hpp"
#include "system/vector_field.hpp"
#include "threading/using.hpp"
#include "threading/workload.hpp"

namespace Ariadne {

/** @brief DiscretisationType: numerical simulation of continuous vector-field systems. */
enum class DiscretisationType {
    None,
    Mince,
    Recombine
};

/** @brief Perform the specified evolution operation. */
inline OutputStream& operator<<(OutputStream& os, const DiscretisationType& dtype) {
    if (dtype == DiscretisationType::Mince) {
        os << "Mince";
    } else if (dtype == DiscretisationType::Recombine) {
        os << "Recombine";
    } else {
        os << "None";
    }
    return os;
}

/** @brief VectorField: numerical simulation of continuous vector-field systems. */
class VectorField;
/** @brief VectorFieldSimulatorConfiguration: numerical simulation of continuous vector-field systems. */
class VectorFieldSimulatorConfiguration;

/** @brief RealExpressionBoundedConstraintSet: numerical simulation of continuous vector-field systems. */
class RealExpressionBoundedConstraintSet;

template<class T>
class Orbit;

//! \brief The subdivision depth associated with each variable of a subspace grid.
using subspace = Map<RealVariable, Nat>;

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
    Nat _num_subdivisions;
    DiscretisationType _discretisation_type;
    subspace _subspace;

  public:
    Void set_step_size(double h) {
        ARIADNE_PRECONDITION(h > 0);
        _step_size = h;
    }
    const FloatDPApproximation& step_size() const { return _step_size; }

    Void set_discretisation_type(DiscretisationType discretisation_type) { _discretisation_type = discretisation_type; }
    const DiscretisationType& discretisation_type() const { return _discretisation_type; }

    Void set_num_subdivisions(Nat num_subdivisions) { _num_subdivisions = num_subdivisions; }
    const Nat& num_subdivisions() const { return _num_subdivisions; }

    /** @brief Perform the specified evolution operation. */
    Void insert_subspace(RealVariable v, Nat n) {
        /** @brief Construct the object. */
        ARIADNE_PRECONDITION(n > 0u);
        /** @brief Construct the object. */
        ARIADNE_PRECONDITION(!_subspace.has_key(v));
        _subspace.insert(v, n);
    }

    /** @brief Perform the specified evolution operation. */
    Void remove_keys() {
        _subspace.remove_keys(_subspace.keys());
    }

    /** @brief Perform the specified evolution operation. */
    Void set_subspace(subspace _subspace_) {
        for (const auto& entry: _subspace_) {
            ARIADNE_PRECONDITION(entry.second > 0u);
        }
        _subspace = _subspace_;
    }

    Bool is_using_subspace() const { return !_subspace.empty(); }

    const subspace& get_subspace() const { return _subspace; }

    /** @brief Perform the specified evolution operation. */
    virtual OutputStream& _write(OutputStream& os) const;
};

} // namespace Ariadne

#endif // ARIADNE_VECTOR_FIELD_SIMULATOR_HPP
