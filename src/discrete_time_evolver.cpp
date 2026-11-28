/***************************************************************************
 *            discrete_time_evolver.cpp
 *
 *  Copyright  2008-20  Alberto Casagrande, Pieter Collins
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

/** @file src/discrete_time_evolver.cpp
 *  @brief Implement difference equation evolver operations.
 */

#include "evolution/discrete_time_evolver.hpp"

#include "algebra/vector.hpp"
#include "evolution/enclosure.hpp"
#include "evolution/orbit.hpp"
#include "function/constraint.hpp"
#include "function/function.hpp"
#include "function/functional.hpp"
#include "function/taylor_function.decl.hpp"
#include "logging/logging.hpp"
#include "symbolic/assignment.hpp"
#include "system/difference_algebraic_equation.hpp"
#include "utility/array.hpp"
#include "utility/macros.hpp"
#include "utility/stlio.hpp"
#include "utility/tuple.hpp"

namespace Ariadne {

namespace {

template<class ES>
List<ES> subdivide(const ES& enclosure) {
    List<ES> result;
    Pair<ES, ES> split = enclosure.split();
    result.append(split.first);
    result.append(split.second);
    return result;
}

} // namespace

DiscreteTimeEvolver<ValidatedTag>::DiscreteTimeEvolver(const SystemType& system)
    : _system(system.clone()), _configuration(new ConfigurationType()) {
}

const typename DiscreteTimeEvolver<ValidatedTag>::FunctionFactoryType DiscreteTimeEvolver<ValidatedTag>::function_factory() const {
    return ValidatedFunctionPatchFactory(make_taylor_function_patch_factory());
}

typename DiscreteTimeEvolver<ValidatedTag>::EnclosureType DiscreteTimeEvolver<ValidatedTag>::enclosure(const ExactBoxType& box) const {
    EnclosureType result(box, this->system().state_space(), EnclosureConfiguration(this->function_factory()));
    result.set_auxiliary(this->system().auxiliary_space(), this->system().auxiliary_function());
    return result;
}

Orbit<DiscreteTimeEvolver<ValidatedTag>::EnclosureType>
    DiscreteTimeEvolver<ValidatedTag>::
        orbit(const RealExpressionBoundedConstraintSet& initial_set,
              const TerminationType& termination,
              Semantics semantics) const {
    return orbit(EnclosureType(initial_set.euclidean_set(this->system().state_space()), this->system().state_space(), EnclosureConfiguration(this->function_factory())), termination, semantics);
}

Orbit<DiscreteTimeEvolver<ValidatedTag>::EnclosureType>
    DiscreteTimeEvolver<ValidatedTag>::
        orbit(const EnclosureType& initial_set,
              const TerminationType& termination,
              Semantics semantics) const {
    ARIADNE_PRECONDITION(this->system().state_space() == initial_set.state_space())
    ARIADNE_PRECONDITION(initial_set.auxiliary_space().size() == 0u ||
                         this->system().auxiliary_space() == initial_set.auxiliary_space())
    // Outputs are the declared correlated map of the current state, including
    // at zero horizon and when resuming a caller-provided state enclosure.
    EnclosureType prepared = initial_set;
    prepared.set_auxiliary(this->system().auxiliary_space(), this->system().auxiliary_function());
    Orbit<EnclosureType> orbit(prepared);
    EnclosureListType final;
    EnclosureListType reachable;
    EnclosureListType intermediate;
    this->_evolution(final, reachable, intermediate,
                     prepared, termination, semantics);
    orbit.adjoin_intermediate(intermediate);
    orbit.adjoin_reach(reachable);
    orbit.adjoin_final(final);
    return orbit;
}

Void DiscreteTimeEvolver<ValidatedTag>::
    _evolution(EnclosureListType& final_sets,
               EnclosureListType& reach_sets,
               EnclosureListType& intermediate_sets,
               const EnclosureType& initial_set,
               const TerminationType& maximum_time,
               Semantics semantics) const {
    LOGGING_SCOPE_CREATE;

    List<TimedEnclosureType> working_sets;

    {
        // Set up initial timed set models
        LOGGING_PRINTLN_AT(1, "initial_set = " << initial_set);
        EnclosureType initial_enclosure(initial_set);
        LOGGING_PRINTLN_AT(1, "initial_enclosure = " << initial_enclosure);
        TimeType initial_time = 0;
        LOGGING_PRINTLN_AT(1, "initial_time = " << initial_time);
        working_sets.push_back(make_pair(initial_time, initial_enclosure));
    }

    while (!working_sets.empty()) {
        TimedEnclosureType current_set = working_sets.back();
        working_sets.pop_back();
        EnclosureType initial_enclosure = current_set.second;
        TimeType initial_time = current_set.first;
        FloatDPUpperBound initial_set_radius = initial_enclosure.euclidean_set().bounding_box().radius();
        if (initial_time >= maximum_time) {
            final_sets.adjoin(EnclosureType(initial_enclosure));
        } else if (semantics == Semantics::UPPER && this->_configuration->enable_subdivisions() && decide(initial_set_radius > this->_configuration->maximum_enclosure_radius())) {
            // Subdivide
            List<EnclosureType> subdivisions = subdivide(initial_enclosure);
            for (SizeType i = 0; i != subdivisions.size(); ++i) {
                const EnclosureType& subdivided_enclosure = subdivisions[i];
                working_sets.push_back(make_pair(initial_time, subdivided_enclosure));
            }
        } else if (semantics == Semantics::LOWER && this->_configuration->enable_premature_termination() && decide(initial_set_radius > this->_configuration->maximum_enclosure_radius())) {
            ARIADNE_WARN("Terminating lower evolution at time " << initial_time << " and set " << initial_enclosure << " due to maximum radius being exceeded.");
        } else {
            // Compute evolution
            this->_evolution_step(working_sets, final_sets, reach_sets, intermediate_sets, current_set, maximum_time, semantics);
        }
    }
}

Void DiscreteTimeEvolver<ValidatedTag>::
    _evolution_step(List<TimedEnclosureType>& working_sets,
                    EnclosureListType&,
                    EnclosureListType& reach_sets,
                    EnclosureListType& intermediate_sets,
                    const TimedEnclosureType& current_set,
                    const TimeType&,
                    Semantics) const {
    LOGGING_SCOPE_CREATE;

    EnclosureType initial_enclosure;
    TimeType initial_time;

    LOGGING_PRINTLN_AT(1, "working_set = " << current_set);
    make_lpair(initial_time, initial_enclosure) = current_set;

    LOGGING_PRINTLN_AT(1, "initial_time = " << initial_time);
    LOGGING_PRINTLN_AT(1, "initial_enclosure = " << initial_enclosure);

    LOGGING_PRINTLN_AT(1, "box = " << initial_enclosure.bounding_box() << " ");
    LOGGING_PRINTLN_AT(1, "radius = " << initial_enclosure.euclidean_set().bounding_box().radius() << "\n");

    /////////////// Main Evolution ////////////////////////////////

    // Compute the map model
    EnclosureType& final_enclosure = initial_enclosure;
    final_enclosure.apply_discrete_time_map_step(_system->function());
    TimeType final_time = initial_time + 1;
    LOGGING_PRINTLN("final_enclosure = " << final_enclosure);

    reach_sets.adjoin(final_enclosure);

    intermediate_sets.adjoin(final_enclosure);
    working_sets.push_back(make_pair(final_time, final_enclosure));
}

DiscreteTimeEvolverConfiguration<ValidatedTag>::DiscreteTimeEvolverConfiguration() {
    set_maximum_enclosure_radius(100);
    set_enable_subdivisions(false);
    set_enable_premature_termination(false);
}

OutputStream&
    DiscreteTimeEvolverConfiguration<ValidatedTag>::_write(OutputStream& os) const {
    os << "DiscreteTimeEvolverConfiguration<ValidatedTag>"
       << ",\n  maximum_enclosure_radius=" << maximum_enclosure_radius()
       << ",\n  enable_subdivisions=" << enable_subdivisions()
       << ",\n  enable_premature_termination=" << enable_premature_termination()
       << "\n)\n";
    return os;
}

} // namespace Ariadne

