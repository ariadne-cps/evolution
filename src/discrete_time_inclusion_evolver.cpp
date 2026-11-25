/***************************************************************************
 *            discrete_time_inclusion_evolver.cpp
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

/** @file src/discrete_time_inclusion_evolver.cpp
 *  @brief Conservative whole-input images at successive integer ticks.
 */

#include "evolution/discrete_time_inclusion_evolver.hpp"

#include <stdexcept>

#include "evolution/enclosure.hpp"
#include "function/function.hpp"
#include "function/functional.hpp"
#include "function/taylor_function.decl.hpp"
#include "utility/tuple.hpp"

namespace Ariadne {

DiscreteTimeInclusionEvolver::DiscreteTimeInclusionEvolver(const SystemType& system)
    : _system(system.clone()), _configuration(std::make_shared<ConfigurationType>()) {}

DiscreteTimeInclusionEvolver::DiscreteTimeInclusionEvolver(const DiscreteTimeInclusionEvolver& other)
    : _system(other._system->clone()),
      _configuration(std::make_shared<ConfigurationType>(*other._configuration)) {}

DiscreteTimeInclusionEvolver::FunctionFactoryType
DiscreteTimeInclusionEvolver::function_factory() const {
    return ValidatedFunctionPatchFactory(make_taylor_function_patch_factory());
}

DiscreteTimeInclusionEvolver::EnclosureType
DiscreteTimeInclusionEvolver::enclosure(const ExactBoxType& box) const {
    if (box.dimension() != this->system().dimension())
        throw std::invalid_argument("Ordinary discrete inclusion initial box dimension mismatch");
    return EnclosureType(box, this->system().state_space(),
                         EnclosureConfiguration(this->function_factory()));
}

DiscreteTimeInclusionEvolver::OrbitType
DiscreteTimeInclusionEvolver::orbit(const RealExpressionBoundedConstraintSet& initial,
                                    const TerminationType& stop, Semantics semantics) const {
    return orbit(EnclosureType(initial.euclidean_set(this->system().state_space()),
                               this->system().state_space(),
                               EnclosureConfiguration(this->function_factory())),
                 stop, semantics);
}

DiscreteTimeInclusionEvolver::OrbitType
DiscreteTimeInclusionEvolver::orbit(const EnclosureType& initial,
                                    const TerminationType& stop, Semantics semantics) const {
    if (semantics != Semantics::UPPER)
        throw std::invalid_argument("Validated ordinary discrete inclusions support only UPPER");
    if (stop < Integer(0))
        throw std::invalid_argument("Ordinary discrete inclusion horizon must be nonnegative");
    if (this->system().auxiliary_space().dimension() != 0u)
        throw std::invalid_argument("Validated ordinary discrete inclusion auxiliaries are not yet supported");
    if (initial.state_space() != this->system().state_space() ||
        initial.state_auxiliary_space() != this->system().state_space())
        throw std::invalid_argument("Ordinary discrete inclusion initial enclosure state space mismatch");

    OrbitType result(initial);
    ListSet<EnclosureType> final_sets, reach_sets, intermediate_sets;
    List<Pair<TimeType, EnclosureType>> working;
    working.push_back(make_pair(TimeType(0), initial));

    while (!working.empty()) {
        const auto current = working.back();
        working.pop_back();
        const TimeType tick = current.first;
        const EnclosureType& source = current.second;
        if (tick >= stop) {
            final_sets.adjoin(source);
            continue;
        }

        const FloatDPUpperBound radius = source.euclidean_set().bounding_box().radius();
        if (_configuration->enable_subdivisions() &&
            decide(radius > _configuration->maximum_enclosure_radius())) {
            const auto pieces = source.split();
            working.push_back(make_pair(tick, pieces.first));
            working.push_back(make_pair(tick, pieces.second));
            continue;
        }

        // Reuse the validated UPPER image construction of the hybrid
        // discrete inclusion progression: each tick receives independent
        // full-box input coordinates; no sample or unproved LOWER witness.
        Enclosure evolved = source.euclidean_set();
        const SizeType input_offset = evolved.number_of_parameters();
        evolved.adjoin_parameters(this->system().inputs(), EnclosureVariableKind::INPUT);
        const ExactBoxType domain = evolved.parameter_domain();
        const auto& factory = evolved.configuration().function_factory();
        auto inputs = factory.create_zeros(this->system().number_of_inputs(), domain);
        for (SizeType j = 0u; j != this->system().number_of_inputs(); ++j)
            inputs[j] = factory.create_coordinate(domain, input_offset + j);

        const auto mapped = compose(this->system().update_function(),
                                    join(evolved.state_function(), inputs));
        evolved.apply_parameterised_fixed_evolve_step(mapped, StepSizeType(1, 0u));
        const EnclosureType successor(evolved, this->system().state_space());
        reach_sets.adjoin(successor);
        intermediate_sets.adjoin(successor);
        working.push_back(make_pair(tick + TimeType(1), successor));
    }

    result.adjoin_reach(reach_sets);
    result.adjoin_intermediate(intermediate_sets);
    result.adjoin_final(final_sets);
    return result;
}

} // namespace Ariadne
