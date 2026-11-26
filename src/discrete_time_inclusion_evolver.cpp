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
    if (semantics != Semantics::UPPER && semantics != Semantics::LOWER)
        throw std::invalid_argument("Unknown ordinary discrete inclusion evolution semantics");
    if (stop < Integer(0))
        throw std::invalid_argument("Ordinary discrete inclusion horizon must be nonnegative");
    const auto& model = this->system();
    const Bool has_auxiliary = model.auxiliary_space().dimension() != 0u;
    if (initial.state_space() != model.state_space() ||
        (initial.auxiliary_space().dimension() != 0u &&
         initial.auxiliary_space() != model.auxiliary_space()))
        throw std::invalid_argument("Ordinary discrete inclusion initial enclosure spaces mismatch");

    EnclosureType seed = initial;
    if (has_auxiliary) {
        if (seed.auxiliary_space().dimension() == 0u) {
            // At the initial tick, a[0] depends on an as-yet unchosen u[0].
            Enclosure extended = seed.euclidean_set();
            extended.adjoin_parameters(model.inputs(), EnclosureVariableKind::INPUT);
            seed = EnclosureType(extended, model.state_space());
        } else {
            // A supplied auxiliary-bearing enclosure (e.g. a resumed orbit)
            // must already bind the input choices for its current index.
            if (seed.auxiliary_arguments() != EnclosureFunctionArguments::STATE_INPUTS ||
                seed.auxiliary_parameter_count() != model.number_of_inputs())
                throw std::invalid_argument("Discrete inclusion resume requires current-input auxiliary binding");
            SizeType remaining = model.number_of_inputs();
            const auto domain = seed.parameter_domain();
            const auto& kinds = seed.variable_kinds();
            for (SizeType i = kinds.size(); i != 0u && remaining != 0u; --i) {
                if (kinds[i - 1u] == EnclosureVariableKind::INPUT) {
                    --remaining;
                    if (!definitely(subset(domain[i - 1u], model.inputs()[remaining])))
                        throw std::invalid_argument("Resumed inclusion input is outside its declared box");
                }
            }
            if (remaining != 0u)
                throw std::invalid_argument("Resumed inclusion is missing current input coordinates");
        }
        seed.set_auxiliary(model.auxiliary_space(), model.auxiliary_function(),
                           EnclosureFunctionArguments::STATE_INPUTS, model.number_of_inputs());
    } else if (seed.auxiliary_space().dimension() != 0u) {
        throw std::invalid_argument("Ordinary inclusion without outputs cannot accept auxiliary labels");
    }

    OrbitType result(seed);
    ListSet<EnclosureType> final_sets, reach_sets, intermediate_sets;
    List<Pair<TimeType, EnclosureType>> working;
    working.push_back(make_pair(TimeType(0), seed));

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
        if (semantics == Semantics::UPPER &&
            _configuration->enable_subdivisions() &&
            decide(radius > _configuration->maximum_enclosure_radius())) {
            const auto pieces = source.split();
            working.push_back(make_pair(tick, pieces.first));
            working.push_back(make_pair(tick, pieces.second));
            continue;
        }

        // For modes without auxiliaries, the current tick's full input box
        // is introduced now. With auxiliaries, u[n] was already installed
        // at the source tick, where a[n]=g(x[n],u[n]) refers to it.
        Enclosure evolved = source.euclidean_set();
        if (!has_auxiliary)
            evolved.adjoin_parameters(model.inputs(), EnclosureVariableKind::INPUT);
        const auto inputs = evolved.input_parameter_function(model.number_of_inputs());
        const auto mapped = compose(model.update_function(),
                                    join(evolved.state_function(), inputs));
        evolved.apply_parameterised_fixed_evolve_step(mapped, StepSizeType(1, 0u));
        if (has_auxiliary) {
            // u[n+1] is independent of u[n] and all earlier tick choices.
            // STATE_INPUTS resolves to these newest coordinates; retaining
            // the old INPUT axes keeps state/input/history correlations.
            evolved.adjoin_parameters(model.inputs(), EnclosureVariableKind::INPUT);
        }
        const EnclosureType successor(evolved, model.state_space(), model.auxiliary_space());
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
