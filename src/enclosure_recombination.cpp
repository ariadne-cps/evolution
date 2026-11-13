/***************************************************************************
 *            enclosure_recombination.cpp
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

/** @file enclosure_recombination.cpp
 *  @brief Outer recombination of constrained continuous enclosures.
 */

#include <algorithm>

#include "evolution/enclosure.hpp"
#include "function/formula.hpp"
#include "function/formula.tpl.hpp"
#include "function/function.hpp"
#include "utility/macros.hpp"

#include "constraint-solving/constraint_solver.hpp"

namespace Ariadne {
namespace {

// Bound an image on the feasible parameter set by contracting its graph.
// Keeping each source's own constraints here avoids intersecting alternatives.
UpperBoxType constrained_image_bounds(const Enclosure& source,
                                      const ValidatedVectorMultivariateFunctionPatch& image) {
    const auto parameter_domain = source.reduced_domain();
    UpperBoxType bounds = apply(image, parameter_domain);
    // A source predicate g in C also bounds any image f by C+range(f-g).
    // In particular, this retains equalities shared with the reference even
    // when an oblique constraint cannot contract the parameter box itself.
    for (SizeType k = 0; k != image.result_size(); ++k) {
        for (const auto& constraint: source.constraint_models()) {
            bounds[k] = intersection(bounds[k], UpperIntervalType(constraint.bounds()) + apply(image.get(k) - constraint.function(), parameter_domain));
        }
    }
    if (definitely(bounds.is_empty()))
        return bounds;
    ARIADNE_PRECONDITION_MSG(definitely(bounds.is_bounded()), "Unbounded recombination image: " << bounds);
    if (source.number_of_constraints() == 0u || parameter_domain.size() == 0u)
        return bounds;

    const auto image_domain = cast_exact_box(bounds);
    const auto graph_domain = product(parameter_domain, image_domain);
    const auto& factory = source.configuration().function_factory();
    const auto parameter_count = parameter_domain.size();
    const auto constraint_count = source.number_of_constraints();
    auto graph = factory.create_zeros(constraint_count + image.result_size(), graph_domain);
    ExactBoxType graph_bounds(graph.result_size());
    for (SizeType k = 0; k != constraint_count; ++k) {
        graph.set(k, embed(restriction(source.constraint(k).function(), parameter_domain), image_domain));
        graph_bounds[k] = source.constraint(k).bounds();
    }
    for (SizeType k = 0; k != image.result_size(); ++k) {
        graph.set(constraint_count + k,
                  embed(restriction(image.get(k), parameter_domain), image_domain) - factory.create_coordinate(graph_domain, parameter_count + k));
        graph_bounds[constraint_count + k] = ExactIntervalType(0, 0);
    }
    UpperBoxType feasible_graph(graph_domain);
    ConstraintSolver contractor;
    const PreparedConstraintSolverFunction prepared(graph);
    contractor.reduce(feasible_graph, prepared, graph_bounds);
    if (definitely(feasible_graph.is_empty())) {
        bounds[0] = UpperIntervalType(1, -1);
        return bounds;
    }
    return project(feasible_graph, range(parameter_count, parameter_count + image.result_size()));
}

} // namespace

Enclosure recombine(const List<Enclosure>& enclosures) {
    ARIADNE_PRECONDITION_MSG(!enclosures.empty(), "Recombination requires a source enclosure.");
    ARIADNE_PRECONDITION_MSG(std::all_of(enclosures.begin(), enclosures.end(), [](const Enclosure& source) {
                                 return source._state_function.raw_pointer() != nullptr;
                             }),
                             "Recombination requires initialised enclosure models.");
    const auto state_dimension = enclosures.front().state_dimension();
    List<Enclosure> sources;
    sources.reserve(enclosures.size());
    SizeType parameter_count = 0u;
    for (auto source: enclosures) {
        ARIADNE_PRECONDITION_MSG(source.state_dimension() == state_dimension,
                                 "Recombination requires matching state dimensions.");
        // is_empty contracts the copy's parameter domain and only discards
        // sources whose infeasibility has actually been established.
        if (definitely(source.is_empty()))
            continue;
        ARIADNE_PRECONDITION_MSG(source.reduced_domain().is_bounded(),
                                 "Recombination requires bounded parameter domains.");
        parameter_count = std::max(parameter_count, source.number_of_parameters());
        sources.append(source);
    }
    if (sources.empty())
        return enclosures.front();
    if (sources.size() == 1u)
        return sources.front();

    const auto& reference = sources.front();
    const auto auxiliary_formula = make_formula(reference.auxiliary_mapping());
    for (const auto& source: sources) {
        ARIADNE_PRECONDITION_MSG(identical(auxiliary_formula, make_formula(source.auxiliary_mapping())),
                                 "Recombination requires identical auxiliary mappings.");
    }
    const auto& factory = reference.configuration().function_factory();
    const ExactBoxType common_domain(parameter_count, ExactIntervalType(0, 1));
    const auto normalise = [&](const Enclosure& source) {
        const auto reduced = source.reduced_domain();
        auto parameters = factory.create_zeros(reduced.size(), common_domain);
        for (SizeType k = 0; k != reduced.size(); ++k) {
            const auto lower = reduced[k].lower_bound();
            const auto upper = reduced[k].upper_bound();
            const auto coordinate = factory.create_coordinate(common_domain, k);
            parameters.set(k, FloatDPBounds(lower) * (1 - coordinate) + FloatDPBounds(upper) * coordinate);
        }
        Enclosure transformed(reference.configuration());
        transformed._domain = common_domain;
        transformed._reduced_domain = common_domain;
        // Each exact affine parameter map covers its own reduced source box.
        // Unchecked composition avoids rejecting outward-rounding spill at
        // the boundary of that box; no model is extrapolated mathematically.
        transformed._state_function = unchecked_compose(source.state_function(), parameters);
        transformed._time_function = unchecked_compose(source.time_function(), parameters);
        transformed._dwell_time_function = unchecked_compose(source.dwell_time_function(), parameters);
        transformed._auxiliary_mapping = reference.auxiliary_mapping();
        transformed._variable_kinds = List<EnclosureVariableKind>(parameter_count, EnclosureVariableKind::UNKNOWN);
        for (SizeType k = 0; k != source.number_of_parameters(); ++k)
            transformed._variable_kinds[k] = source.variable_kinds()[k];
        for (const auto& constraint: source.constraint_models()) {
            transformed._constraints.append(ValidatedConstraintModel(
                constraint.lower_bound(), unchecked_compose(constraint.function(), parameters), constraint.upper_bound()));
        }
        return transformed;
    };

    Enclosure result = normalise(sources.front());
    const auto reference_map = join(join(result.state_function(), result.time_function()), result.dwell_time_function());
    const auto component_count = reference_map.result_size();
    const auto constraint_count = result.number_of_constraints();
    UpperBoxType combined(component_count + constraint_count);
    Bool first = true;
    for (SizeType index = 0; index != sources.size(); ++index) {
        // Retain only the reference and the current transformed source, rather
        // than a second full collection of normalised function models.
        const auto source = index == 0u ? result : normalise(sources[index]);
        const auto source_map = join(join(source.state_function(), source.time_function()), source.dwell_time_function());
        auto image = factory.create_zeros(component_count + constraint_count, common_domain);
        for (SizeType k = 0; k != component_count; ++k)
            image.set(k, index == 0u ? factory.create_zero(common_domain) : source_map.get(k) - reference_map.get(k));
        for (SizeType k = 0; k != constraint_count; ++k)
            image.set(component_count + k, result.constraint(k).function());
        const auto bounds = constrained_image_bounds(source, image);
        if (definitely(bounds.is_empty()))
            continue;
        combined = first ? bounds : hull(combined, bounds);
        first = false;
    }
    if (first) {
        result.new_zero_parameter_constraint(factory.create_constant(common_domain, Real(1)));
        return result;
    }

    // Retain reference predicates with bounds covering every feasible source
    // alternative. The contractor uses each source's full constraint system;
    // it does not require matching constraint counts or ordering.
    for (SizeType k = 0; k != constraint_count; ++k) {
        const auto bounds = cast_exact_box(UpperBoxType({combined[component_count + k]}))[0];
        result._constraints[k] = ValidatedConstraintModel(bounds.lower_bound(), result.constraint(k).function(), bounds.upper_bound());
    }
    List<SizeType> residual_components;
    for (SizeType k = 0; k != component_count; ++k) {
        if (!same(combined[k], UpperIntervalType(0, 0)))
            residual_components.append(k);
    }
    ExactBoxType residual_domain(residual_components.size());
    for (SizeType k = 0; k != residual_components.size(); ++k)
        residual_domain[k] = cast_exact_box(UpperBoxType({combined[residual_components[k]]}))[0];

    result._domain = product(common_domain, residual_domain);
    result._reduced_domain = result._domain;
    result._state_function = embed(result._state_function, residual_domain);
    result._time_function = embed(result._time_function, residual_domain);
    result._dwell_time_function = embed(result._dwell_time_function, residual_domain);
    for (auto& constraint: result._constraints)
        constraint.function() = embed(constraint.function(), residual_domain);
    for (SizeType k = 0; k != residual_components.size(); ++k) {
        const auto correction = factory.create_coordinate(result._domain, parameter_count + k);
        const auto component = residual_components[k];
        if (component < state_dimension)
            result._state_function[component] = result._state_function.get(component) + correction;
        else if (component == state_dimension)
            result._time_function = result._time_function + correction;
        else
            result._dwell_time_function = result._dwell_time_function + correction;
    }
    result._variable_kinds.concatenate(List<EnclosureVariableKind>(residual_components.size(), EnclosureVariableKind::ERROR));
    result.reduce();
    return result;
}

} // namespace Ariadne
