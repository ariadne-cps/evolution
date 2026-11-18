/***************************************************************************
 *            enclosure_recombiner.cpp
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

/** @file enclosure_recombiner.cpp
 *  @brief Reference strategy and handle for continuous enclosure recombination.
 */

#include "evolution/enclosure_recombiner.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include "function/formula.hpp"
#include "function/formula.tpl.hpp"
#include "function/function.hpp"
#include "function/procedure.hpp"
#include "utility/macros.hpp"

#include "constraint-solving/constraint_solver.hpp"

namespace Ariadne {
namespace {

using RecombinationClock = std::chrono::steady_clock;

std::uint64_t phase_elapsed_ns(RecombinationClock::time_point start) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(RecombinationClock::now() - start).count());
}

// In ordinary (unprofiled) calls this guard never reads a clock.
class ScopedRecombinationTimer {
    std::uint64_t* _elapsed;
    RecombinationClock::time_point _start{};

  public:
    explicit ScopedRecombinationTimer(std::uint64_t* elapsed) : _elapsed(elapsed) {
        if (_elapsed)
            _start = RecombinationClock::now();
    }
    ~ScopedRecombinationTimer() {
        if (_elapsed)
            *_elapsed += phase_elapsed_ns(_start);
    }
    ScopedRecombinationTimer(const ScopedRecombinationTimer&) = delete;
    ScopedRecombinationTimer& operator=(const ScopedRecombinationTimer&) = delete;
};

// Bound an image on the feasible parameter set by contracting its graph.
// Keeping each source's own constraints here avoids intersecting alternatives.
UpperBoxType constrained_image_bounds(const Enclosure& source,
                                      const ValidatedVectorMultivariateFunctionPatch& image,
                                      EnclosureRecombinationProfile* profile,
                                      EnclosureRecombinationSourceProfile* source_profile) {
    const auto parameter_domain = source.reduced_domain();
    UpperBoxType bounds = [&] {
        ScopedRecombinationTimer timer(profile ? &profile->image_bounds_ns : nullptr);
        ScopedRecombinationTimer source_timer(source_profile ? &source_profile->image_bounds_ns : nullptr);
        UpperBoxType image_bounds = apply(image, parameter_domain);
        // A source predicate g in C also bounds any image f by C+range(f-g).
        // In particular, this retains equalities shared with the reference even
        // when an oblique constraint cannot contract the parameter box itself.
        for (SizeType k = 0; k != image.result_size(); ++k) {
            const auto component = image.get(k);
            for (const auto& constraint: source.constraint_models()) {
                image_bounds[k] = intersection(image_bounds[k], UpperIntervalType(constraint.bounds()) + apply(component - constraint.function(), parameter_domain));
            }
        }
        return image_bounds;
    }();
    if (definitely(bounds.is_empty()))
        return bounds;
    ARIADNE_PRECONDITION_MSG(definitely(bounds.is_bounded()), "Unbounded recombination image: " << bounds);
    if (source.number_of_constraints() == 0u || parameter_domain.size() == 0u)
        return bounds;

    const auto graph_setup_start = (profile || source_profile) ? RecombinationClock::now() : RecombinationClock::time_point{};
    // Singleton image bounds are already certified by the range calculation.
    // Adding a fixed graph variable cannot tighten that image,
    // but would enlarge every embedded model and the solver's shaving loops.
    List<SizeType> active_components;
    for (SizeType k = 0; k != image.result_size(); ++k) {
        if (bounds[k].lower_bound().raw() != bounds[k].upper_bound().raw())
            active_components.append(k);
    }
    if (source_profile)
        source_profile->graph_active_components = active_components.size();
    if (active_components.empty()) {
        if (profile)
            profile->graph_setup_ns += phase_elapsed_ns(graph_setup_start);
        if (source_profile)
            source_profile->graph_setup_ns += phase_elapsed_ns(graph_setup_start);
        return bounds;
    }
    ExactBoxType image_domain(active_components.size());
    const auto exact_bounds = cast_exact_box(bounds);
    for (SizeType k = 0; k != active_components.size(); ++k)
        image_domain[k] = exact_bounds[active_components[k]];
    const auto graph_domain = product(parameter_domain, image_domain);
    if (source_profile) {
        source_profile->graph_dimension = graph_domain.size();
        source_profile->graph_constraints = source.number_of_constraints() + image.result_size();
    }
    const auto& factory = source.configuration().function_factory();
    const auto parameter_count = parameter_domain.size();
    const auto constraint_count = source.number_of_constraints();
    auto graph = factory.create_zeros(constraint_count + image.result_size(), graph_domain);
    ExactBoxType graph_bounds(graph.result_size());
    const auto restrict_to_parameters = [&](const ValidatedScalarMultivariateFunctionPatch& function) {
        // Normalised sources already use parameter_domain. Restricting a Taylor
        // patch to its own domain still performs a full identity composition.
        return same(function.domain(), parameter_domain) ? function : restriction(function, parameter_domain);
    };
    for (SizeType k = 0; k != constraint_count; ++k) {
        graph.set(k, embed(restrict_to_parameters(source.constraint(k).function()), image_domain));
        graph_bounds[k] = source.constraint(k).bounds();
    }
    SizeType active_index = 0u;
    for (SizeType k = 0; k != image.result_size(); ++k) {
        const auto component = embed(restrict_to_parameters(image.get(k)), image_domain);
        if (active_index != active_components.size() && active_components[active_index] == k) {
            graph.set(constraint_count + k,
                      component - factory.create_coordinate(graph_domain, parameter_count + active_index));
            graph_bounds[constraint_count + k] = ExactIntervalType(0, 0);
            ++active_index;
        } else {
            // Retain the implied predicate even for a fixed image: it can
            // help contract the parameters and therefore other image bounds.
            graph.set(constraint_count + k, component);
            graph_bounds[constraint_count + k] = exact_bounds[k];
        }
    }
    UpperBoxType feasible_graph(graph_domain);
    ConstraintSolver contractor;
    std::vector<ConstraintPropagationConstraint> prepared;
    prepared.reserve(graph.result_size());
    for (SizeType k = 0; k != graph.result_size(); ++k) {
        const ValidatedScalarMultivariateFunction function = graph.get(k);
        prepared.push_back({function, graph_bounds[k], {}, false, false, std::make_shared<const ValidatedProcedure>(function)});
    }
    if (profile)
        profile->graph_setup_ns += phase_elapsed_ns(graph_setup_start);
    if (source_profile)
        source_profile->graph_setup_ns += phase_elapsed_ns(graph_setup_start);
    Bool graph_is_empty;
    ConstraintPropagationStatistics statistics;
    {
        ScopedRecombinationTimer timer(profile ? &profile->graph_reduce_ns : nullptr);
        ScopedRecombinationTimer source_timer(source_profile ? &source_profile->graph_reduce_ns : nullptr);
        graph_is_empty = contractor.propagate(feasible_graph, prepared, false, statistics);
    }
    if (source_profile) {
        source_profile->hull_rounds = statistics.hull_rounds;
        source_profile->hull_contractor_calls = statistics.hull_contractor_calls;
        source_profile->hull_effective = statistics.hull_effective;
        source_profile->shaving_rounds = statistics.shaving_rounds;
        source_profile->shaving_coordinate_attempts = statistics.shaving_coordinate_attempts;
        source_profile->shaving_function_evaluations = statistics.shaving_function_evaluations;
        source_profile->shaving_effective = statistics.shaving_effective;
        source_profile->monotone_rounds = statistics.monotone_rounds;
        source_profile->monotone_effective = statistics.monotone_effective;
    }
    if (graph_is_empty || definitely(feasible_graph.is_empty())) {
        bounds[0] = UpperIntervalType(1, -1);
        return bounds;
    }
    for (SizeType k = 0; k != active_components.size(); ++k)
        bounds[active_components[k]] = feasible_graph[parameter_count + k];
    return bounds;
}

} // namespace

EnclosureRecombinationContext::EnclosureRecombinationContext(ExactDouble maximum_enclosure_radius)
    : _maximum_enclosure_radius(maximum_enclosure_radius) {
    ARIADNE_PRECONDITION_MSG(std::isfinite(maximum_enclosure_radius.get_d()) && maximum_enclosure_radius >= 0.0_x,
                             "The requested maximum enclosure radius must be finite and nonnegative.");
}

EnclosureRecombiner::EnclosureRecombiner(const EnclosureRecombinerInterface& recombiner)
    : Handle<const EnclosureRecombinerInterface>(recombiner.clone()) {}

EnclosureRecombinationResult EnclosureRecombiner::recombine(const List<Enclosure>& enclosures,
                                                            const EnclosureRecombinationContext& context) const {
    return this->_ptr->recombine(enclosures, context);
}

ReferenceEnclosureRecombiner* ReferenceEnclosureRecombiner::clone() const {
    return new ReferenceEnclosureRecombiner(*this);
}

OutputStream& ReferenceEnclosureRecombiner::_write(OutputStream& os) const {
    return os << "ReferenceEnclosureRecombiner()";
}

EnclosureRecombinationResult ReferenceEnclosureRecombiner::recombine(const List<Enclosure>& enclosures,
                                                                     const EnclosureRecombinationContext& context) const {
    return recombine(enclosures, context, nullptr);
}

EnclosureRecombinationResult ReferenceEnclosureRecombiner::recombine(const List<Enclosure>& enclosures,
                                                                     const EnclosureRecombinationContext& context,
                                                                     EnclosureRecombinationProfile* profile) const {
    return recombine(enclosures, context, profile, nullptr);
}

EnclosureRecombinationResult ReferenceEnclosureRecombiner::recombine(
    const List<Enclosure>& enclosures, const EnclosureRecombinationContext& context,
    EnclosureRecombinationProfile* profile,
    std::vector<EnclosureRecombinationSourceProfile>* source_profiles) const {
    if (profile)
        *profile = {};
    if (source_profiles) {
        source_profiles->clear();
        source_profiles->resize(enclosures.size());
        for (SizeType index = 0u; index != enclosures.size(); ++index)
            (*source_profiles)[index].source_index = index;
    }
    EnclosureRecombinationResult batch;
    if (enclosures.empty())
        return batch;
    ARIADNE_PRECONDITION_MSG(std::all_of(enclosures.begin(), enclosures.end(), [](const Enclosure& source) {
                                 return source._state_function.raw_pointer() != nullptr;
                             }),
                             "Recombination requires initialised enclosure models.");
    const auto state_dimension = enclosures.front().state_dimension();
    List<Enclosure> sources;
    sources.reserve(enclosures.size());
    List<SizeType> source_indices;
    source_indices.reserve(enclosures.size());
    SizeType parameter_count = 0u;
    {
        ScopedRecombinationTimer timer(profile ? &profile->filter_ns : nullptr);
        for (SizeType index = 0; index != enclosures.size(); ++index) {
            auto* individual = source_profiles ? &(*source_profiles)[index] : nullptr;
            ScopedRecombinationTimer source_timer(individual ? &individual->filter_ns : nullptr);
            auto source = enclosures[index];
            if (individual) {
                individual->input_parameters = source.number_of_parameters();
                individual->input_constraints = source.number_of_constraints();
            }
            ARIADNE_PRECONDITION_MSG(source.state_dimension() == state_dimension,
                                     "Recombination requires matching state dimensions.");
            // is_empty contracts the copy's parameter domain and only discards
            // sources whose infeasibility has actually been established.
            if (definitely(source.is_empty())) {
                if (individual)
                    individual->discarded = true;
                batch.discarded_sources.append(index);
                continue;
            }
            ARIADNE_PRECONDITION_MSG(source.reduced_domain().is_bounded(),
                                     "Recombination requires bounded parameter domains.");
            parameter_count = std::max(parameter_count, source.number_of_parameters());
            sources.append(source);
            source_indices.append(index);
        }
    }
    if (sources.empty())
        return batch;
    const auto discard_survivors = [&] {
        batch.discarded_sources.insert(batch.discarded_sources.end(), source_indices.begin(), source_indices.end());
        std::sort(batch.discarded_sources.begin(), batch.discarded_sources.end());
        return std::move(batch);
    };
    const auto finish = [&](Enclosure enclosure) {
        if (definitely(enclosure.reduced_domain().is_empty()))
            return discard_survivors();
        if (context.maximum_enclosure_radius())
            batch.radius_limit_satisfied = definitely(enclosure.radius() <= *context.maximum_enclosure_radius());
        batch.outputs.push_back({std::move(enclosure), std::move(source_indices)});
        return std::move(batch);
    };
    if (sources.size() == 1u)
        return finish(std::move(sources.front()));

    const auto& reference = sources.front();
    const auto auxiliary_formula = make_formula(reference.auxiliary_mapping());
    for (const auto& source: sources) {
        ARIADNE_PRECONDITION_MSG(identical(auxiliary_formula, make_formula(source.auxiliary_mapping())),
                                 "Recombination requires identical auxiliary mappings.");
    }
    const auto& factory = reference.configuration().function_factory();
    const ExactBoxType common_domain(parameter_count, ExactIntervalType(0, 1));
    const auto normalise = [&](const Enclosure& source,
                               EnclosureRecombinationSourceProfile* individual) {
        ScopedRecombinationTimer timer(profile ? &profile->normalise_ns : nullptr);
        ScopedRecombinationTimer source_timer(individual ? &individual->normalise_ns : nullptr);
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
        // Every exact source parameter in reduced maps from some u in [0,1]^m
        // by p_k = lower_k*(1-u_k) + upper_k*u_k. The rounded patch is an
        // enclosure of this exact affine map, not the map defining feasibility.
        // Thus the exact map remains in the source patch's valid domain even
        // if the *enclosing interval evaluation* extends slightly beyond it.
        // unchecked_compose bypasses the latter over-strict domain test; the
        // validated Taylor-model composition must still enclose evaluation
        // along the exact admissible affine map. See the constructive coverage
        // proof in doc/evolution_tree_design.md.
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

    Enclosure result = normalise(sources.front(),
                                source_profiles ? &(*source_profiles)[source_indices.front()] : nullptr);
    const auto reference_map = join(join(result.state_function(), result.time_function()), result.dwell_time_function());
    const auto component_count = reference_map.result_size();
    const auto constraint_count = result.number_of_constraints();
    UpperBoxType combined(component_count + constraint_count);
    Bool first = true;
    for (SizeType index = 0; index != sources.size(); ++index) {
        // Retain only the reference and the current transformed source, rather
        // than a second full collection of normalised function models.
        auto* individual = source_profiles ? &(*source_profiles)[source_indices[index]] : nullptr;
        const auto source = index == 0u ? result : normalise(sources[index], individual);
        const auto image = [&] {
            ScopedRecombinationTimer timer(profile ? &profile->image_setup_ns : nullptr);
            ScopedRecombinationTimer source_timer(individual ? &individual->image_setup_ns : nullptr);
            const auto source_map = join(join(source.state_function(), source.time_function()), source.dwell_time_function());
            auto residual_image = factory.create_zeros(component_count + constraint_count, common_domain);
            for (SizeType k = 0; k != component_count; ++k)
                residual_image.set(k, index == 0u ? factory.create_zero(common_domain) : source_map.get(k) - reference_map.get(k));
            for (SizeType k = 0; k != constraint_count; ++k)
                residual_image.set(component_count + k, result.constraint(k).function());
            return residual_image;
        }();
        const auto bounds = constrained_image_bounds(source, image, profile, individual);
        if (definitely(bounds.is_empty()))
            continue;
        combined = first ? bounds : hull(combined, bounds);
        first = false;
    }
    if (first) {
        return discard_survivors();
    }

    // Constructive witness for each feasible normalised source parameter u:
    // let r = source_map(u) - reference_map(u). Every component of r lies
    // in combined[0:component_count], and each reference predicate at u
    // lies in the corresponding combined constraint interval. Therefore
    // the result contains (state,time,dwell) = reference_map(u) + r.
    // Contractor reductions below must preserve this witness. Source-only
    // constraints restrict the bounds calculation, never the other alternatives.
    // See the constructive coverage proof in doc/evolution_tree_design.md.
    // The contractor uses each source's full constraint system and permits
    // different source constraint lists and ordering.
    {
        ScopedRecombinationTimer timer(profile ? &profile->output_assembly_ns : nullptr);
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
    }
    {
        ScopedRecombinationTimer timer(profile ? &profile->output_reduce_ns : nullptr);
        result.reduce();
    }
    return finish(std::move(result));
}

Enclosure recombine(const List<Enclosure>& enclosures) {
    ARIADNE_PRECONDITION_MSG(!enclosures.empty(), "Recombination requires a source enclosure.");
    auto result = ReferenceEnclosureRecombiner().recombine(enclosures);
    if (result.outputs.empty()) {
        auto empty = enclosures.front();
        if (!definitely(empty.is_empty()))
            empty.new_zero_parameter_constraint(EffectiveScalarMultivariateFunction::constant(empty.number_of_parameters(), Real(1)));
        return empty;
    }
    return std::move(result.outputs.front().enclosure);
}

} // namespace Ariadne
