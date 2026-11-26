/***************************************************************************
 *            src/hybrid_discrete_time_evolver.cpp
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

/** @file src/hybrid_discrete_time_evolver.cpp
 *  @brief Public deterministic discrete hybrid UPPER evolver.
 */

#include "evolution/hybrid_discrete_time_evolver.hpp"
#include "evolution/hybrid_discrete_time_step_image.hpp"
#include "evolution/hybrid_discrete_time_predicate_evaluator.hpp"
#include "evolution/hybrid_discrete_time_step_planner.hpp"
#include "hybrid_discrete_time_predicate_function.hpp"
#include "hybrid_discrete_time_inclusion_capability.hpp"

#include <deque>
#include <optional>
#include <stdexcept>
#include <utility>

namespace Ariadne {
namespace {

// This is an existential admission certificate, not the negation of
// Enclosure::is_empty(): that one-sided predicate can prove emptiness
// but its indeterminate result cannot establish a LOWER witness.
//
// A nonempty exact parameter domain provides an initial witness when
// no constraints are present. With constraints, certify either the entire
// reduced domain or one explicit point in it. Both checks use validated
// ranges: an UNKNOWN constraint result never becomes a LOWER witness.
// This is a sufficient test, not a general existential solver.
Bool has_certified_lower_source_on_box(
    const Enclosure& source, const ExactBoxType& box) {
    if (!definitely(!box.is_empty()) ||
        !definitely(Ariadne::subset(box, source.reduced_domain())))
        return false;
    if (source.number_of_constraints() == 0u)
        return true;

    const auto function = source.constraint_function();
    const ExactBoxType bounds = source.constraint_bounds();
    if (definitely(Ariadne::subset(Ariadne::apply(function, box), bounds)))
        return true;
    if (!definitely(box.is_bounded()))
        return false;

    const auto midpoint = box.midpoint();
    const ExactBoxType point(box.dimension(), [&](SizeType i) {
        return ExactIntervalType(midpoint[i]);
    });
    return definitely(Ariadne::subset(Ariadne::apply(function, point), bounds));
}

Bool has_certified_lower_source(const Enclosure& source) {
    return has_certified_lower_source_on_box(source, source.reduced_domain());
}

constexpr SizeType lower_refinement_depth = 4u;

struct ParameterBoxSearch {
    ExactBoxType box;
    SizeType depth;
};

Void enqueue_parameter_halves(
    const ExactBoxType& box, SizeType depth, std::deque<ParameterBoxSearch>& pending) {
    // Rotate the starting axis with depth: always splitting the first
    // nondegenerate parameter can miss a decidable guard on later axes.
    // This changes only search completeness, not branch certification.
    const SizeType dimension = box.dimension();
    for (SizeType offset = 0u; offset < dimension; ++offset) {
        const SizeType axis = (depth + offset) % dimension;
        if (!(box[axis].lower_bound().raw() < box[axis].upper_bound().raw()))
            continue;
        const auto children = box.split(axis);
        if (!same(children.first, box) && definitely(!children.first.is_empty()))
            pending.push_back({children.first, depth + 1u});
        if (!same(children.second, box) && definitely(!children.second.is_empty()))
            pending.push_back({children.second, depth + 1u});
        return;
    }
}

Bool has_refined_lower_witness(const Enclosure& source) {
    std::deque<ParameterBoxSearch> pending{{source.reduced_domain(), 0u}};
    while (!pending.empty()) {
        const ParameterBoxSearch candidate = pending.front();
        pending.pop_front();
        if (has_certified_lower_source_on_box(source, candidate.box))
            return true;
        if (candidate.depth < lower_refinement_depth)
            enqueue_parameter_halves(candidate.box, candidate.depth, pending);
    }
    return false;
}

// Reuse the same complete branch planner as validated UPPER, but require a
// single branch whose precondition is certified for the entire source.
// Unlike UPPER, UNKNOWN evidence never yields a LOWER witness. In particular
// we do not impose an arbitrary priority on competing transitions.
std::optional<HybridDiscreteTimePlannedBranch> unique_certified_lower_branch(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source, const Integer& index) {
    if (!automaton.has_mode(source.location()) ||
        (!automaton.has_equation(source.location()) &&
         !automaton.has_inclusion(source.location())))
        throw std::invalid_argument(
            "Discrete hybrid LOWER requires an equation or bounded-input inclusion mode");
    const auto evaluate = [&](DiscreteLocation, DiscreteEvent,
                              DiscreteTimePredicateRole,
                              const ContinuousPredicate& predicate) {
        return HybridDiscreteTimePredicateEvaluator::classify(
            automaton, source, predicate);
    };
    const auto plan = HybridDiscreteTimeStepPlanner::plan(
        automaton, source.location(), index, evaluate);
    if (plan.branches.size() != 1u || !plan.branches.front().certain)
        return std::nullopt;
    const auto branch = plan.branches.front();
    if (branch.kind == HybridDiscreteTimeBranchKind::TRANSITION &&
        !automaton.has_equation(branch.target) &&
        !automaton.has_inclusion(branch.target))
        return std::nullopt;
    return branch;
}

// Certify a predicate on an exact parameter subbox using the original
// validated state patch. Boundary/undecided evidence stays UNKNOWN.
DiscreteTimePredicateEvidence classify_on_box(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source, const ContinuousPredicate& predicate,
    const ExactBoxType& box) {
    const auto indicator =
        EvolutionDetail::discrete_source_indicator(automaton, source, predicate);
    const auto args = source.continuous_set().argument_function(
        indicator.arguments, 0u);
    const auto range = Ariadne::apply(compose(indicator.function, args), box);
    if (definitely(range.lower_bound() > 0))
        return DiscreteTimePredicateEvidence::TRUE;
    if (definitely(range.upper_bound() < 0))
        return DiscreteTimePredicateEvidence::FALSE;
    return DiscreteTimePredicateEvidence::UNKNOWN;
}

std::optional<HybridDiscreteTimePlannedBranch> unique_lower_branch_on_box(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source, const Integer& index,
    const ExactBoxType& box) {
    const auto evaluate = [&](DiscreteLocation, DiscreteEvent,
                              DiscreteTimePredicateRole,
                              const ContinuousPredicate& predicate) {
        return classify_on_box(automaton, source, predicate, box);
    };
    const auto plan = HybridDiscreteTimeStepPlanner::plan(
        automaton, source.location(), index, evaluate);
    if (plan.branches.size() != 1u || !plan.branches.front().certain)
        return std::nullopt;
    const auto branch = plan.branches.front();
    if (branch.kind == HybridDiscreteTimeBranchKind::TRANSITION &&
        !automaton.has_equation(branch.target) &&
        !automaton.has_inclusion(branch.target))
        return std::nullopt;
    return branch;
}

// A certified exact subbox becomes parameter constraints on a COPY of
// the source, not a restriction of the function-patch domain. The latter
// produced a reproducible crash on a refined reset in v7.05/v7.06.
// Composition preserves all existing and newly added constraints.
HybridEnclosure constrain_lower_source(
    const HybridEnclosure& source, const ExactBoxType& box) {
    HybridEnclosure restricted = source;
    Enclosure set = source.continuous_set();
    const auto domain = set.domain();
    for (SizeType i = 0u; i < box.dimension(); ++i) {
        if (same(box[i], domain[i]))
            continue;
        const auto coordinate =
            EffectiveScalarMultivariateFunction::coordinate(box.dimension(), i);
        // FloatDP box endpoints are exact binary values, but are not
        // implicitly convertible to ValidatedNumber. Convert through an
        // exact dyadic number before constructing the interval constraint.
        const ValidatedNumber lower(ExactNumber(Dyadic(box[i].lower_bound())));
        const ValidatedNumber upper(ExactNumber(Dyadic(box[i].upper_bound())));
        set.new_parameter_constraint(ValidatedConstraint(lower, coordinate, upper));
    }
    restricted.replace_continuous_set(set);
    return restricted;
}

// Both a parameter witness and the unique branch are certified on the
// SAME subbox; the ensuing reset retains that box as parameter constraints.
struct CertifiedLowerStep {
    HybridEnclosure source;
    HybridDiscreteTimePlannedBranch branch;
};

CertifiedLowerStep certified_lower_step(
    const DiscreteTimeHybridAutomatonInterface& automaton,
    const HybridEnclosure& source, const Integer& index) {
    const Bool inclusion = automaton.has_inclusion(source.location());
    if (!automaton.has_mode(source.location()) ||
        (!automaton.has_equation(source.location()) && !inclusion))
        throw std::invalid_argument(
            "Discrete hybrid LOWER requires a supported equation or inclusion mode");

    if (has_certified_lower_source(source.continuous_set())) {
        const auto branch = unique_certified_lower_branch(automaton, source, index);
        if (branch)
            return {source, *branch};
    }

    // Inclusion LOWER must preserve ALL current inputs, not select a
    // favourable input subbox to resolve an uncertain event.
    if (inclusion)
        throw std::invalid_argument(
            "Discrete hybrid inclusion LOWER cannot certify a unique branch over the full input box");

    std::deque<ParameterBoxSearch> pending{{source.continuous_set().reduced_domain(), 0u}};
    while (!pending.empty()) {
        const ParameterBoxSearch candidate = pending.front();
        pending.pop_front();
        if (has_certified_lower_source_on_box(source.continuous_set(), candidate.box)) {
            const auto branch = unique_lower_branch_on_box(
                automaton, source, index, candidate.box);
            if (branch)
                return {constrain_lower_source(source, candidate.box), *branch};
        }
        if (candidate.depth < lower_refinement_depth)
            enqueue_parameter_halves(candidate.box, candidate.depth, pending);
    }
    throw std::invalid_argument(
        "Discrete hybrid LOWER could not certify a feasible subbox with one branch");
}

} // namespace

HybridDiscreteTimeTerminationCriterion::HybridDiscreteTimeTerminationCriterion(
    Integer initial_index, Integer termination_index)
    : _initial_index(std::move(initial_index)),
      _termination_index(std::move(termination_index)) {
    if (_termination_index < _initial_index)
        throw std::invalid_argument("Discrete hybrid termination index precedes initial index");
}

HybridDiscreteTimeEvolverConfiguration<ValidatedTag>::HybridDiscreteTimeEvolverConfiguration()
    : _maximum_nodes(10000u) {}

Void HybridDiscreteTimeEvolverConfiguration<ValidatedTag>::set_maximum_nodes(SizeType value) {
    if (value == 0u)
        throw std::invalid_argument("Discrete hybrid maximum node budget must be positive");
    _maximum_nodes = value;
}

OutputStream& HybridDiscreteTimeEvolverConfiguration<ValidatedTag>::_write(OutputStream& os) const {
    return os << "HybridDiscreteTimeEvolverConfiguration<ValidatedTag>(maximum_nodes="
              << _maximum_nodes << ")";
}

HybridDiscreteTimeEvolver<ValidatedTag>::HybridDiscreteTimeEvolver(const SystemType& system)
    : _system(system.clone()), _configuration(std::make_shared<ConfigurationType>()) {}

HybridDiscreteTimeEvolver<ValidatedTag>::HybridDiscreteTimeEvolver(
    const HybridDiscreteTimeEvolver<ValidatedTag>& other)
    : _system(other._system->clone()),
      _configuration(std::make_shared<ConfigurationType>()) {
    _configuration->set_maximum_nodes(other.configuration().maximum_nodes());
}

HybridDiscreteTimeEvolver<ValidatedTag>* HybridDiscreteTimeEvolver<ValidatedTag>::clone() const {
    return new HybridDiscreteTimeEvolver<ValidatedTag>(*this);
}

HybridDiscreteTimeUpperExploration HybridDiscreteTimeEvolver<ValidatedTag>::explore(
    const EnclosureType& initial, const TerminationType& termination, Semantics semantics) const {
    if (semantics != Semantics::UPPER)
        throw std::invalid_argument("Discrete hybrid validated LOWER evolution is not implemented");
    return HybridDiscreteTimeUpperScheduler::evolve(
        this->system(), initial, termination.initial_index(),
        termination.termination_index(), this->configuration().maximum_nodes());
}

HybridDiscreteTimeUpperOrbitProjection HybridDiscreteTimeEvolver<ValidatedTag>::orbit_with_blocked(
    const EnclosureType& initial, const TerminationType& termination, Semantics semantics) const {
    return HybridDiscreteTimeUpperScheduler::project_orbit(
        this->explore(initial, termination, semantics));
}

Orbit<HybridEnclosure> HybridDiscreteTimeEvolver<ValidatedTag>::orbit(
    const EnclosureType& initial, const TerminationType& termination, Semantics semantics) const {
    if (semantics == Semantics::UPPER)
        return this->orbit_with_blocked(initial, termination, semantics).orbit;
    if (semantics != Semantics::LOWER)
        throw std::invalid_argument("Unknown discrete hybrid evolution semantics");

    // Both semantics use the same numerical step map and entire
    // inclusion input box. Only discrete branch certification differs.
    const DiscreteLocation location = initial.location();
    const Bool equation_source = this->system().has_equation(location);
    const Bool inclusion_source = this->system().has_inclusion(location);
    if (!this->system().has_mode(location) ||
        (!equation_source && !inclusion_source))
        throw std::invalid_argument(
            "Discrete hybrid LOWER requires an equation or bounded-input inclusion source");
    if (initial.state_space() != this->system().state_space(location))
        throw std::invalid_argument("Discrete hybrid LOWER initial state space mismatch");

    HybridEnclosure current = initial;
    if (equation_source) {
        const auto& equation = this->system().equation(location);
        if (initial.auxiliary_space().dimension() != 0u &&
            initial.auxiliary_space() != equation.auxiliary_space())
            throw std::invalid_argument("Discrete hybrid LOWER initial auxiliary space mismatch");
        current.set_auxiliary(equation.auxiliary_space().variables(),
                              equation.auxiliary_function());
    } else {
        // The current full input box is already required by a(x,u[n])
        // at zero horizon and when certifying auxiliary guards.
        current = EvolutionDetail::prepare_discrete_inclusion_source(this->system(), initial);
    }
    if (termination.initial_index() == termination.termination_index() &&
        !has_refined_lower_witness(current.continuous_set()))
        throw std::invalid_argument(
            "Discrete hybrid LOWER requires a certified feasible initial parameter witness");

    const auto check_index = [](const HybridEnclosure& set, const Integer& index) {
        const Real expected_time(index);
        const UpperIntervalType range = set.time_range();
        if (!definitely(range.lower_bound() >= expected_time) ||
            !definitely(range.upper_bound() <= expected_time))
            throw std::invalid_argument(
                "Discrete hybrid LOWER enclosure time does not match the absolute index");
    };
    check_index(initial, termination.initial_index());

    OrbitType orbit(current);
    Integer index = termination.initial_index();
    SizeType nodes = 1u;
    while (index < termination.termination_index()) {
        if (nodes >= this->configuration().maximum_nodes())
            throw std::length_error("Discrete hybrid LOWER evolution exceeds node budget");
        const auto certified = certified_lower_step(this->system(), current, index);
        // The one-tick map is the same as UPPER. Inclusion sources
        // retain their full input box; equations retain their existing
        // state-parameter guard refinement. No input is sampled.
        current = HybridDiscreteTimeStepImage::apply(
            this->system(), certified.source, certified.branch);
        index = index + Integer(1);
        ++nodes;
        check_index(current, index);
        orbit.adjoin_reach(current);
        if (index < termination.termination_index())
            orbit.adjoin_intermediate(current);
    }
    orbit.adjoin_final(current);
    return orbit;
}

} // namespace Ariadne
