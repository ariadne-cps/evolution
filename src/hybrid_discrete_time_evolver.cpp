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

#include <stdexcept>
#include <utility>
#include <vector>

#include "evolution/hybrid_discrete_time_predicate_evaluator.hpp"
#include "evolution/hybrid_discrete_time_step_image.hpp"
#include "evolution/hybrid_discrete_time_step_planner.hpp"
#include "hybrid_discrete_time_inclusion_capability.hpp"
#include "hybrid_discrete_time_predicate_function.hpp"

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

// LOWER may not independently select a parameter child of a
// witness-bearing outer enclosure, in line with hybrid continuous LOWER.

// Reuse the same complete branch planner as validated UPPER, but require a
// single branch whose precondition is certified for the entire source.
// Unlike UPPER, UNKNOWN evidence never yields a LOWER witness. In particular
// we do not impose an arbitrary priority on competing transitions.
// Each branch that is certain on the WHOLE current LOWER enclosure
// carries a real reachable source subset to a legal successor.
// Aggregate all such branches, preserving provenance and without ever
// splitting the numerical enclosure or prioritising one event.
List<HybridDiscreteTimePlannedBranch> certified_lower_branches(
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
    List<HybridDiscreteTimePlannedBranch> certain;
    for (const auto& candidate: plan.branches) {
        if (!candidate.certain)
            continue;
        if (candidate.kind == HybridDiscreteTimeBranchKind::TRANSITION &&
            !automaton.has_equation(candidate.target) &&
            !automaton.has_inclusion(candidate.target))
            continue;
        certain.append(candidate);
    }
    // UNKNOWN guards/invariants never become a certified LOWER branch.
    return certain;
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
    // Prove nonempty source once; no repeated witness selection in time.
    if (!has_certified_lower_source(current.continuous_set()))
        throw std::invalid_argument(
            "Discrete hybrid LOWER requires a certified feasible initial source");

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
    std::vector<HybridEnclosure> frontier{current};
    Integer index = termination.initial_index();
    SizeType nodes = 1u;
    while (index < termination.termination_index()) {
        const Integer successor_index = index + Integer(1);
        std::vector<HybridEnclosure> next_frontier;
        // Every alternative uses the complete original source and the
        // same numerical map as UPPER. No subboxes or input sampling.
        for (const HybridEnclosure& source: frontier) {
            const auto choices =
                certified_lower_branches(this->system(), source, index);
            for (const auto& choice: choices) {
                if (nodes >= this->configuration().maximum_nodes())
                    throw std::length_error("Discrete hybrid LOWER evolution exceeds node budget");
                HybridEnclosure successor =
                    HybridDiscreteTimeStepImage::apply(this->system(), source, choice);
                ++nodes;
                check_index(successor, successor_index);
                orbit.adjoin_reach(successor);
                if (successor_index < termination.termination_index())
                    orbit.adjoin_intermediate(successor);
                next_frontier.push_back(std::move(successor));
            }
        }
        if (next_frontier.empty())
            throw std::invalid_argument(
                "Discrete hybrid LOWER has no certain continuation to the requested horizon");
        frontier = std::move(next_frontier);
        index = successor_index;
    }
    for (const HybridEnclosure& final: frontier)
        orbit.adjoin_final(final);
    return orbit;
}

} // namespace Ariadne
