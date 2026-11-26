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
Bool has_certified_lower_source(const Enclosure& source) {
    const ExactBoxType domain = source.reduced_domain();
    if (!definitely(!domain.is_empty()))
        return false;
    if (source.number_of_constraints() == 0u)
        return true;

    const auto constraint_function = source.constraint_function();
    const ExactBoxType constraints = source.constraint_bounds();
    const auto range = Ariadne::apply(constraint_function, domain);
    if (definitely(Ariadne::subset(range, constraints)))
        return true;

    // A bounded box has an exactly representable midpoint. Proving that
    // every constraint contains its validated value at this singleton
    // proves the existence of at least one feasible initial parameter.
    // Preserve the original constrained enclosure for the actual orbit:
    // the witness certificate must not weaken or delete its constraints.
    if (!definitely(domain.is_bounded()))
        return false;
    const auto centre = domain.midpoint();
    const ExactBoxType point(domain.dimension(), [&](SizeType i) {
        return ExactIntervalType(centre[i]);
    });
    const auto point_values = Ariadne::apply(constraint_function, point);
    return definitely(Ariadne::subset(point_values, constraints));
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

    // A single deterministic mode with no events has exactly one feasible
    // execution from each nonempty initial state. Its validated image map
    // is witness-preserving, as in ordinary discrete validated LOWER.
    // Do not infer a LOWER witness from an UPPER guard candidate: the
    // event-bearing and bounded-input cases have separate proof obligations.
    const DiscreteLocation location = initial.location();
    if (!this->system().has_mode(location) ||
        !this->system().has_equation(location) ||
        !this->system().events(location).empty())
        throw std::invalid_argument(
            "Discrete hybrid LOWER currently requires an event-free difference equation mode");
    if (initial.state_space() != this->system().state_space(location))
        throw std::invalid_argument("Discrete hybrid LOWER initial state space mismatch");
    const auto& equation = this->system().equation(location);
    if (initial.auxiliary_space().dimension() != 0u &&
        initial.auxiliary_space() != equation.auxiliary_space())
        throw std::invalid_argument("Discrete hybrid LOWER initial auxiliary space mismatch");

    // Enclosure::is_empty() is an emptiness semi-decision and cannot
    // certify a nonempty source by negation. Admit only a validated
    // existential witness: a nonempty exact parameter domain, with
    // all declared constraints provably satisfied on it when present.
    if (!has_certified_lower_source(initial.continuous_set()))
        throw std::invalid_argument(
            "Discrete hybrid LOWER requires a certified feasible initial parameter domain");

    const auto check_index = [](const HybridEnclosure& set, const Integer& index) {
        const Real expected_time(index);
        const UpperIntervalType range = set.time_range();
        if (!definitely(range.lower_bound() >= expected_time) ||
            !definitely(range.upper_bound() <= expected_time))
            throw std::invalid_argument(
                "Discrete hybrid LOWER enclosure time does not match the absolute index");
    };
    check_index(initial, termination.initial_index());

    // At zero horizon, expose the same equation auxiliary mapping as
    // after progression. This is a view of the current state, not a
    // new witness or an additional parameter.
    HybridEnclosure current = initial;
    current.set_auxiliary(equation.auxiliary_space().variables(),
                          equation.auxiliary_function());
    OrbitType orbit(current);
    Integer index = termination.initial_index();
    SizeType nodes = 1u;
    while (index < termination.termination_index()) {
        if (nodes >= this->configuration().maximum_nodes())
            throw std::length_error("Discrete hybrid LOWER evolution exceeds node budget");
        const HybridDiscreteTimePlannedBranch step{
            HybridDiscreteTimeBranchKind::PROGRESSION, location, std::nullopt, true};
        current = HybridDiscreteTimeStepImage::apply(this->system(), current, step);
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
