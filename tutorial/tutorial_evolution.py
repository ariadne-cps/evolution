#!/usr/bin/env python3
"""Bounded discrete and hybrid Evolution workflows using installed pyariadne_evolution."""

from pyariadne_evolution import (
    DiscreteEvent, DiscreteLocation, ContinuousTimeHybridAutomaton, HybridContinuousTimeEvolver, ApproximateTag,
    HybridTime, Integer, Rational, Real, Semantics, RealVariable,
    IteratedMap, DiscreteTimeEvolver, ValidatedTag, RealExpressionBoundedConstraintSet,
    next, dec_,
)


def main():
    x = RealVariable("x")
    contraction = IteratedMap([next(x) << x / 2])
    initial = RealExpressionBoundedConstraintSet(
        [(dec_(0.9) <= x) & (x <= dec_(1.1))], [])
    discrete_evolver = DiscreteTimeEvolver[ValidatedTag](contraction)
    discrete_evolver.configuration().set_maximum_enclosure_radius(1.0)
    orbit = discrete_evolver.orbit(initial, Integer(1), Semantics.UPPER)
    print(f"DiscreteTimeEvolver[ValidatedTag]: {len(orbit.final())} final enclosure(s)")

    idle, active = DiscreteLocation(0), DiscreteLocation(1)
    activate = DiscreteEvent("activate")
    deactivate = DiscreteEvent("deactivate")
    automaton = ContinuousTimeHybridAutomaton()
    automaton.new_mode(idle)
    automaton.new_mode(active)
    automaton.new_transition(idle, activate, active)
    automaton.new_transition(active, deactivate, idle)
    assert automaton.has_transition(idle, activate)
    assert automaton.target(active, deactivate) == idle

    simulator = HybridContinuousTimeEvolver[ApproximateTag](automaton)
    horizon = HybridTime(Real(Rational(1)), Integer(2))
    assert horizon.discrete_time() == Integer(2)
    assert Semantics.UPPER != Semantics.LOWER
    print(f"Hybrid simulation model: {automaton.number_of_modes()} modes")
    print(f"Hybrid budget: {horizon.continuous_time()}, {horizon.discrete_time()}")
    print(f"Simulator: {type(simulator).__name__}")


if __name__ == "__main__":
    main()
