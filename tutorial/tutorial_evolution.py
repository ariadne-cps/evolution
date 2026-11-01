#!/usr/bin/env python3
"""Short Evolution example using the types exposed by the Python bindings.

The discrete-only hybrid model matches the final section of the C++ tutorial.
The continuous VectorField and IteratedMap examples require System's C++
modelling API, which is not yet exposed as constructors in Python.
"""

from pyariadne import (
    DiscreteEvent, DiscreteLocation, HybridAutomaton, HybridSimulator,
    HybridTime, Integer, Rational, Real, Semantics,
    VectorFieldSimulator, IteratedMapEvolver,
)


def main():
    idle, active = DiscreteLocation(0), DiscreteLocation(1)
    activate = DiscreteEvent("activate")
    deactivate = DiscreteEvent("deactivate")

    automaton = HybridAutomaton()
    automaton.new_mode(idle)
    automaton.new_mode(active)
    automaton.new_transition(idle, activate, active)
    automaton.new_transition(active, deactivate, idle)

    if not automaton.has_transition(idle, activate):
        raise RuntimeError("Missing activation transition")
    if automaton.target(active, deactivate) != idle:
        raise RuntimeError("Incorrect deactivation transition")

    # Evolution owns the simulator and its hybrid time horizon.
    simulator = HybridSimulator(automaton)
    horizon = HybridTime(Real(Rational(1)), Integer(2))
    if horizon.discrete_time() != Integer(2):
        raise RuntimeError("Incorrect hybrid jump budget")
    if Semantics.UPPER == Semantics.LOWER:
        raise RuntimeError("Upper and lower semantics must differ")

    print(f"Hybrid evolution: {automaton.number_of_modes()} modes")
    print(f"Budget: continuous={horizon.continuous_time()}, "
          f"discrete={horizon.discrete_time()}")
    print(f"Simulator: {type(simulator).__name__}")
    print("Continuous and discrete evolver interfaces available:",
          VectorFieldSimulator.__name__, IteratedMapEvolver.__name__)


if __name__ == "__main__":
    main()
