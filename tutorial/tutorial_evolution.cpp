#include <iostream>
#include <stdexcept>

#include "ariadne-evolution.hpp"

using namespace Ariadne;

int main() {
    // Continuous-time system: a one-dimensional stable vector field.
    RealVariable x("x");
    VectorField decay({dot(x)=-x});
    VectorFieldSimulator simulator(decay);
    simulator.configuration().set_step_size(0.05);
    simulator.configuration().set_num_subdivisions(0);
    simulator.configuration().insert_subspace(x,1);

    RealExpressionBoundedConstraintSet initial({0.9_dec<=x<=1.1_dec});
    auto simulated = simulator.orbit(initial,0.1_dec);
    const auto ncurves = simulated.curves().size();
    if (ncurves != 1u) {
        throw std::runtime_error("Expected one short continuous trajectory");
    }
    std::cout << "Continuous evolution: " << ncurves << " trajectory\n";

    // Discrete-time system: one step of a contracting iterated map.
    IteratedMap contraction({next(x)=0.5_dec*x});
    IteratedMapEvolver evolver(contraction);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    auto discrete = evolver.orbit(initial,IteratedMap::TimeType(1),Semantics::UPPER);
    if (discrete.final().size() == 0u) {
        throw std::runtime_error("Discrete orbit has no final enclosure");
    }
    std::cout << "Discrete evolution: " << discrete.final().size()
              << " final enclosure(s)\n";

    // Discrete-only hybrid system: the same two-mode automaton as the
    // Python tutorial. The simulation budget counts time and jumps.
    const DiscreteLocation idle(0), active(1);
    const DiscreteEvent activate("activate"), deactivate("deactivate");
    HybridAutomaton automaton;
    automaton.new_mode(idle);
    automaton.new_mode(active);
    automaton.new_transition(idle,activate,active);
    automaton.new_transition(active,deactivate,idle);

    if (!automaton.has_transition(idle,activate) ||
        automaton.target(active,deactivate) != idle) {
        throw std::runtime_error("Hybrid transition relation is incorrect");
    }
    HybridSimulator hybrid_simulator(automaton);
    HybridTime budget(1.0_x,2);
    if (budget.discrete_time() != Integer(2)) {
        throw std::runtime_error("Incorrect hybrid jump budget");
    }
    std::cout << "Hybrid evolution: " << automaton.locations().size()
              << " modes, budget " << budget << "\n";
    (void)hybrid_simulator;
    return 0;
}
