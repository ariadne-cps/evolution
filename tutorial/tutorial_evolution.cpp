/***************************************************************************
 *            tutorial_evolution.cpp
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

/** @file tutorial/tutorial_evolution.cpp
 *  @brief Demonstrate continuous, discrete and hybrid evolution of dynamical systems.
 */

#include <iostream>
#include <stdexcept>

#include "function/function.hpp"
#include "hybrid/hybrid_set.hpp"
#include "integration/integrator.hpp"
#include "numeric/float_bounds.hpp"

#include "ariadne-evolution.hpp"

using namespace Ariadne;

int main() {
    // Continuous-time system: a one-dimensional stable vector field.
    RealVariable x("x");
    List<DottedRealAssignment> dynamics;
    dynamics.append(dot(x) = -x);
    OrdinaryDifferentialAlgebraicEquation decay(dynamics);
    ContinuousTimeEvolver<ApproximateTag> simulator(decay);
    simulator.configuration().set_step_size(0.05);
    simulator.configuration().set_num_subdivisions(0);
    simulator.configuration().insert_subspace(x, 1);

    RealExpressionBoundedConstraintSet initial({0.9_dec <= x <= 1.1_dec});
    auto simulated = simulator.orbit(initial, 0.1_dec);
    const auto ncurves = simulated.curves().size();
    if (ncurves != 1u) {
        throw std::runtime_error("Expected one short continuous trajectory");
    }
    std::cout << "Continuous evolution: " << ncurves << " trajectory\n";

    // Validated continuous evolution: an enclosure, not a sampled trajectory.
    // The short time horizon keeps the demonstration inexpensive.
    TaylorPicardIntegrator integrator(StepMaximumError(1e-5));
    ContinuousTimeEvolver<ValidatedTag> validated(decay, integrator);
    validated.configuration().set_maximum_step_size(0.05);
    validated.configuration().set_maximum_enclosure_radius(1.0);
    auto validated_orbit = validated.orbit(initial, Real(0.05_dec), Semantics::UPPER);
    std::cout << "ContinuousTimeEvolver<ValidatedTag>: " << validated_orbit.final().size()
              << " final enclosure(s)\n";

    // The same evolver and Orbit type handle measurable uncertain inputs.
    // Exactly one tube-rigorous reach method and one endpoint method are chosen.
    RealVariable u("u");
    OrdinaryDifferentialAlgebraicInclusion disturbed(
        {dot(x) = u - x}, RealVariablesBox({-1 <= u <= 1}));
    ContinuousTimeEvolver<ValidatedTag> inclusion(disturbed, integrator, ZeroApproximation(), ConstantApproximation());
    inclusion.configuration().set_maximum_step_size(0.03125);
    const auto inclusion_orbit = inclusion.orbit(initial, Real(0.0625_x), Semantics::UPPER);
    std::cout << "Continuous inclusion: " << inclusion_orbit.final().size() << " final enclosure(s)\n";

    // Discrete-time system: one step of a contracting difference equation.
    List<PrimedRealAssignment> updates;
    updates.append(next(x) = 0.5_dec * x);
    DifferenceAlgebraicEquation contraction(updates);
    DiscreteTimeEvolver<ValidatedTag> evolver(contraction);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    auto discrete = evolver.orbit(initial, DifferenceAlgebraicEquation::TimeType(1), Semantics::UPPER);
    std::cout << "Discrete evolution: " << discrete.final().size()
              << " final enclosure(s)\n";

    // Discrete-only hybrid system: the same two-mode automaton as the
    // Python tutorial. The simulation budget counts time and jumps.
    const DiscreteLocation idle(0), active(1);
    const DiscreteEvent activate("activate"), deactivate("deactivate");
    ContinuousTimeHybridAutomaton automaton;
    automaton.new_mode(idle);
    automaton.new_mode(active);
    automaton.new_transition(idle, activate, active);
    automaton.new_transition(active, deactivate, idle);

    if (!automaton.has_transition(idle, activate) ||
        automaton.target(active, deactivate) != idle) {
        throw std::runtime_error("Hybrid transition relation is incorrect");
    }
    HybridContinuousTimeEvolver<ApproximateTag> hybrid_continuous_time_simulator(automaton);
    HybridTime budget(1.0_x, 2);
    if (budget.discrete_time() != Integer(2)) {
        throw std::runtime_error("Incorrect hybrid jump budget");
    }
    std::cout << "Hybrid evolution: " << automaton.locations().size()
              << " modes, budget " << budget << "\n";
    (void)hybrid_continuous_time_simulator;

    // Hybrid validated evolution with continuous dynamics in a single mode.
    // No events are enabled: this isolates the continuous flow in a hybrid
    // enclosure, with a short maximum time and at most one discrete jump.
    ContinuousTimeHybridAutomaton flow_automaton;
    List<DottedRealAssignment> hybrid_dynamics;
    hybrid_dynamics.append(dot(x) = -x);
    flow_automaton.new_mode(idle, hybrid_dynamics);
    HybridContinuousTimeEvolver<ValidatedTag> hybrid_continuous_time_evolver(flow_automaton);
    hybrid_continuous_time_evolver.set_integrator(TaylorPicardIntegrator(StepMaximumError(1e-5)));
    hybrid_continuous_time_evolver.configuration().set_maximum_step_size(0.05);
    hybrid_continuous_time_evolver.configuration().set_maximum_enclosure_radius(1.0);
    // ExactBoxType uses double-exact bounds: use dyadic values 7/8 and 9/8,
    // not non-dyadic decimal fractions.
    HybridExactBox hybrid_initial(idle, flow_automaton.continuous_state_space(idle),
                                  ExactBoxType{
                                      {0.875_x, 1.125_x}
    });
    auto hybrid_orbit = hybrid_continuous_time_evolver.orbit(hybrid_initial, HybridTime(Real(1 / 20_q), Integer(1)), Semantics::UPPER);
    std::cout << "HybridContinuousTimeEvolver<ValidatedTag>: " << hybrid_orbit.final().size()
              << " final enclosure(s)\n";

    // Inclusion modes use the same hybrid evolver and explicit numerical choices.
    // This mode has no outgoing events; its inputs may vary measurably in time.
    ContinuousTimeHybridAutomaton inclusion_automaton;
    inclusion_automaton.new_mode(idle, disturbed);
    HybridContinuousTimeEvolver<ValidatedTag> hybrid_inclusion(inclusion_automaton);
    hybrid_inclusion.set_integrator(integrator, ZeroApproximation(), ConstantApproximation());
    hybrid_inclusion.configuration().set_maximum_step_size(0.03125);
    auto hybrid_included = hybrid_inclusion.orbit(hybrid_initial, HybridTime(0.0625_x, 1), Semantics::UPPER);
    std::cout << "Hybrid inclusion: " << hybrid_included.final().size() << " final enclosure(s)\n";

    // First-order PDE: solve a small two-dimensional constant-coefficient
    // three-component hyperbolic system with smooth initial conditions.
    const Real rho = 1 / 2_q, c = 1;
    Matrix<Real> A = {
        {rho,   0, 0},
        {  0, rho, 0},
        {  0,   0, 1}
    };
    Matrix<Real> B0 = {
        {          0, 0, 1},
        {          0, 0, 0},
        {rho * c * c, 0, 0}
    };
    Matrix<Real> B1 = {
        {0,           0, 0},
        {0,           0, 1},
        {0, rho * c * c, 0}
    };
    Array<Matrix<Real>> Bs = {B0, B1};
    Array<DiagonalMatrix<Real>> Ds = {
        DiagonalMatrix<Real>(Array<Real>{c, 0, -c}),
        DiagonalMatrix<Real>(Array<Real>{c, 0, -c})};
    Real det = sqrt(1 + sqr(rho) * sqr(c));
    Real c1 = 1 / det, cr = (rho * c) / det;
    Matrix<Real> T0 = {
        {c1, 0,  c1},
        { 0, 1,   0},
        {cr, 0, -cr}
    };
    Matrix<Real> T1 = {
        { 0, 1,   0},
        {c1, 0,  c1},
        {cr, 0, -cr}
    };
    Array<Matrix<Real>> Ts = {T0, T1};
    auto coords = EffectiveVectorMultivariateFunction::identity(EuclideanDomain(2));
    auto zero = EffectiveScalarMultivariateFunction::zero(EuclideanDomain(2));
    EffectiveVectorMultivariateFunction forcing =
        EffectiveVectorMultivariateFunction::zeros(3, EuclideanDomain(3));
    EffectiveVectorMultivariateFunction profile = {sin(coords[0]), sin(coords[1]), zero};
    FirstOrderPDE pde(A, forcing);
    pde.add_direction(B0, Ds[0], T0).add_direction(B1, Ds[1], T1);
    pde.validate();
    auto solution = first_order_pde(pde, profile, double_precision);
    std::cout << "FirstOrderPDE: " << solution.uts.size(0) << " x "
              << solution.uts.size(1) << " x " << solution.uts.size(2)
              << " grid samples\n";
    return 0;
}
