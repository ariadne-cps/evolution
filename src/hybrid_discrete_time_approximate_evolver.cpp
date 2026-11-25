/***************************************************************************
 *            src/hybrid_discrete_time_approximate_evolver.cpp
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

/** @file src/hybrid_discrete_time_approximate_evolver.cpp
 *  @brief Branching sampled iteration of native hybrid difference dynamics.
 */

#include "evolution/hybrid_discrete_time_approximate_evolver.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

#include "evolution/hybrid_discrete_time_step_planner.hpp"
#include "evolution/initial_point_sampling.hpp"
#include "function/functional.hpp"
#include "function/function.hpp"
#include "system/system_function.hpp"

namespace Ariadne {
namespace {

using Evolver = HybridDiscreteTimeEvolver<ApproximateTag>;
using ApproximateHybridPath = HybridDiscreteTimeApproximatePath;
using State = Point<FloatDPApproximation>;
using Values = Vector<FloatDPApproximation>;
using Evidence = DiscreteTimePredicateEvidence;

State canonical_state(const HybridApproximatePoint& point, const RealSpace& space) {
    State state(space.dimension(), FloatDPApproximation(0, dp));
    for (SizeType i = 0u; i != space.dimension(); ++i) {
        if (!point.space().contains(space[i]))
            throw std::invalid_argument("Approximate hybrid point is missing a state variable");
        state[i] = point.point()[point.space()[space[i]]];
    }
    return state;
}

Values arguments(const State& state, const Values& inputs) {
    Values result(state.dimension() + inputs.size(), FloatDPApproximation(0, dp));
    for (SizeType i = 0u; i != state.dimension(); ++i)
        result[i] = state[i];
    for (SizeType i = 0u; i != inputs.size(); ++i)
        result[state.dimension() + i] = inputs[i];
    return result;
}

Values sample_inputs(const DiscreteTimeHybridAutomatonInterface& system,
                     const InputSamplingPolicy* policy, DiscreteLocation location,
                     SizeType trajectory, SizeType tick) {
    if (!system.has_inclusion(location))
        return Values(0u, FloatDPApproximation(0, dp));
    if (policy == nullptr)
        throw std::invalid_argument("Approximate hybrid inclusion requires an explicit input policy");
    const auto& inclusion = system.inclusion(location);
    const FloatDP time(Integer(tick), near, dp);
    const InputSamplingContext context{trajectory, tick,
        FloatDPApproximation(ExactDouble(time.get_d()), dp), FloatDPApproximation(1, dp)};
    const Values sample = policy->sample(inclusion.inputs(), context);
    if (sample.size() != inclusion.number_of_inputs())
        throw std::invalid_argument("Hybrid approximate input dimension mismatch");
    for (SizeType i = 0u; i != sample.size(); ++i) {
        const double x = sample[i].raw().get_d();
        if (!std::isfinite(x) ||
            x < inclusion.inputs()[i].lower_bound().get_d() ||
            x > inclusion.inputs()[i].upper_bound().get_d())
            throw std::invalid_argument("Hybrid approximate input outside declared box");
    }
    return sample;
}

HybridApproximatePoint labelled_point(
    const DiscreteTimeHybridAutomatonInterface& system,
    DiscreteLocation location, const State& state, const Values& inputs) {
    const RealSpace state_space = system.state_space(location);
    const RealSpace output_space = system.has_equation(location)
        ? system.equation(location).state_auxiliary_space()
        : system.inclusion(location).state_auxiliary_space();
    const RealSpace auxiliary_space = system.has_equation(location)
        ? system.equation(location).auxiliary_space()
        : system.inclusion(location).auxiliary_space();
    const Values aux = evaluate(system.auxiliary_function(location), arguments(state, inputs));
    State output(output_space.dimension(), FloatDPApproximation(0, dp));
    for (SizeType i = 0u; i != output_space.dimension(); ++i) {
        const RealVariable variable = output_space[i];
        output[i] = state_space.contains(variable) ? state[state_space[variable]]
                                                  : aux[auxiliary_space[variable]];
    }
    return HybridApproximatePoint(location, output_space, output);
}

Evidence classify(const DiscreteTimeHybridAutomatonInterface& system,
                  DiscreteLocation location, const State& state,
                  const ContinuousPredicate& predicate) {
    const auto& auxiliaries = system.has_equation(location)
        ? system.equation(location).auxiliary_assignments()
        : system.inclusion(location).auxiliary_assignments();
    const auto fn = make_constraint_function(system.state_space(location),
                                              algebraic_sort(auxiliaries),
                                              predicate, Sign::POSITIVE);
    if (fn.argument_size() != state.dimension())
        throw std::invalid_argument("Hybrid predicate depends on unsupported free input variables");
    const double x = evaluate(fn, Values(state)).raw().get_d();
    if (!std::isfinite(x))
        throw std::invalid_argument("Hybrid predicate produced a nonfinite numerical value");
    if (x > 0.0)
        return Evidence::TRUE;
    if (x < 0.0)
        return Evidence::FALSE;
    // Numerical equality is not a proof of strict predicate polarity.
    // Both alternatives remain possible; no arbitrary event priority.
    return Evidence::UNKNOWN;
}

List<PrimedRealAssignment> ordered_reset(
    const List<PrimedRealAssignment>& assignments,
    const RealSpace& destination) {
    List<PrimedRealAssignment> result;
    for (const RealVariable& variable : destination.variables()) {
        Bool found = false;
        for (const PrimedRealAssignment& assignment : assignments) {
            if (assignment.lhs.base() == variable) {
                result.append(assignment);
                found = true;
                break;
            }
        }
        if (!found)
            throw std::invalid_argument("Hybrid approximate reset missing destination state");
    }
    return result;
}

State apply_reset(const DiscreteTimeHybridAutomatonInterface& system,
                  DiscreteLocation location, DiscreteEvent event, const State& state) {
    const auto destination = system.target(location, event);
    const auto aux_assignments = system.has_equation(location)
        ? system.equation(location).auxiliary_assignments()
        : system.inclusion(location).auxiliary_assignments();
    const auto map = make_reset_function(system.state_space(location),
        algebraic_sort(aux_assignments),
        ordered_reset(system.reset_assignments(location, event), system.state_space(destination)));
    return State(evaluate(map, Values(state)));
}

struct Frontier {
    ApproximateHybridPath path;
    State state;
    Values inputs;
};

} // namespace

HybridDiscreteTimeEvolver<ApproximateTag>::HybridDiscreteTimeEvolver(const SystemType& system)
    : _system(system.clone()), _configuration(std::make_shared<ConfigurationType>()),
      _input_sampling(nullptr) {}

HybridDiscreteTimeEvolver<ApproximateTag>::HybridDiscreteTimeEvolver(
    const SystemType& system, const InputSamplingPolicy& policy)
    : _system(system.clone()), _configuration(std::make_shared<ConfigurationType>()),
      _input_sampling(policy.clone()) {}

HybridDiscreteTimeEvolver<ApproximateTag>::HybridDiscreteTimeEvolver(
    const HybridDiscreteTimeEvolver& other)
    : _system(other._system->clone()),
      _configuration(std::make_shared<ConfigurationType>(*other._configuration)),
      _input_sampling(other._input_sampling ? other._input_sampling->clone() : nullptr) {}

auto HybridDiscreteTimeEvolver<ApproximateTag>::orbit(
    const InitialType& initial, const TerminationType& horizon) const -> OrbitType {
    InitialListType points(1u, initial);
    return orbit(points, horizon);
}

auto HybridDiscreteTimeEvolver<ApproximateTag>::orbit(
    const UpperBoxType& box, DiscreteLocation location,
    const TerminationType& horizon) const -> OrbitType {
    if (!system().has_mode(location))
        throw std::invalid_argument("Hybrid approximate initial location does not exist");
    return orbit(EvolutionDetail::sample_hybrid_initial_points(
        box, location, system().state_space(location), configuration().sampling_counts()), horizon);
}

auto HybridDiscreteTimeEvolver<ApproximateTag>::orbit(
    const InitialListType& initial, const TerminationType& horizon) const -> OrbitType {
    if (horizon < Integer(0))
        throw std::invalid_argument("Hybrid approximate horizon must be nonnegative");
    if (initial.size() > configuration().maximum_paths())
        throw std::overflow_error("Hybrid approximate initial trajectories exceed path budget");
    std::vector<Frontier> pending;
    for (SizeType i = 0u; i != initial.size(); ++i) {
        const auto location = initial[i].location();
        if (!system().has_mode(location))
            throw std::invalid_argument("Hybrid approximate initial mode does not exist");
        const State state = canonical_state(initial[i], system().state_space(location));
        const Values inputs = sample_inputs(system(), _input_sampling.get(), location, i, 0u);
        pending.push_back({ApproximateHybridPath{i, {labelled_point(system(), location, state, inputs)}, {}, false},
                           state, inputs});
    }

    OrbitType result;
    SizeType tick = 0u;
    for (Integer index = 0; index < horizon; ++index) {
        std::vector<Frontier> successors;
        for (const Frontier& entry : pending) {
            const DiscreteLocation location = entry.path.points.back().location();
            const auto evaluator = [&](DiscreteLocation, DiscreteEvent,
                                       DiscreteTimePredicateRole, const ContinuousPredicate& predicate) {
                return classify(system(), location, entry.state, predicate);
            };
            const auto plan = HybridDiscreteTimeStepPlanner::plan(system(), location, index, evaluator);
            if (plan.branches.empty()) {
                ApproximateHybridPath stopped = entry.path;
                stopped.blocked = true;
                result._paths.push_back(std::move(stopped));
                continue;
            }
            for (const auto& branch : plan.branches) {
                if (successors.size() + result._paths.size() >= configuration().maximum_paths())
                    throw std::overflow_error("Hybrid approximate branch budget exceeded");
                const State next_state = branch.kind == HybridDiscreteTimeBranchKind::TRANSITION
                    ? apply_reset(system(), location, *branch.event, entry.state)
                    : State(evaluate(system().update_function(location),
                                     arguments(entry.state, entry.inputs)));
                const SizeType new_path_id = successors.size();
                const Values next_inputs = sample_inputs(system(), _input_sampling.get(),
                    branch.target, new_path_id, tick + 1u);
                ApproximateHybridPath branch_path = entry.path;
                branch_path.points.push_back(labelled_point(system(), branch.target, next_state, next_inputs));
                branch_path.events.push_back(branch.event);
                successors.push_back({std::move(branch_path), next_state, next_inputs});
            }
        }
        pending = std::move(successors);
        ++tick;
    }
    for (auto& survivor : pending)
        result._paths.push_back(std::move(survivor.path));
    return result;
}

OutputStream& HybridDiscreteTimeEvolverConfiguration<ApproximateTag>::_write(OutputStream& os) const {
    return os << "HybridDiscreteTimeEvolverConfiguration<ApproximateTag>(sampling_counts="
              << sampling_counts() << ", maximum_paths=" << maximum_paths() << ")";
}

} // namespace Ariadne
