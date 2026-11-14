/***************************************************************************
 *            benchmark_recombination.cpp
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

/** @file benchmarks/benchmark_recombination.cpp
 *  @brief Measure the cost and geometric inflation of recombining actual hybrid jump enclosures.
 */

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "evolution/enclosure_recombiner.hpp"
#include "evolution/hybrid_evolver.hpp"
#include "evolution/hybrid_orbit.hpp"
#include "evolution/hybrid_time.hpp"
#include "integration/integrator.hpp"
#include "system/hybrid_automata.hpp"
#include "threading/thread_manager.hpp"
#include "threading/threading_command_line_interface.hpp"

#include "ariadne-benchmarking.hpp"

using namespace Ariadne;
namespace Bench = Ariadne::Benchmarking;

namespace {

struct Options {
    double step = 1.0 / 16.0;
    double accuracy = 1e-5;
    unsigned group_size = 4u;
    unsigned samples = 3u;
    bool sweep = false;
    std::vector<double> steps = {1.0 / 16.0, 1.0 / 64.0, 1.0 / 256.0};
    std::vector<unsigned> groups = {2u, 4u, 8u};
    std::string csv = "recombination-sweep.csv";
    std::optional<std::pair<double, double>> time_window;
    std::vector<double> window_widths;
    std::optional<unsigned> distributed_count;
    std::vector<unsigned> distributed_counts;
};

// Parse a bounded-size CLI axis without relying on textual benchmark output.
std::vector<double> parse_steps(const String& value) {
    const std::string input(value);
    if (input.empty() || input.front() == ',' || input.back() == ',')
        throw std::invalid_argument("expected comma-separated positive steps");
    std::vector<double> steps;
    std::size_t start = 0u;
    while (start < input.size()) {
        const auto end = input.find(',', start);
        const std::string token = input.substr(start, end == std::string::npos ? end : end - start);
        std::size_t used = 0u;
        const double step = std::stod(token, &used);
        if (used != token.size() || !std::isfinite(step) || step <= 0.0 ||
            std::find(steps.begin(), steps.end(), step) != steps.end())
            throw std::invalid_argument("invalid or repeated positive step");
        steps.push_back(step);
        if (end == std::string::npos)
            break;
        start = end + 1u;
    }
    return steps;
}


// A fixed absolute birth-time window separates group cardinality from the
// integration-step resolution. Reject malformed endpoints before any evolution.
std::pair<double, double> parse_time_window(const String& value) {
    const std::string text(value);
    const auto comma = text.find(',');
    if (comma == std::string::npos || text.find(',', comma + 1u) != std::string::npos)
        throw std::invalid_argument("expected --time-window START,END");
    const auto parse = [](const std::string& token) {
        std::size_t used = 0u;
        const double number = std::stod(token, &used);
        if (used != token.size() || !std::isfinite(number))
            throw std::invalid_argument("invalid time-window bound");
        return number;
    };
    const double start = parse(text.substr(0u, comma));
    const double end = parse(text.substr(comma + 1u));
    if (start < 0.0 || end <= start)
        throw std::invalid_argument("time-window requires 0 <= START < END");
    return {start, end};
}

List<ArgumentParser> argument_parsers(Options& options) {
    List<ArgumentParser> parsers;
    auto positive = [&](const String& short_id, const String& name, const String& description, double& value) {
        parsers.append(make_valued_argument_parser(short_id, name, description,
                                                   [&value](const String& text) -> VoidFunction {
                                                       std::size_t consumed = 0u;
                                                       double parsed = std::stod(std::string(text), &consumed);
                                                       if (consumed != std::string(text).size() || !std::isfinite(parsed) || parsed <= 0.0)
                                                           throw std::invalid_argument("expected a finite positive number");
                                                       return [&value, parsed] { value = parsed; };
                                                   }));
    };
    auto count = [&](const String& short_id, const String& name, const String& description, unsigned& value) {
        parsers.append(make_valued_argument_parser(short_id, name, description,
                                                   [&value](const String& text) -> VoidFunction {
                                                       unsigned parsed = Bench::parse_unsigned<unsigned>(text);
                                                       if (parsed == 0u)
                                                           throw std::invalid_argument("expected a positive count");
                                                       return [&value, parsed] { value = parsed; };
                                                   }));
    };
    positive("p", "step", "Maximum integration step (default: 1/16)", options.step);
    positive("a", "accuracy", "Integrator step maximum error (default: 1e-5)", options.accuracy);
    count("g", "group-size", "Nearby jump enclosures to recombine (default: 4)", options.group_size);
    count("n", "samples", "Measured runs per scenario (default: 3)", options.samples);
    parsers.append(make_valued_argument_parser("", "time-window",
                                               "Use all jump enclosures intersecting absolute birth-time window START,END (ignores group-size)",
                                               [&options](const String& value) -> VoidFunction {
                                                   const auto bounds = parse_time_window(value);
                                                   return [&options, bounds] { options.time_window = bounds; };
                                               }));
    parsers.append(make_valued_argument_parser("", "distributed-count",
                                               "Fixed number of time-distributed source enclosures (requires --time-window)",
                                               [&options](const String& value) -> VoidFunction {
                                                   const auto count = Bench::parse_unsigned<unsigned>(value);
                                                   if (count < 2u)
                                                       throw std::invalid_argument("distributed-count must be >= 2");
                                                   return [&options, count] { options.distributed_count = count; };
                                               }));
    parsers.append(make_valued_argument_parser("", "distributed-counts",
                                               "Sweep source counts within the same window, comma-separated",
                                               [&options](const String& value) -> VoidFunction {
                                                   auto counts = Bench::parse_unique_unsigned_list(value, 2u);
                                                   return [&options, counts = std::move(counts)] {
                                                       options.distributed_counts = counts;
                                                   };
                                               }));
    parsers.append(make_valued_argument_parser("", "window-widths",
                                               "Sweep comma-separated widths centred on --time-window midpoint",
                                               [&options](const String& value) -> VoidFunction {
                                                   auto widths = parse_steps(value);
                                                   return [&options, widths = std::move(widths)] {
                                                       options.window_widths = widths;
                                                   };
                                               }));
    parsers.append(make_valued_argument_parser("", "sweep",
                                               "Run the built-in C++ step/group sweep: 1 enables (default: 0)",
                                               [&options](const String& value) -> VoidFunction {
                                                   const auto enabled = Bench::parse_unsigned<unsigned>(value);
                                                   if (enabled > 1u)
                                                       throw std::invalid_argument("sweep must be 0 or 1");
                                                   return [&options, enabled] { options.sweep = enabled == 1u; };
                                               }));
    parsers.append(make_valued_argument_parser("", "steps",
                                               "Sweep maximum steps, comma-separated (default: 1/16,1/64,1/256)",
                                               [&options](const String& value) -> VoidFunction {
                                                   auto steps = parse_steps(value);
                                                   return [&options, steps = std::move(steps)] { options.steps = steps; };
                                               }));
    parsers.append(make_valued_argument_parser("", "groups",
                                               "Sweep group sizes, comma-separated (default: 2,4,8)",
                                               [&options](const String& value) -> VoidFunction {
                                                   auto groups = Bench::parse_unique_unsigned_list(value, 2u);
                                                   return [&options, groups = std::move(groups)] { options.groups = groups; };
                                               }));
    parsers.append(make_valued_argument_parser("", "csv",
                                               "CSV filename for the C++ sweep (default: recombination-sweep.csv)",
                                               [&options](const String& value) -> VoidFunction {
                                                   const std::string filename(value);
                                                   if (filename.empty())
                                                       throw std::invalid_argument("CSV filename cannot be empty");
                                                   return [&options, filename] { options.csv = filename; };
                                               }));
    return parsers;
}

struct Scenario {
    const char* name;
    EventKind kind;
    double generation_horizon;
    double continuation_horizon;
};

HybridAutomaton make_system(const Scenario& scenario, DiscreteLocation source,
                            DiscreteLocation target, DiscreteEvent jump,
                            RealVariable x, RealVariable y) {
    HybridAutomaton system(scenario.name);
    system.new_mode(source, {dot(x) = 1, dot(y) = sqr(x)});
    system.new_mode(target, {dot(x) = 1 + sqr(y), dot(y) = 1});
    system.new_transition(source, jump, target, {next(x) = x, next(y) = y},
                          x + y >= 0, scenario.kind);
    return system;
}

void configure(GeneralHybridEvolver& evolver, const Options& options) {
    evolver.set_integrator(TaylorPicardIntegrator(StepMaximumError(options.accuracy)));
    evolver.configuration().set_maximum_step_size(ApproximateDouble(options.step));
    // Prevent radius-driven splitting from dominating the first baseline.
    evolver.configuration().set_maximum_enclosure_radius(ApproximateDouble(100.0));
}

List<Enclosure> generate(const Scenario& scenario, const Options& options,
                               DiscreteLocation source, DiscreteLocation target,
                               DiscreteEvent jump, RealVariable x, RealVariable y,
                               bool verbose, SizeType& generated, SizeType& eligible) {
    GeneralHybridEvolver evolver(make_system(scenario, source, target, jump, x, y));
    configure(evolver, options);
    const HybridBoxSet initial = scenario.kind == EventKind::URGENT
        ? HybridBoxSet(source, {Real(-2) <= x <= Real(-1), Real(-1 / 2_q) <= y <= Real(1 / 2_q)})
        : HybridBoxSet(source, {Real(1) <= x <= Real(2), Real(-1 / 2_q) <= y <= Real(1 / 2_q)});
    const HybridTerminationCriterion stop(Real(ExactDouble(scenario.generation_horizon)), Integer(8),
                                           Set<DiscreteEvent>{jump});
    const auto orbit = evolver.orbit(initial, stop, Semantics::UPPER);
    if (verbose)
        std::cout << "  generation reach=" << orbit.reach().size()
                  << " intermediate=" << orbit.intermediate().size()
                  << " final=" << orbit.final().size()
                  << " (not accepted steps)\n";
    List<HybridEnclosure> candidates;
    for (const auto& enclosure: orbit.final()) {
        if (enclosure.location() == target && enclosure.previous_events().size() == 1u &&
            enclosure.previous_events()[0] == jump)
            candidates.append(enclosure);
    }
    std::sort(candidates.begin(), candidates.end(), [](const HybridEnclosure& a, const HybridEnclosure& b) {
        return a.time_range().lower_bound().raw().get_d() < b.time_range().lower_bound().raw().get_d();
    });
    generated = candidates.size();
    if (candidates.size() < 2u)
        throw std::runtime_error("not enough post-reset checkpoints; adjust integration step");
    List<Enclosure> selection;
    eligible = 0u;
    if (options.time_window) {
        const auto [window_start, window_end] = *options.time_window;
        if (window_end >= scenario.generation_horizon)
            throw std::invalid_argument("time-window END must be before generation horizon");
        const auto absolute_time = EffectiveScalarMultivariateFunction::coordinate(3u, 2u);
        List<Enclosure> window_sources;
        for (const auto& checkpoint: candidates) {
            const auto bounds = checkpoint.time_range();
            if (bounds.upper_bound().raw().get_d() < window_start ||
                bounds.lower_bound().raw().get_d() > window_end)
                continue;
            auto source = checkpoint.continuous_set();
            source.new_state_time_constraint(ValidatedConstraint(
                ExactDouble(window_start), absolute_time, ExactDouble(window_end)));
            auto probe = source;
            if (!definitely(probe.is_empty()))
                window_sources.append(std::move(source));
        }
        eligible = window_sources.size();
        if (!options.distributed_count) {
            selection = std::move(window_sources);
        } else {
            const SizeType count = *options.distributed_count;
            if (count > eligible)
                throw std::runtime_error("distributed-count exceeds available window sources");
            // Evenly distributed deterministic indices, including the first
            // and last source. Only a subset is tested, not the complete union.
            for (SizeType i = 0; i != count; ++i) {
                const SizeType index = i * (eligible - 1u) / (count - 1u);
                selection.append(window_sources[index]);
            }
        }
    } else {
        const SizeType take = std::min<SizeType>(candidates.size(), options.group_size);
        if (take < 2u)
            throw std::runtime_error("group-size must be at least two");
        const SizeType offset = (candidates.size() - take) / 2u;
        for (SizeType i = offset; i != offset + take; ++i)
            selection.append(candidates[i].continuous_set());
        eligible = candidates.size();
    }
    if (selection.size() < 2u)
        throw std::runtime_error("fewer than two nonempty post-reset checkpoints in selection");
    if (verbose) {
        std::cout << scenario.name << ": generated=" << candidates.size()
                  << " post-reset checkpoints, eligible=" << eligible
                  << " selected=" << selection.size();
        if (options.time_window) {
            std::cout << " with constrained time-window=[" << options.time_window->first
                      << "," << options.time_window->second << "]";
            std::cout << (options.distributed_count
                ? " (distributed subset; not full-window coverage)"
                : " (all eligible sources)");
        } else {
            std::cout << " (middle adjacent time window)";
        }
        std::cout << "\n";
    }
    return selection;
}

// Outward-rounded projections: wrapping diagnostics, not set-inclusion proofs.
struct Projection {
    UpperBoxType state;
    UpperIntervalType time, dwell, invariant;
};

void accumulate(std::optional<Projection>& total, const Enclosure& enclosure) {
    const auto domain = enclosure.reduced_domain();
    const auto x = enclosure.state_function().get(0);
    const auto y = enclosure.state_function().get(1);
    // For the target flow x'=1+y^2, y'=1, the quantity
    // I=x-y-y^3/3 is constant along exact trajectories.
    // Evaluation on a function patch (rather than on a state box) retains
    // some nonlinear correlations; it remains an outer bound, not a proof.
    const auto invariant = x - y - FloatDPBounds(Real(1 / 3_q), dp) * y * sqr(y);
    Projection current{apply(enclosure.state_function(), domain),
                       apply(enclosure.time_function(), domain),
                       apply(enclosure.dwell_time_function(), domain),
                       apply(invariant, domain)};
    if (!total) {
        total = std::move(current);
    } else {
        total->state = hull(total->state, current.state);
        total->time = hull(total->time, current.time);
        total->dwell = hull(total->dwell, current.dwell);
        total->invariant = hull(total->invariant, current.invariant);
    }
}

// Closed-form target flow independent of the numerical integrator:
// y(T)=q+(T-tau), x(T)=p+(1+q^2)(T-tau)+q(T-tau)^2+(T-tau)^3/3.
// Arithmetic is validated over each input's complete reduced parameter box;
// no source parameters are sampled or replaced by a manually constructed box.
Projection analytic_target_projection(const Enclosure& source, double final_time) {
    const auto domain = source.reduced_domain();
    const auto p = source.state_function().get(0);
    const auto q = source.state_function().get(1);
    const auto delta = FloatDPBounds(Real(ExactDouble(final_time)), dp) - source.time_function();
    const auto delta2 = sqr(delta);
    const auto one_third = FloatDPBounds(Real(1 / 3_q), dp);
    const auto x = p + (1 + sqr(q)) * delta + q * delta2 + one_third * delta2 * delta;
    const auto y = q + delta;
    const auto t = source.time_function() + delta;
    const auto dwell = source.dwell_time_function() + delta;
    const auto invariant = x - y - one_third * y * sqr(y);
    return {UpperBoxType({apply(x, domain), apply(y, domain)}),
            apply(t, domain), apply(dwell, domain), apply(invariant, domain)};
}

// Recover the source initial coordinates from the post-reset model using
// the *independent* analytic source flow. Interval intersection tests can
// refute gross inconsistencies; they do not certify that all represented
// points obey the source initial-set constraints.
void check_analytic_source_consistency(const Enclosure& source, bool urgent) {
    const auto domain = source.reduced_domain();
    const auto x = source.state_function().get(0);
    const auto y = source.state_function().get(1);
    const auto tau = source.time_function();
    const auto a = x - tau;
    const auto b = y - FloatDPBounds(Real(1 / 3_q), dp) * (x * sqr(x) - a * sqr(a));
    const auto initial_a = urgent ? UpperIntervalType(-2, -1) : UpperIntervalType(1, 2);
    const UpperIntervalType initial_b(-0.5_x, 0.5_x);
    if (definitely(intersection(apply(a, domain), initial_a).is_empty()) ||
        definitely(intersection(apply(b, domain), initial_b).is_empty()))
        throw std::runtime_error("analytical source-flow reconstruction contradicts initial set");
    const auto guard = apply(x + y, domain);
    if (urgent && definitely(intersection(guard, UpperIntervalType(0, 0)).is_empty()))
        throw std::runtime_error("urgent post-reset checkpoint contradicts guard equality");
    if (!urgent && guard.upper_bound().raw().get_d() < 0.0)
        throw std::runtime_error("permissive post-reset checkpoint contradicts guard");
}

void check_reference_overlap(const Projection& reference, const Projection& evolved,
                             const char* path) {
    // A definite disjointness of analytic and numerical *outer* ranges is a
    // necessary-condition failure. Mere overlap is NOT inclusion certification.
    for (SizeType k = 0; k != reference.state.size(); ++k)
        if (definitely(intersection(reference.state[k], evolved.state[k]).is_empty()))
            throw std::runtime_error(std::string(path) + ": analytic and numerical state projections are disjoint");
    if (definitely(intersection(reference.time, evolved.time).is_empty()) ||
        definitely(intersection(reference.dwell, evolved.dwell).is_empty()))
        throw std::runtime_error(std::string(path) + ": analytic and numerical clock projections are disjoint");
}

struct Measurement {
    double contraction_ns = 0, recombination_ns = 0;
    double independent_ns = 0, merged_ns = 0;
    SizeType independent_final = 0, merged_final = 0;
    EnclosureRecombinationResult result;
    std::optional<Projection> source, merged, continued_independent, continued_merged;
    std::optional<Projection> analytical_final;
    std::vector<Projection> analytical_per_source;
    double independent_total_ns() const { return contraction_ns + independent_ns; }
    double merged_total_ns() const { return contraction_ns + recombination_ns + merged_ns; }
};

Measurement measure(const List<Enclosure>& original,
                    const GeneralHybridEvolver& continuation,
                    const HybridTerminationCriterion& stop,
                    DiscreteLocation target, RealSpace space,
                    double final_time, bool urgent) {
    Measurement m;
    List<Enclosure> contracted = original;
    const auto contraction_start = Bench::Clock::now();
    for (const auto& enclosure: contracted)
        enclosure.reduce();
    const auto contraction_end = Bench::Clock::now();
    m.contraction_ns = Bench::elapsed_ns(contraction_start, contraction_end);

    const auto recombination_start = Bench::Clock::now();
    m.result = ReferenceEnclosureRecombiner().recombine(contracted);
    const auto recombination_end = Bench::Clock::now();
    m.recombination_ns = Bench::elapsed_ns(recombination_start, recombination_end);
    if (m.result.outputs.size() != 1u || !m.result.discarded_sources.empty() ||
        m.result.outputs.front().source_indices.size() != contracted.size())
        throw std::runtime_error("unexpected recombination output or provenance");
    for (SizeType i = 0; i != contracted.size(); ++i)
        if (m.result.outputs.front().source_indices[i] != i)
            throw std::runtime_error("incorrect recombination source index");

    // Use the SAME reduced continuous input in both branches. Both wrappers
    // start in the target-only automaton with empty event histories.
    std::vector<HybridEnclosure> independent_inputs;
    independent_inputs.reserve(contracted.size());
    for (const auto& source: contracted)
        independent_inputs.emplace_back(target, space, source);
    const HybridEnclosure merged_input(target, space, m.result.outputs.front().enclosure);

    // Keep orbits alive until BOTH timers end: no projection/validation in
    // a timed region. Independent timings sum only the orbit() call costs.
    std::vector<Orbit<HybridEnclosure>> independent_orbits;
    independent_orbits.reserve(independent_inputs.size());
    for (const auto& source: independent_inputs) {
        const auto start = Bench::Clock::now();
        auto orbit = continuation.orbit(source, stop, Semantics::UPPER);
        const auto finish = Bench::Clock::now();
        m.independent_ns += Bench::elapsed_ns(start, finish);
        independent_orbits.push_back(std::move(orbit));
    }
    const auto merged_start = Bench::Clock::now();
    const auto merged_orbit = continuation.orbit(merged_input, stop, Semantics::UPPER);
    const auto merged_end = Bench::Clock::now();
    m.merged_ns = Bench::elapsed_ns(merged_start, merged_end);

    for (const auto& source: contracted) {
        check_analytic_source_consistency(source, urgent);
        accumulate(m.source, source);
        // Independently evaluate the full exact target polynomial flow,
        // retaining parameter correlations instead of evolving interval boxes.
        const auto projection = analytic_target_projection(source, final_time);
        m.analytical_per_source.push_back(projection);
        if (!m.analytical_final)
            m.analytical_final = projection;
        else {
            m.analytical_final->state = hull(m.analytical_final->state, projection.state);
            m.analytical_final->time = hull(m.analytical_final->time, projection.time);
            m.analytical_final->dwell = hull(m.analytical_final->dwell, projection.dwell);
            m.analytical_final->invariant = hull(m.analytical_final->invariant, projection.invariant);
        }
    }
    accumulate(m.merged, m.result.outputs.front().enclosure);
    for (const auto& orbit: independent_orbits) {
        m.independent_final += orbit.final().size();
        for (const auto& enclosure: orbit.final())
            accumulate(m.continued_independent, enclosure.continuous_set());
    }
    m.merged_final = merged_orbit.final().size();
    for (const auto& enclosure: merged_orbit.final())
        accumulate(m.continued_merged, enclosure.continuous_set());
    if (!m.source || !m.merged || !m.continued_independent ||
        !m.continued_merged || !m.analytical_final)
        throw std::runtime_error("missing source or final projection");
    check_reference_overlap(*m.analytical_final, *m.continued_independent, "independent union");
    check_reference_overlap(*m.analytical_final, *m.continued_merged, "recombined union");
    // Whole-union overlap can mask a source that is disjoint from the output.
    for (SizeType i = 0; i != m.analytical_per_source.size(); ++i) {
        const auto label = std::string("source ") + std::to_string(i);
        check_reference_overlap(m.analytical_per_source[i], *m.continued_independent,
                                (label + " independent").c_str());
        check_reference_overlap(m.analytical_per_source[i], *m.continued_merged,
                                (label + " recombined").c_str());
    }
    return m;
}

double width(const UpperIntervalType& value) {
    return value.upper_bound().raw().get_d() - value.lower_bound().raw().get_d();
}

// Report outward-rounded widths together with endpoint discrepancies. A
// negative hull-containment verdict cannot decide constrained-set coverage.
// Tiny time/dwell widths near a fixed horizon are roundoff, not real spread.
void compare_interval(const std::string& name, const UpperIntervalType& baseline,
                      const UpperIntervalType& candidate) {
    const double lo = baseline.lower_bound().raw().get_d();
    const double hi = baseline.upper_bound().raw().get_d();
    const double mlo = candidate.lower_bound().raw().get_d();
    const double mhi = candidate.upper_bound().raw().get_d();
    const double before = width(baseline), after = width(candidate);
    const double scale = std::max({1.0, std::abs(lo), std::abs(hi)});
    const double negligible = 32.0 * std::numeric_limits<double>::epsilon() * scale;

    std::cout << "    " << name << " independent-width=" << std::setprecision(9) << before
              << " merged-width=" << after << " ratio=";
    if (before > negligible)
        std::cout << after / before;
    else
        std::cout << "n/a (degenerate baseline)";
    // Endpoint shortfalls are comparisons of validated interval *bounds*
    // represented as doubles. They are not errors or counterexamples to set
    // inclusion, because the input interval hull may contain infeasible points.
    std::cout << "\n      independent=[" << std::setprecision(17) << lo << "," << hi
              << "] merged=[" << mlo << "," << mhi
              << "] hull-endpoint-shortfall lower=" << std::max(0.0, mlo - lo)
              << " upper=" << std::max(0.0, hi - mhi) << "\n";
}

void compare_quality(const char* phase, const Projection& independent,
                     const Projection& merged) {
    std::cout << "  " << phase << " (validated coordinate/invariant outer bounds):\n";
    for (SizeType k = 0; k != independent.state.size(); ++k)
        compare_interval("state[" + std::to_string(k) + "]",
                         independent.state[k], merged.state[k]);
    compare_interval("time", independent.time, merged.time);
    compare_interval("dwell", independent.dwell, merged.dwell);
    compare_interval("invariant(x-y-y^3/3)", independent.invariant, merged.invariant);
    std::cout << "    These bounds do not certify or refute constrained-set inclusion.\n";
}

void describe_quality(const Measurement& m) {
    const auto& output = m.result.outputs.front();
    std::cout << std::setprecision(6)
              << "  provenance=" << output.source_indices.size()
              << " discarded=" << m.result.discarded_sources.size()
              << " parameters=" << output.enclosure.number_of_parameters()
              << " constraints=" << output.enclosure.number_of_constraints()
              << " radius-limit-satisfied=" << m.result.radius_limit_satisfied << "\n";
    compare_quality("at recombination", *m.source, *m.merged);
    compare_quality("after continuation", *m.continued_independent, *m.continued_merged);
    compare_quality("analytic target flow vs independent continuation",
                    *m.analytical_final, *m.continued_independent);
    compare_quality("analytic target flow vs recombined continuation",
                    *m.analytical_final, *m.continued_merged);
    std::cout << "  Analytic overlap/reconstruction checks passed; "
                 "these are necessary conditions, NOT set-inclusion certificates.\n";
}

struct ScenarioSummary {
    SizeType generated = 0, eligible = 0, selected = 0;
    SizeType independent_final = 0, recombined_final = 0;
    double contraction_mean_ns = 0.0, recombination_mean_ns = 0.0;
    double contract_recombine_mean_ns = 0.0;
    double continue_independent_mean_ns = 0.0, continue_recombined_mean_ns = 0.0;
    double independent_mean_ns = 0.0, recombined_mean_ns = 0.0;
    double x_final_ratio = 0.0, y_final_ratio = 0.0;
    double invariant_initial_ratio = 0.0, invariant_final_ratio = 0.0;
    double analytic_x_ratio = 0.0, analytic_y_ratio = 0.0, analytic_invariant_ratio = 0.0;
};

double width_ratio(const UpperIntervalType& baseline, const UpperIntervalType& candidate) {
    const double denominator = width(baseline);
    return denominator > 0.0 ? width(candidate) / denominator
                             : std::numeric_limits<double>::quiet_NaN();
}

// This is the single numerical execution path for interactive and sweep runs.
// Sweeps receive typed data directly; no log scraping, subprocess or scripting.
ScenarioSummary run_scenario(const Scenario& scenario, const Options& options,
                             const Bench::CommandLine& cli, bool verbose) {
    const DiscreteLocation source(1), target(2);
    const DiscreteEvent jump("jump");
    const RealVariable x("x"), y("y");
    if (verbose) {
        std::cout << "\nscenario=" << scenario.name << " step=" << options.step
                  << " accuracy=" << options.accuracy;
        if (options.time_window) {
            std::cout << " time-window=[" << options.time_window->first
                      << "," << options.time_window->second << "]";
            if (options.distributed_count)
                std::cout << " distributed-count=" << *options.distributed_count;
        } else {
            std::cout << " group-size=" << options.group_size;
        }
        std::cout << "\n";
    }
    ScenarioSummary summary;
    const auto original = generate(scenario, options, source, target, jump, x, y,
                                   verbose, summary.generated, summary.eligible);
    summary.selected = original.size();

    HybridAutomaton target_system("recombination-continuation");
    target_system.new_mode(target, {dot(x) = 1 + sqr(y), dot(y) = 1});
    GeneralHybridEvolver continuation(target_system);
    configure(continuation, options);
    const HybridTerminationCriterion stop(Real(ExactDouble(scenario.continuation_horizon)), Integer(8));

    Bench::Runner runner(cli);
    Bench::Samples contractions, recombinations, totals, independents, merged;
    Bench::Samples independent_totals, merged_totals;
    std::optional<Measurement> last;
    unsigned round = 0u;
    const auto completed = runner.run(scenario.name, options.samples, [] {}, [&] {
        last = measure(original, continuation, stop, target, RealSpace({x, y}),
                       scenario.continuation_horizon, scenario.kind == EventKind::URGENT);
    }, [&] {
        if (!last)
            throw std::runtime_error("missing measurement");
        if (round >= cli.warmups()) {
            contractions.add(last->contraction_ns);
            recombinations.add(last->recombination_ns);
            totals.add(last->contraction_ns + last->recombination_ns);
            independents.add(last->independent_ns);
            merged.add(last->merged_ns);
            independent_totals.add(last->independent_total_ns());
            merged_totals.add(last->merged_total_ns());
            if (verbose)
                std::cout << "  sample=" << round - cli.warmups() + 1
                          << " independent-final=" << last->independent_final
                          << " recombined-final=" << last->merged_final
                          << " net-speedup=" << last->independent_total_ns() / last->merged_total_ns() << "\n";
        }
        ++round;
    });
    if (completed.size() != options.samples ||
        independent_totals.size() != options.samples ||
        merged_totals.size() != options.samples)
        throw std::runtime_error("incorrect measured sample count");
    if (!last)
        throw std::runtime_error("missing final sample");

    summary.contraction_mean_ns = Bench::summarize(contractions).mean;
    summary.recombination_mean_ns = Bench::summarize(recombinations).mean;
    summary.contract_recombine_mean_ns = Bench::summarize(totals).mean;
    summary.continue_independent_mean_ns = Bench::summarize(independents).mean;
    summary.continue_recombined_mean_ns = Bench::summarize(merged).mean;
    summary.independent_mean_ns = Bench::summarize(independent_totals).mean;
    summary.recombined_mean_ns = Bench::summarize(merged_totals).mean;
    summary.independent_final = last->independent_final;
    summary.recombined_final = last->merged_final;
    summary.x_final_ratio = width_ratio(last->continued_independent->state[0],
                                        last->continued_merged->state[0]);
    summary.y_final_ratio = width_ratio(last->continued_independent->state[1],
                                        last->continued_merged->state[1]);
    summary.invariant_initial_ratio = width_ratio(last->source->invariant,
                                                  last->merged->invariant);
    summary.invariant_final_ratio = width_ratio(last->continued_independent->invariant,
                                                last->continued_merged->invariant);
    summary.analytic_x_ratio = width_ratio(last->analytical_final->state[0],
                                           last->continued_merged->state[0]);
    summary.analytic_y_ratio = width_ratio(last->analytical_final->state[1],
                                           last->continued_merged->state[1]);
    summary.analytic_invariant_ratio = width_ratio(last->analytical_final->invariant,
                                                   last->continued_merged->invariant);
    if (!verbose)
        return summary;

    describe_quality(*last);
    Bench::Report report(std::string("Recombination timings: ") + scenario.name + " (nanoseconds)",
                         cli.summary_fields());
    report.set_parameter_header("One sample = fresh constraint reduction, reference recombination, "
                                "independent and merged target-only continuation. "
                                "Dataset generation, copying, quality evaluation and reporting are untimed. "
                                "Totals: contraction+independent vs contraction+recombination+merged.");
    report.set_dimension_columns({{"Scenario"}, {"Phase"}, {"Inputs", 0u, Bench::TableAlignment::Right}});
    report.set_table_options(Bench::TableOptions::ruled());
    const std::string count = std::to_string(original.size());
    const auto row = [&](const char* phase, const Bench::Samples& samples) {
        report.add_row({scenario.name, phase, count}, Bench::summarize(samples));
    };
    row("contract", contractions);
    row("recombine", recombinations);
    row("contract+recombine", totals);
    row("continue-independent", independents);
    row("continue-recombined", merged);
    row("total-independent", independent_totals);
    row("total-recombined", merged_totals);
    report.print();
    return summary;
}

std::string csv_escape(const std::string& value) {
    std::string escaped = "\"";
    for (char character: value) {
        if (character == '"')
            escaped += "\"\"";
        else if (character == '\r' || character == '\n')
            escaped += ' ';
        else
            escaped += character;
    }
    escaped += '"';
    return escaped;
}

void csv_number(std::ostream& stream, double value) {
    if (std::isfinite(value))
        stream << value;
}

// Both success and failure rows use the same typed configuration prefix.
// In time-window mode group_size is deliberately blank, not spuriously zero.
void csv_configuration(std::ostream& stream, const Options& options) {
    stream << options.step << ',' << options.accuracy << ',';
    if (options.time_window) {
        stream << (options.distributed_count ? "window-distributed" : "time-window") << ',';
        if (options.distributed_count)
            stream << *options.distributed_count;
        stream << ',' << options.time_window->first << ','
               << options.time_window->second << ','
               << options.time_window->second - options.time_window->first << ',';
    } else {
        stream << "adjacent," << options.group_size << ",,,,";
    }
}

void csv_success(std::ostream& stream, const Options& options,
                 const Scenario& scenario, const ScenarioSummary& result) {
    csv_configuration(stream, options);
    stream << scenario.name << ",ok," << result.generated << ','
           << result.eligible << ',' << result.selected << ',' << options.samples << ',' << result.independent_final << ','
           << result.recombined_final << ',';
    csv_number(stream, result.contraction_mean_ns);
    stream << ',';
    csv_number(stream, result.recombination_mean_ns);
    stream << ',';
    csv_number(stream, result.contract_recombine_mean_ns);
    stream << ',';
    csv_number(stream, result.continue_independent_mean_ns);
    stream << ',';
    csv_number(stream, result.continue_recombined_mean_ns);
    stream << ',';
    csv_number(stream, result.independent_mean_ns);
    stream << ',';
    csv_number(stream, result.recombined_mean_ns);
    stream << ',';
    csv_number(stream, result.independent_mean_ns / result.recombined_mean_ns);
    stream << ',';
    csv_number(stream, result.x_final_ratio);
    stream << ',';
    csv_number(stream, result.y_final_ratio);
    stream << ',';
    csv_number(stream, result.invariant_initial_ratio);
    stream << ',';
    csv_number(stream, result.invariant_final_ratio);
    stream << ',';
    csv_number(stream, result.analytic_x_ratio);
    stream << ',';
    csv_number(stream, result.analytic_y_ratio);
    stream << ',';
    csv_number(stream, result.analytic_invariant_ratio);
    stream << ",\"\"\n";
}

void csv_failure(std::ostream& stream, const Options& options,
                 const Scenario& scenario, const std::string& message) {
    csv_configuration(stream, options);
    stream << scenario.name << ",failed,";
    // Six count fields, five phase means, two totals, one speedup
    // and seven quality ratios.
    for (unsigned index = 0u; index != 21u; ++index)
        stream << ',';
    stream << csv_escape(message) << '\n';
}

int run_sweep(const Scenario (&scenarios)[2], const Options& options,
              const Bench::CommandLine& cli) {
    if (!options.window_widths.empty() && !options.time_window)
        throw std::invalid_argument("--window-widths requires --time-window");
    if ((!options.distributed_counts.empty() || options.distributed_count) && !options.time_window)
        throw std::invalid_argument("--distributed-count(s) requires --time-window");
    if (!options.distributed_counts.empty() && options.distributed_count)
        throw std::invalid_argument("choose --distributed-count or --distributed-counts");
    std::vector<std::optional<std::pair<double, double>>> windows;
    if (!options.time_window) {
        windows.push_back(std::nullopt);
    } else if (options.window_widths.empty()) {
        windows.push_back(options.time_window);
    } else {
        const double center = options.time_window->first / 2.0 +
                              options.time_window->second / 2.0;
        for (const auto width: options.window_widths) {
            const auto lo = center - width / 2.0;
            const auto hi = center + width / 2.0;
            if (lo < 0.0 || hi >= 1.0)
                throw std::invalid_argument("width exceeds common generation horizon");
            windows.push_back(std::make_pair(lo, hi));
        }
    }

    const std::filesystem::path path(options.csv);
    if (path.has_parent_path())
        std::filesystem::create_directories(path.parent_path());
    std::ofstream csv(path, std::ios::out | std::ios::trunc);
    if (!csv)
        throw std::runtime_error("cannot open CSV file: " + options.csv);
    csv << std::setprecision(17);
    csv << "step,accuracy,selection_mode,group_size,window_start,window_end,window_width,"
           "scenario,status,generated,eligible,selected,samples,independent_final,recombined_final,"
           "contraction_mean_ns,recombination_mean_ns,contract_recombine_mean_ns,"
           "continue_independent_mean_ns,continue_recombined_mean_ns,"
           "independent_mean_ns,recombined_mean_ns,net_speedup,"
           "x_final_ratio,y_final_ratio,invariant_initial_ratio,invariant_final_ratio,"
           "analytic_x_ratio,analytic_y_ratio,analytic_invariant_ratio,error\n";
    bool failed = false;
    for (const auto step: options.steps) {
        for (const auto& window: windows) {
            const auto groups = !window ? options.groups
                : !options.distributed_counts.empty() ? options.distributed_counts
                : std::vector<unsigned>{options.distributed_count.value_or(options.group_size)};
            for (const auto size: groups) {
                auto configuration = options;
                configuration.step = step;
                configuration.group_size = size;
                configuration.time_window = window;
                if (window && !options.distributed_counts.empty())
                    configuration.distributed_count = size;
                for (const auto& scenario: scenarios) {
                    std::cout << "RUN " << scenario.name << " step=" << step;
                    if (window) {
                        std::cout << " window=[" << window->first << "," << window->second << "]";
                        if (configuration.distributed_count)
                            std::cout << " distributed-count=" << *configuration.distributed_count;
                        else
                            std::cout << " all-sources";
                    } else {
                        std::cout << " group-size=" << size;
                    }
                    std::cout << " accuracy=" << options.accuracy << std::endl;
                    try {
                        const auto result = run_scenario(scenario, configuration, cli, false);
                        csv_success(csv, configuration, scenario, result);
                        std::cout << "  ok generated=" << result.generated
                                  << " eligible=" << result.eligible
                                  << " selected=" << result.selected
                                  << " contract-ms=" << result.contraction_mean_ns / 1.0e6
                                  << " recombine-ms=" << result.recombination_mean_ns / 1.0e6
                                  << " independent-ms=" << result.continue_independent_mean_ns / 1.0e6
                                  << " merged-ms=" << result.continue_recombined_mean_ns / 1.0e6
                                  << " speedup=" << result.independent_mean_ns / result.recombined_mean_ns
                                  << " invariant-ratio=" << result.invariant_final_ratio << '\n';
                    } catch (const std::exception& error) {
                        failed = true;
                        csv_failure(csv, configuration, scenario, error.what());
                        std::cerr << "  FAILED " << scenario.name << " step=" << step
                                  << ": " << error.what() << '\n';
                    }
                    csv.flush();
                    if (!csv)
                        throw std::runtime_error("cannot write CSV file: " + options.csv);
                }
            }
        }
    }
    std::cout << "CSV: " << path.string() << '\n';
    return failed ? 1 : 0;
}

} // namespace

int main(int argc, const char* argv[]) {
    try {
        Options options;
        Bench::CommandLine cli;
        ThreadManager::instance();
        cli.add(threading_command_line_contribution());
        cli.add(CommandLineContribution("evolution-recombination-benchmark", argument_parsers(options)));
        if (!cli.acquire(argc, argv))
            return 0;
        const Scenario scenarios[] = {
            {"urgent-crossing", EventKind::URGENT, 3.0, 4.0},
            {"permissive-jump", EventKind::PERMISSIVE, 1.0, 2.0}
        };
        if (options.sweep)
            return run_sweep(scenarios, options, cli);
        if (!options.window_widths.empty() || !options.distributed_counts.empty())
            throw std::invalid_argument("--window-widths / --distributed-counts require --sweep 1");
        if (options.distributed_count && !options.time_window)
            throw std::invalid_argument("--distributed-count requires --time-window");
        for (const auto& scenario: scenarios)
            run_scenario(scenario, options, cli, true);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "benchmark_recombination: " << error.what() << "\n";
        return 1;
    }
}
