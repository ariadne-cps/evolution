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
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
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
};

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

List<HybridEnclosure> generate(const Scenario& scenario, const Options& options,
                               DiscreteLocation source, DiscreteLocation target,
                               DiscreteEvent jump, RealVariable x, RealVariable y) {
    GeneralHybridEvolver evolver(make_system(scenario, source, target, jump, x, y));
    configure(evolver, options);
    const HybridBoxSet initial = scenario.kind == EventKind::URGENT
        ? HybridBoxSet(source, {Real(-2) <= x <= Real(-1), Real(-1 / 2_q) <= y <= Real(1 / 2_q)})
        : HybridBoxSet(source, {Real(1) <= x <= Real(2), Real(-1 / 2_q) <= y <= Real(1 / 2_q)});
    const HybridTerminationCriterion stop(Real(ExactDouble(scenario.generation_horizon)), Integer(8),
                                           Set<DiscreteEvent>{jump});
    const auto orbit = evolver.orbit(initial, stop, Semantics::UPPER);
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
    if (candidates.size() < 2u)
        throw std::runtime_error("not enough post-reset checkpoints; adjust integration step");
    // Select adjacent crossing/jump times rather than recombining a distant
    // collection whose geometric inflation would be predetermined.
    const SizeType take = std::min<SizeType>(candidates.size(), options.group_size);
    if (take < 2u)
        throw std::runtime_error("group-size must be at least two");
    const SizeType offset = (candidates.size() - take) / 2u;
    List<HybridEnclosure> selection;
    for (SizeType i = offset; i != offset + take; ++i)
        selection.append(candidates[i]);
    std::cout << scenario.name << ": generated=" << candidates.size()
              << " post-reset checkpoints, selected=" << selection.size()
              << " (middle adjacent time window)\n";
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
    check_reference_overlap(*m.analytical_final, *m.continued_independent, "independent");
    check_reference_overlap(*m.analytical_final, *m.continued_merged, "recombined");
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

void run_scenario(const Scenario& scenario, const Options& options, const Bench::CommandLine& cli) {
    const DiscreteLocation source(1), target(2);
    const DiscreteEvent jump("jump");
    const RealVariable x("x"), y("y");
    std::cout << "\nscenario=" << scenario.name << " step=" << options.step
              << " accuracy=" << options.accuracy << " group-size=" << options.group_size << "\n";
    const auto checkpoints = generate(scenario, options, source, target, jump, x, y);
    List<Enclosure> original;
    for (const auto& checkpoint: checkpoints)
        original.append(checkpoint.continuous_set());

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
            std::cout << "  sample=" << round - cli.warmups() + 1
                      << " independent-final=" << last->independent_final
                      << " recombined-final=" << last->merged_final
                      << " net-speedup=" << last->independent_total_ns() / last->merged_total_ns() << "\n";
        }
        ++round;
    });
    if (completed.size() != options.samples || totals.size() != options.samples)
        throw std::runtime_error("incorrect measured sample count");
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
        for (const auto& scenario: scenarios)
            run_scenario(scenario, options, cli);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "benchmark_recombination: " << error.what() << "\n";
        return 1;
    }
}
