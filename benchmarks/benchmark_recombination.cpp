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
    UpperIntervalType time, dwell;
};

void accumulate(std::optional<Projection>& total, const Enclosure& enclosure) {
    const auto domain = enclosure.reduced_domain();
    Projection current{apply(enclosure.state_function(), domain),
                       apply(enclosure.time_function(), domain),
                       apply(enclosure.dwell_time_function(), domain)};
    if (!total) {
        total = std::move(current);
    } else {
        total->state = hull(total->state, current.state);
        total->time = hull(total->time, current.time);
        total->dwell = hull(total->dwell, current.dwell);
    }
}

struct Measurement {
    double contraction_ns = 0, recombination_ns = 0;
    double independent_ns = 0, merged_ns = 0;
    SizeType independent_final = 0, merged_final = 0;
    EnclosureRecombinationResult result;
    std::optional<Projection> source, merged, continued_independent, continued_merged;
    double independent_total_ns() const { return contraction_ns + independent_ns; }
    double merged_total_ns() const { return contraction_ns + recombination_ns + merged_ns; }
};

Measurement measure(const List<Enclosure>& original,
                    const GeneralHybridEvolver& continuation,
                    const HybridTerminationCriterion& stop,
                    DiscreteLocation target, RealSpace space) {
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

    for (const auto& source: contracted)
        accumulate(m.source, source);
    accumulate(m.merged, m.result.outputs.front().enclosure);
    for (const auto& orbit: independent_orbits) {
        m.independent_final += orbit.final().size();
        for (const auto& enclosure: orbit.final())
            accumulate(m.continued_independent, enclosure.continuous_set());
    }
    m.merged_final = merged_orbit.final().size();
    for (const auto& enclosure: merged_orbit.final())
        accumulate(m.continued_merged, enclosure.continuous_set());
    if (!m.source || !m.merged || !m.continued_independent || !m.continued_merged)
        throw std::runtime_error("missing source or final projection");
    return m;
}

double width(const UpperIntervalType& value) {
    return value.upper_bound().raw().get_d() - value.lower_bound().raw().get_d();
}

void compare_quality(const char* phase, const Projection& independent,
                     const Projection& merged) {
    std::cout << "  " << phase << " (validated projection hulls):\n";
    for (SizeType k = 0; k != independent.state.size(); ++k) {
        const double before = width(independent.state[k]), after = width(merged.state[k]);
        std::cout << "    state[" << k << "] independent-width=" << before
                  << " merged-width=" << after << " ratio="
                  << (before > 0 ? after / before : 0) << "\n";
    }
    const double bt = width(independent.time), at = width(merged.time);
    const double bd = width(independent.dwell), ad = width(merged.dwell);
    std::cout << "    time independent-width=" << bt << " merged-width=" << at
              << " ratio=" << (bt > 0 ? at / bt : 0) << "\n";
    std::cout << "    dwell independent-width=" << bd << " merged-width=" << ad
              << " ratio=" << (bd > 0 ? ad / bd : 0) << "\n";
    // Hull containment is necessary for exact hulls, but neither necessary
    // nor sufficient for constrained-source inclusion with outer ranges.
    std::cout << "    hull-containment (diagnostic only): state="
              << definitely(subset(independent.state, merged.state))
              << " time=" << definitely(subset(independent.time, merged.time))
              << " dwell=" << definitely(subset(independent.dwell, merged.dwell))
              << "\n";
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
        last = measure(original, continuation, stop, target, RealSpace({x, y}));
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
