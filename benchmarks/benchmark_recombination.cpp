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

struct Measurement {
    double contraction_ns = 0, recombination_ns = 0;
    double independent_ns = 0, merged_ns = 0;
    SizeType independent_final = 0, merged_final = 0;
    EnclosureRecombinationResult result;
};

Measurement measure(const List<Enclosure>& original,
                    const List<HybridEnclosure>& checkpoints,
                    const GeneralHybridEvolver& continuation,
                    const HybridTerminationCriterion& stop,
                    DiscreteLocation target, RealSpace space) {
    Measurement m;
    // Copying the immutable baseline is preparation, not measured contraction.
    List<Enclosure> contracted = original;
    const auto contraction_start = Bench::Clock::now();
    for (const auto& enclosure: contracted)
        enclosure.reduce();
    const auto contraction_end = Bench::Clock::now();
    const auto recombination_start = Bench::Clock::now();
    m.result = ReferenceEnclosureRecombiner().recombine(contracted);
    const auto recombination_end = Bench::Clock::now();
    m.contraction_ns = Bench::elapsed_ns(contraction_start, contraction_end);
    m.recombination_ns = Bench::elapsed_ns(recombination_start, recombination_end);
    if (m.result.outputs.size() != 1u)
        throw std::runtime_error("expected one nonempty reference recombination output");

    // The continuation system contains only the target dynamics: neither
    // collection can jump again, and no termination event history is erased.
    // Recombination creates a new trajectory; its provenance remains in result.
    const HybridEnclosure recombined(target, space, m.result.outputs.front().enclosure);
    const auto independent_start = Bench::Clock::now();
    for (const auto& checkpoint: checkpoints) {
        const auto orbit = continuation.orbit(checkpoint, stop, Semantics::UPPER);
        m.independent_final += orbit.final().size();
    }
    const auto independent_end = Bench::Clock::now();
    const auto merged_start = Bench::Clock::now();
    const auto merged_orbit = continuation.orbit(recombined, stop, Semantics::UPPER);
    const auto merged_end = Bench::Clock::now();
    m.independent_ns = Bench::elapsed_ns(independent_start, independent_end);
    m.merged_ns = Bench::elapsed_ns(merged_start, merged_end);
    m.merged_final = merged_orbit.final().size();
    return m;
}

double width(const UpperIntervalType& value) {
    return value.upper_bound().raw().get_d() - value.lower_bound().raw().get_d();
}

void describe_quality(const List<Enclosure>& original, const EnclosureRecombinationResult& result) {
    const auto& output = result.outputs.front();
    UpperBoxType input_state = apply(original.front().state_function(), original.front().reduced_domain());
    UpperIntervalType input_time = apply(original.front().time_function(), original.front().reduced_domain());
    UpperIntervalType input_dwell = apply(original.front().dwell_time_function(), original.front().reduced_domain());
    for (SizeType i = 1; i != original.size(); ++i) {
        const auto domain = original[i].reduced_domain();
        input_state = hull(input_state, apply(original[i].state_function(), domain));
        input_time = hull(input_time, apply(original[i].time_function(), domain));
        input_dwell = hull(input_dwell, apply(original[i].dwell_time_function(), domain));
    }
    const auto& merged = output.enclosure;
    const auto domain = merged.reduced_domain();
    const auto output_state = apply(merged.state_function(), domain);
    const auto output_time = apply(merged.time_function(), domain);
    const auto output_dwell = apply(merged.dwell_time_function(), domain);
    std::cout << std::setprecision(6)
              << "  provenance=" << output.source_indices.size()
              << " discarded=" << result.discarded_sources.size()
              << " parameters=" << merged.number_of_parameters()
              << " constraints=" << merged.number_of_constraints()
              << " radius-limit-satisfied=" << result.radius_limit_satisfied << "\n";
    for (SizeType k = 0; k != input_state.size(); ++k) {
        std::cout << "  state[" << k << "] hull-width=" << width(input_state[k])
                  << " result-width=" << width(output_state[k]) << "\n";
    }
    std::cout << "  time hull-width=" << width(input_time)
              << " result-width=" << width(output_time)
              << " dwell hull-width=" << width(input_dwell)
              << " result-width=" << width(output_dwell) << "\n";
    // These interval widths measure wrapping/inflation, not certified set
    // inclusion. Soundness is the recombiner's contract, tested separately.
}

void run_scenario(const Scenario& scenario, const Options& options, const Bench::CommandLine& cli) {
    const DiscreteLocation source(1), target(2);
    const DiscreteEvent jump("jump");
    const RealVariable x("x"), y("y");
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
    std::optional<Measurement> last;
    unsigned round = 0u;
    const auto completed = runner.run(scenario.name, options.samples, [] {}, [&] {
        last = measure(original, checkpoints, continuation, stop, target, RealSpace({x, y}));
    }, [&] {
        if (!last)
            throw std::runtime_error("missing measurement");
        if (round >= cli.warmups()) {
            contractions.add(last->contraction_ns);
            recombinations.add(last->recombination_ns);
            totals.add(last->contraction_ns + last->recombination_ns);
            independents.add(last->independent_ns);
            merged.add(last->merged_ns);
            std::cout << "  sample=" << round - cli.warmups() + 1
                      << " independent-final=" << last->independent_final
                      << " recombined-final=" << last->merged_final << "\n";
        }
        ++round;
    });
    if (completed.size() != options.samples || totals.size() != options.samples)
        throw std::runtime_error("incorrect measured sample count");
    describe_quality(original, last->result);

    Bench::Report report(std::string("Recombination timings: ") + scenario.name + " (nanoseconds)",
                         cli.summary_fields());
    report.set_parameter_header("One sample = fresh constraint reduction, reference recombination, "
                                "independent and merged target-only continuation. "
                                "Dataset generation, copying, verification and reporting are untimed.");
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
