/***************************************************************************
 *            benchmark_watertank.cpp
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

/** @file benchmarks/benchmark_watertank.cpp
 *  @brief Benchmark numerical reachability of the hybrid watertank model.
 */

#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

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
    double duration = 80.0, step = 1.0, radius = 0.5, accuracy = 1e-5;
    unsigned events = 8u, samples = 5u;
};
List<ArgumentParser> argument_parsers(Options& options) {
    List<ArgumentParser> parsers;
    auto positive = [&](const String& short_id, const String& name, const String& description, double& setting) {
        parsers.append(make_valued_argument_parser(short_id, name, description,
                                                   [&setting](const String& value) -> VoidFunction {
                                                       std::size_t used = 0u;
                                                       const std::string text(value);
                                                       const double number = std::stod(text, &used);
                                                       if (used != text.size() || !std::isfinite(number) || number <= 0.0)
                                                           throw std::invalid_argument("expected a finite positive number");
                                                       return [&setting, number] { setting = number; };
                                                   }));
    };
    positive("d", "duration", "Continuous time horizon (default: 80)", options.duration);
    positive("p", "step", "Maximum integration step (default: 1)", options.step);
    positive("r", "radius", "Maximum enclosure radius (default: 0.5)", options.radius);
    positive("a", "accuracy", "Integrator step maximum error (default: 1e-5)", options.accuracy);
    auto count = [&](const String& short_id, const String& name, const String& description, unsigned& setting, bool allow_zero) {
        parsers.append(make_valued_argument_parser(short_id, name, description,
                                                   [&setting, allow_zero](const String& value) -> VoidFunction {
                                                       const unsigned number = Bench::parse_unsigned<unsigned>(value);
                                                       if (!allow_zero && number == 0u)
                                                           throw std::invalid_argument("expected a positive count");
                                                       return [&setting, number] { setting = number; };
                                                   }));
    };
    count("e", "events", "Maximum discrete events (default: 8)", options.events, true);
    count("n", "samples", "Measured runs (default: 5)", options.samples, false);
    return parsers;
}
HybridAutomaton make_watertank() {
    RealVariable height("height"), aperture("aperture");
    StringVariable valve("valve");
    DiscreteLocation open(valve | "open"), opening(valve | "opening"), closed(valve | "closed"), closing(valve | "closing");
    HybridAutomaton tank("watertank");
    RealExpression flow = -Real(1 / 50_q) * height + Real(3 / 10_q) * aperture;
    tank.new_mode(open, {dot(height) = flow, dot(aperture) = 0});
    tank.new_mode(closed, {dot(height) = flow, dot(aperture) = 0});
    tank.new_mode(opening, {dot(height) = flow, dot(aperture) = Real(1 / 4_q)});
    tank.new_mode(closing, {dot(height) = flow, dot(aperture) = -Real(1 / 4_q)});
    tank.new_transition(closed, DiscreteEvent("start_opening"), opening, {next(height) = height, next(aperture) = aperture}, height <= Real(11 / 2_q), EventKind::URGENT);
    tank.new_transition(open, DiscreteEvent("start_closing"), closing, {next(height) = height, next(aperture) = aperture}, height >= 8, EventKind::URGENT);
    tank.new_transition(opening, DiscreteEvent("finished_opening"), open, {next(height) = height, next(aperture) = aperture}, aperture >= 1, EventKind::URGENT);
    tank.new_transition(closing, DiscreteEvent("finished_closing"), closed, {next(height) = height, next(aperture) = aperture}, aperture <= 0, EventKind::URGENT);
    return tank;
}
// Keep timing boundaries identical: only evolver.orbit is measured.
// The Runner executes warm-ups, but the separate phase clock excludes
// verification and Logging just as the original experiment did.
struct Measurement {
    double nanoseconds;
    SizeType reach, intermediate, final;
};

Measurement measure(const GeneralHybridEvolver& evolver,
                    const HybridSet& initial, const HybridTime& horizon,
                    bool describe) {
    const auto start = Bench::Clock::now();
    auto orbit = evolver.orbit(initial, horizon, Semantics::UPPER);
    const auto stop = Bench::Clock::now();
    if (orbit.final().size() == 0u)
        throw std::runtime_error("evolution produced no final enclosures");
    if (describe) {
        for (const auto& enclosure: orbit.final()) {
            LOGGING_PRINTLN_AT(1, "final location=" << enclosure.location()
                                                    << " time=" << enclosure.time_range()
                                                    << " events=" << enclosure.previous_events()
                                                    << " bounds=" << enclosure.state_bounding_box())
        }
    }
    return {Bench::elapsed_ns(start, stop), orbit.reach().size(),
            orbit.intermediate().size(), orbit.final().size()};
}

} // namespace
int main(int argc, const char* argv[]) {
    try {
        Options options;
        Bench::CommandLine cli;
        ThreadManager::instance();
        cli.add(threading_command_line_contribution());
        cli.add(CommandLineContribution("evolution-watertank-benchmark",
                                        argument_parsers(options)));
        if (!cli.acquire(argc, argv))
            return 0;
        const HybridAutomaton tank = make_watertank();
        RealVariable height("height"), aperture("aperture");
        StringVariable valve("valve");
        HybridSet initial(DiscreteLocation(valve | "opening"),
                          {0 <= height <= Real(1 / 16_q), 0 <= aperture <= Real(1 / 64_q)});
        HybridTime horizon(ExactDouble(options.duration), options.events);
        GeneralHybridEvolver evolver(tank);
        evolver.set_integrator(
            TaylorPicardIntegrator(StepMaximumError(options.accuracy)));
        evolver.configuration().set_maximum_step_size(
            ApproximateDouble(options.step));
        evolver.configuration().set_maximum_enclosure_radius(
            ApproximateDouble(options.radius));

        LOGGING_PRINTLN("Watertank benchmark | duration=" << options.duration
                                                          << " | events=" << options.events
                                                          << " | step=" << options.step << " | radius=" << options.radius
                                                          << " | accuracy=" << options.accuracy
                                                          << " | concurrency=" << ThreadManager::instance().concurrency()
                                                          << " | samples=" << options.samples << " | warmups=" << cli.warmups())
        LOGGING_PRINTLN("Initial: " << initial)

        Bench::Runner runner(cli);
        Bench::Samples durations;
        std::optional<Measurement> last;
        unsigned round = 0u;
        // A Runner round is one complete orbit. The sample's elapsed time is
        // recorded by the inner phase clock, not the Runner's encompassing
        // timer (which includes validation and diagnostic logging).
        const auto completed = runner.run("watertank", options.samples, [] {}, [&] {
                const bool measured=round>=cli.warmups();
                last=measure(evolver,initial,horizon,measured); }, [&] {
                if(not last)throw std::runtime_error("missing watertank result");
                if(round>=cli.warmups()) {
                    durations.add(last->nanoseconds);
                    LOGGING_PRINTLN_AT(1,"sample="<<round-cli.warmups()+1u
                        <<" | ms="<<Bench::Table::number(last->nanoseconds/1e6)
                        <<" | reach="<<last->reach
                        <<" | intermediate="<<last->intermediate
                        <<" | final="<<last->final)
                }
                ++round; });
        if (completed.size() != options.samples ||
            durations.size() != options.samples)
            throw std::runtime_error("incorrect watertank sample count");

        Bench::Report report("Watertank orbit duration (nanoseconds)",
                             cli.summary_fields());
        report.set_parameter_header(
            "One sample = one complete orbit; timing excludes correctness "
            "checks and enclosure diagnostic output.");
        report.set_dimension_columns({
            {"Scenario"},
            {"Duration", 0u, Bench::TableAlignment::Right},
            {"Events", 0u, Bench::TableAlignment::Right},
            {"Samples", 0u, Bench::TableAlignment::Right}
        });
        report.set_table_options(Bench::TableOptions::ruled());
        report.add_row({"watertank", Bench::Table::number(options.duration),
                        std::to_string(options.events),
                        std::to_string(options.samples)},
                       Bench::summarize(durations));
        report.print();
        return 0;
    } catch (const std::exception& error) {
        LOGGING_PRINTLN("benchmark_watertank: " << error.what())
        return 1;
    }
}
