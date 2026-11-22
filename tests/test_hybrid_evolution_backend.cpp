/***************************************************************************
 *            test_hybrid_evolution_backend.cpp
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

/** @file tests/test_hybrid_evolution_backend.cpp
 *  @brief Test backend independence, publication boundaries and Orbit roles.
 */

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "evolution/hybrid_continuous_time_evolver.hpp"
#include "evolution/hybrid_evolution_backend.hpp"
#include "evolution/hybrid_evolution_tree.hpp"
#include "logging/test_runner.hpp"
#include "system/hybrid_automata.hpp"
#include "threading/thread_manager.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

class ConcurrencyScope {
    std::size_t _previous = ThreadManager::instance().concurrency();

  public:
    explicit ConcurrencyScope(std::size_t count) { ThreadManager::instance().set_concurrency(count); }
    ~ConcurrencyScope() { ThreadManager::instance().set_concurrency(_previous); }
};

class StepProbe : public HybridContinuousTimeEvolver<ValidatedTag> {
  public:
    using HybridContinuousTimeEvolver<ValidatedTag>::HybridContinuousTimeEvolver;
    using HybridContinuousTimeEvolverBase<ValidatedTag>::_process_working_set;
};

// Independent LIFO executor: a test consumer of the same numerical callback.
class StackBackend : public HybridEvolutionBackendInterface {
    Orbit<HybridEnclosure> _orbit;
    std::vector<HybridEvolutionWorkItem> _pending;

  public:
    explicit StackBackend(const HybridEnclosure& initial)
        : _orbit(initial), _pending{
                               {initial, HybridEvolutionPhase::PROCESS_STARTING_EVENTS}
    } {}
    Void run(const StepFunction& advance) override {
        while (!_pending.empty()) {
            auto item = std::move(_pending.back());
            _pending.pop_back();
            auto result = advance(item, {_orbit.reach().size(), _orbit.intermediate().size(), _orbit.final().size()});
            for (const auto& contribution: result.contributions) {
                switch (contribution.kind) {
                    case HybridEvolutionSetKind::REACH:
                        _orbit.adjoin_reach(contribution.enclosure);
                        break;
                    case HybridEvolutionSetKind::INTERMEDIATE:
                        _orbit.adjoin_intermediate(contribution.enclosure);
                        break;
                    case HybridEvolutionSetKind::FINAL:
                        _orbit.adjoin_final(contribution.enclosure);
                        break;
                    default:
                        throw std::logic_error("Unexpected contribution");
                }
            }
            for (auto& successor: result.successors)
                _pending.push_back(std::move(successor.work));
        }
    }
    Orbit<HybridEnclosure> release_orbit() override { return std::move(_orbit); }
};

class StackFactory : public HybridEvolutionBackendFactoryInterface {
    std::shared_ptr<std::atomic<std::size_t>> _created;

  public:
    explicit StackFactory(std::shared_ptr<std::atomic<std::size_t>> created) : _created(std::move(created)) {}
    StackFactory* clone() const override { return new StackFactory(*this); }
    std::unique_ptr<HybridEvolutionBackendInterface> create(const HybridEnclosure& initial, const Real&) const override {
        ++*_created;
        return std::make_unique<StackBackend>(initial);
    }
};

std::vector<std::string> descriptions(const ListSet<HybridEnclosure>& sets) {
    std::vector<std::string> result;
    for (const auto& set: sets) {
        std::ostringstream stream;
        stream << set;
        result.push_back(stream.str());
    }
    std::sort(result.begin(), result.end());
    return result;
}

} // namespace

class TestHybridEvolutionBackend : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    void test_progress_revises_remaining_work_without_recounting_history() {
        DiscreteLocation q(1);
        RealVariable x("x");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1});
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        const auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                                     x
        }),
                                                              ExactBoxType({{0, 0}})));
        auto at_four = initial;
        at_four.set_time_function(evolver.function_factory().create_constant(initial.parameter_domain(), Real(4)));
        HybridEvolutionProgressTracker tracker(10.0);
        tracker.add_pending(initial);
        ARIADNE_TESTING_EXPECT_EXCEPTION(tracker.finish());
        HybridEvolutionStepResult flow;
        flow.continuous_evolution_time = 4.0;
        flow.continue_with(at_four);
        flow.record_reach(at_four);
        flow.record_intermediate(at_four);
        tracker.publish({initial, HybridEvolutionPhase::EVOLVE_CONTINUOUSLY}, flow);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().processed_time, 4.0);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().remaining_time.value(), 6.0);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().completion_rate(), 0.4);

        HybridEvolutionStepResult split;
        split.continue_with(at_four, HybridEvolutionSuccessorKind::SUBDIVISION);
        split.continue_with(at_four, HybridEvolutionSuccessorKind::SUBDIVISION);
        tracker.publish({at_four, HybridEvolutionPhase::EVOLVE_CONTINUOUSLY}, split);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().processed_time, 4.0);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().pending, 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().remaining_time.value(), 12.0);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().completion_rate(), 0.4);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().credited_time, 8.0);

        HybridEvolutionStepResult recondition;
        recondition.continue_with(at_four, HybridEvolutionSuccessorKind::RECONDITIONING);
        tracker.publish({at_four, HybridEvolutionPhase::EVOLVE_CONTINUOUSLY}, recondition);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().processed_time, 4.0);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().remaining_time.value(), 12.0);
        // A committed 2-to-1 replacement changes only outstanding work.
        tracker.add_pending(at_four);
        tracker.remove_pending(at_four);
        tracker.remove_pending(at_four);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().pending, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().processed_time, 4.0);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().remaining_time.value(), 6.0);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().completion_rate(), 0.4);

        HybridEvolutionStepResult empty;
        empty.stop_reason = HybridEvolutionStopReason::EMPTY;
        tracker.publish({at_four, HybridEvolutionPhase::EVOLVE_CONTINUOUSLY}, empty);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().processed_time, 4.0);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().remaining_time.value(), 0.0);
        ARIADNE_TESTING_PROVE(!tracker.snapshot().complete);
        tracker.finish();
        ARIADNE_TESTING_PROVE(tracker.snapshot().complete == true);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().completion_rate(), 1.0);
    }

    void test_progress_retains_the_120_second_prefix_when_permissive_branches_accumulate() {
        DiscreteLocation q(1);
        RealVariable x("x");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1});
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        const auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                                     x
        }),
                                                              ExactBoxType({{0, 0}})));
        const auto at_time = [&](Real time) {
            auto set = initial;
            set.set_time_function(evolver.function_factory().create_constant(initial.parameter_domain(), time));
            return set;
        };
        HybridEvolutionProgressTracker tracker(200.0);
        tracker.add_pending(initial);
        auto current = at_time(Real(120));
        HybridEvolutionStepResult first;
        first.continuous_evolution_time = 120.0;
        first.continue_with(current);
        tracker.publish({initial, HybridEvolutionPhase::EVOLVE_CONTINUOUSLY}, first);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().completion_rate(), 0.6);
        for (unsigned n = 1u; n <= 500u; ++n) {
            auto next = at_time(Real(120) + Real(n) / 1024);
            HybridEvolutionStepResult branching;
            branching.continuous_evolution_time = 1.0 / 1024.0;
            branching.transition_to(next, DiscreteEvent("timeout"), EventKind::PERMISSIVE);
            branching.continue_with(next);
            tracker.publish({current, HybridEvolutionPhase::EVOLVE_CONTINUOUSLY}, branching);
            current = std::move(next);
        }
        const auto progress = tracker.snapshot();
        ARIADNE_TESTING_EXPECT_EQUAL(progress.pending, 501u);
        ARIADNE_TESTING_EXPECT_EQUAL(static_cast<unsigned>(std::round(progress.completion_rate() * 100)), 60u);
        ARIADNE_TESTING_PROVE(progress.completion_rate() >= 0.6);
        ARIADNE_TESTING_PROVE(progress.completion_rate() < 0.603);
    }

    void test_temporal_progress_uses_the_initial_time_and_terminal_branch_credit() {
        DiscreteLocation q(1);
        RealVariable x("x");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1});
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                               x
        }),
                                                        ExactBoxType({{0, 0}})));
        initial.set_time_function(evolver.function_factory().create_constant(initial.parameter_domain(), Real(40)));
        auto later = initial;
        later.set_time_function(evolver.function_factory().create_constant(initial.parameter_domain(), Real(120)));
        HybridEvolutionProgressTracker tracker(200.0);
        tracker.add_pending(initial);
        HybridEvolutionStepResult split;
        split.continuous_evolution_time = 80.0;
        split.continue_with(later, HybridEvolutionSuccessorKind::SUBDIVISION);
        split.continue_with(later, HybridEvolutionSuccessorKind::SUBDIVISION);
        tracker.publish({initial, HybridEvolutionPhase::EVOLVE_CONTINUOUSLY}, split);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().completion_rate(), 0.5);
        tracker.publish({later, HybridEvolutionPhase::EVOLVE_CONTINUOUSLY}, {});
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().processed_time, 80.0);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().credited_time, 240.0);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().remaining_time.value(), 80.0);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().completion_rate(), 0.75);
    }

    void test_invalid_cost_metadata_does_not_reset_valid_temporal_progress() {
        DiscreteLocation q(1);
        RealVariable x("x");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1});
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        const auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                                     x
        }),
                                                              ExactBoxType({{0, 0}})));
        auto next = initial;
        next.set_time_function(evolver.function_factory().create_constant(initial.parameter_domain(), Real(120)));
        HybridEvolutionProgressTracker tracker(200.0);
        tracker.add_pending(initial);
        HybridEvolutionStepResult result;
        result.continuous_evolution_time = std::numeric_limits<double>::quiet_NaN();
        result.continue_with(next);
        tracker.publish({initial, HybridEvolutionPhase::EVOLVE_CONTINUOUSLY}, result);
        ARIADNE_TESTING_PROVE(std::isnan(tracker.snapshot().processed_time));
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().completion_rate(), 0.6);
    }

    void test_progress_handles_zero_and_unknown_horizons() {
        DiscreteLocation q(1);
        RealVariable x("x");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1});
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        const auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                                     x
        }),
                                                              ExactBoxType({{0, 0}})));
        for (auto horizon: {std::optional<double>(0.0), std::optional<double>{},
                            std::optional<double>(std::numeric_limits<double>::infinity())}) {
            HybridEvolutionProgressTracker tracker(horizon);
            tracker.add_pending(initial);
            ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().completion_rate(), 0.0);
            ARIADNE_TESTING_PROVE(tracker.snapshot().remaining_time.has_value() == (horizon && *horizon == 0.0));
            tracker.publish({initial, HybridEvolutionPhase::PROCESS_STARTING_EVENTS}, {});
            tracker.finish();
            ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().processed_time, 0.0);
            ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().completion_rate(), 1.0);
            ARIADNE_TESTING_EXPECT_EXCEPTION(tracker.remove_pending(initial));
        }
    }

    void test_progress_uses_lower_clock_bounds_for_source_and_replacement_estimates() {
        DiscreteLocation q(1);
        RealVariable x("x");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1});
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        const auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                                     x
        }),
                                                              ExactBoxType({{0, 1}})));
        auto first = initial, second = initial, merged = initial;
        first.set_time_function(initial.state_function().get(0) + 4);
        second.set_time_function(initial.state_function().get(0) + 7);
        merged.set_time_function(4 * initial.state_function().get(0) + 4);
        HybridEvolutionProgressTracker tracker(10.0);
        tracker.add_pending(first);  // Time [4,5], remaining 6.
        tracker.add_pending(second); // Time [7,8], remaining 3.
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().remaining_time.value(), 9.0);
        tracker.add_pending(merged); // Time [4,8], remaining 6, including uncertainty introduced by merging.
        tracker.remove_pending(first);
        tracker.remove_pending(second);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().remaining_time.value(), 6.0);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().processed_time, 0.0);
        ARIADNE_TESTING_EXPECT_EQUAL(tracker.snapshot().pending, 1u);
    }

    void test_flat_and_tree_report_the_same_time_work_after_splitting() {
        ConcurrencyScope concurrency(0u);
        DiscreteLocation q(1);
        RealVariable x("x");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1});
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        const auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                                     x
        }),
                                                              ExactBoxType({{0, 0}})));
        auto at_four = initial, at_ten = initial;
        at_four.set_time_function(evolver.function_factory().create_constant(initial.parameter_domain(), Real(4)));
        at_ten.set_time_function(evolver.function_factory().create_constant(initial.parameter_domain(), Real(10)));
        std::vector<std::shared_ptr<const HybridEvolutionBackendFactoryInterface>> factories{
            std::make_shared<FlatHybridEvolutionBackendFactory>(),
            std::make_shared<TreeHybridEvolutionBackendFactory>(),
            std::make_shared<TreeHybridEvolutionBackendFactory>(ImmediateCheckpointAdmissionPolicy{}),
            std::make_shared<TreeHybridEvolutionBackendFactory>(LocationCheckpointAdmissionPolicy{})};
        for (const auto& factory: factories) {
            auto backend = factory->create(initial, Real(10));
            bool split = false;
            std::size_t pieces = 0u;
            backend->run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
                const auto progress = backend->progress();
                HybridEvolutionStepResult result;
                if (item.phase == HybridEvolutionPhase::PROCESS_STARTING_EVENTS) {
                    result.continue_with(initial);
                } else if (item.enclosure.time_range().lower_bound().get_d() == 0.0) {
                    result.continuous_evolution_time = 4.0;
                    result.record_intermediate(at_ten); // Output roles must not inflate processed time.
                    result.continue_with(at_four);
                } else if (!split) {
                    ARIADNE_TESTING_EXPECT_EQUAL(progress.completion_rate(), 0.4);
                    split = true;
                    result.continue_with(at_four, HybridEvolutionSuccessorKind::SUBDIVISION);
                    result.continue_with(at_four, HybridEvolutionSuccessorKind::SUBDIVISION);
                } else {
                    ARIADNE_TESTING_EXPECT_EQUAL(progress.completion_rate(), pieces == 0u ? 0.4 : 0.7);
                    ARIADNE_TESTING_EXPECT_EQUAL(progress.pending, 2u - pieces);
                    ARIADNE_TESTING_EXPECT_EQUAL(progress.remaining_time.value(), pieces == 0u ? 12.0 : 6.0);
                    ARIADNE_TESTING_EXPECT_EQUAL(progress.processed_time, pieces == 0u ? 4.0 : 10.0);
                    ++pieces;
                    result.continuous_evolution_time = 6.0;
                    result.record_final(at_ten);
                }
                return result;
            });
            const auto progress = backend->progress();
            ARIADNE_TESTING_EXPECT_EQUAL(pieces, 2u);
            ARIADNE_TESTING_EXPECT_EQUAL(progress.processed_time, 16.0);
            ARIADNE_TESTING_EXPECT_EQUAL(progress.remaining_time.value(), 0.0);
            ARIADNE_TESTING_EXPECT_EQUAL(progress.pending, 0u);
            ARIADNE_TESTING_PROVE(progress.complete);
            ARIADNE_TESTING_EXPECT_EQUAL(progress.completion_rate(), 1.0);
        }
    }

    void test_numerical_progress_uses_the_selected_continuous_duration() {
        DiscreteLocation q(1);
        RealVariable x("x");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1});
        StepProbe probe(system);
        probe.configuration().set_maximum_step_size(0.25_x);
        const auto initial = probe.enclosure(HybridExactBox(q, RealSpace({
                                                                   x
        }),
                                                            ExactBoxType({{0, 0}})));
        const auto start = probe._process_working_set({initial, HybridEvolutionPhase::PROCESS_STARTING_EVENTS}, HybridTime(2, 1), Semantics::UPPER, {});
        ARIADNE_TESTING_EXPECT_EQUAL(start.continuous_evolution_time, 0.0);
        const auto step = probe._process_working_set({initial, HybridEvolutionPhase::EVOLVE_CONTINUOUSLY}, HybridTime(2, 1), Semantics::UPPER, {});
        ARIADNE_TESTING_EXPECT_EQUAL(step.continuous_evolution_time, 0.25);
        const auto last = probe._process_working_set({initial, HybridEvolutionPhase::EVOLVE_CONTINUOUSLY}, HybridTime(0.125_x, 1), Semantics::UPPER, {});
        ARIADNE_TESTING_EXPECT_EQUAL(last.continuous_evolution_time, 0.125);
    }

    void test_output_and_work_are_distinct() {
        DiscreteLocation q(1);
        RealVariable x("x"), y("y");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0});
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                               x, y
        }),
                                                        ExactBoxType({{0, 0}, {0, 0}})));
        auto next = initial, output = initial;
        next.set_time_function(evolver.function_factory().create_constant(initial.parameter_domain(), Real(1)));
        output.set_time_function(evolver.function_factory().create_constant(initial.parameter_domain(), Real(9)));
        for (std::size_t concurrency: {0u, 2u}) {
            ConcurrencyScope scope(concurrency);
            FlatHybridEvolutionBackendFactory factory;
            auto backend = factory.create(initial, Real(10));
            ARIADNE_TESTING_EXPECT_EXCEPTION(backend->release_orbit());
            std::atomic<std::size_t> calls{0};
            HybridEvolutionStatistics observed;
            double observed_time = -1;
            backend->run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics& stats) {
                ++calls;
                HybridEvolutionStepResult result;
                if (item.phase == HybridEvolutionPhase::PROCESS_STARTING_EVENTS) {
                    result.record_reach(initial);
                    result.record_intermediate(output);
                    result.record_final(output);
                    result.continue_with(next);
                } else {
                    observed = stats;
                    observed_time = item.enclosure.time_range().lower_bound().get_d();
                    result.stop_reason = HybridEvolutionStopReason::EMPTY;
                }
                return result;
            });
            auto orbit = backend->release_orbit();
            ARIADNE_TESTING_EXPECT_EQUAL(calls.load(), 2u);
            ARIADNE_TESTING_EXPECT_EQUAL(observed.reach, 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(observed.intermediate, 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(observed.final, 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(observed_time, 1.0);
            ARIADNE_TESTING_EXPECT_EQUAL(orbit.reach().size(), 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(orbit.intermediate().size(), 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(orbit.initial().time_range().lower_bound().get_d(), 0.0);
            ARIADNE_TESTING_EXPECT_EQUAL(orbit.intermediate().begin()->time_range().lower_bound().get_d(), 9.0);
            ARIADNE_TESTING_EXPECT_EXCEPTION(backend->release_orbit());
            ARIADNE_TESTING_EXPECT_EXCEPTION(backend->run({}));
        }
    }

    void test_failed_step_cannot_release_an_orbit() {
        DiscreteLocation q(1);
        RealVariable x("x"), y("y");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0});
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                               x, y
        }),
                                                        ExactBoxType({{0, 0}, {0, 0}})));
        FlatHybridEvolutionBackendFactory factory;
        for (std::size_t concurrency: {0u, 2u}) {
            ConcurrencyScope scope(concurrency);
            auto backend = factory.create(initial, Real(1));
            ARIADNE_TESTING_EXPECT_EXCEPTION(backend->run([](const HybridEvolutionWorkItem&, const HybridEvolutionStatistics&) -> HybridEvolutionStepResult {
                throw std::runtime_error("Numerical step failed");
            }));
            ARIADNE_TESTING_EXPECT_EXCEPTION(backend->release_orbit());
            ARIADNE_TESTING_PROVE(!backend->progress().complete);
            ARIADNE_TESTING_EXPECT_EQUAL(backend->progress().pending, 1u);
        }
    }

    void test_evolver_accepts_an_independent_scheduler() {
        ConcurrencyScope scope(0);
        for (auto kind: {EventKind::URGENT, EventKind::PERMISSIVE}) {
            DiscreteLocation source(1), target(2);
            DiscreteEvent event("jump");
            RealVariable x("x"), y("y");
            HybridAutomaton system;
            system.new_mode(source, {dot(x) = 1, dot(y) = 0});
            system.new_mode(target, {dot(x) = 0, dot(y) = 0});
            system.new_transition(source, event, target, {next(x) = x, next(y) = y}, x >= 0, kind);
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            evolver.configuration().set_maximum_step_size(0.25_x);
            auto initial = evolver.enclosure(HybridExactBox(source, RealSpace({
                                                                        x, y
            }),
                                                            ExactBoxType({{-1, -0.75_x}, {0, 0}})));
            HybridTime horizon(2, 1);
            auto expected = evolver.orbit(initial, horizon, Semantics::UPPER);
            auto created = std::make_shared<std::atomic<std::size_t>>(0);
            {
                StackFactory factory(created);
                evolver.set_evolution_backend_factory(factory);
            }
            for (unsigned int run = 0; run != 2; ++run) {
                auto actual = evolver.orbit(initial, horizon, Semantics::UPPER);
                ARIADNE_TESTING_PROVE(descriptions(actual.reach()) == descriptions(expected.reach()));
                ARIADNE_TESTING_PROVE(descriptions(actual.intermediate()) == descriptions(expected.intermediate()));
                ARIADNE_TESTING_PROVE(descriptions(actual.final()) == descriptions(expected.final()));
            }
            ARIADNE_TESTING_EXPECT_EQUAL(created->load(), 2u);
        }
    }

    void test_lower_initial_activation_requires_the_whole_source() {
        const DiscreteLocation source(1), target(2);
        const DiscreteEvent event("activate");
        const RealVariable x("x");
        const RealSpace space({x});
        for (auto kind: {EventKind::URGENT, EventKind::PERMISSIVE}) {
            HybridAutomaton system;
            system.new_mode(source, {dot(x) = 0});
            system.new_mode(target, {dot(x) = 0});
            system.new_transition(source, event, target, {next(x) = x}, x >= 0, kind);
            StepProbe probe(system);
            HybridTime horizon(1, 1);
            const auto fully_active = probe.enclosure(
                HybridExactBox(source, space, ExactBoxType({
                                                  {0.25_x, 0.5_x}
            })));
            const auto partially_active = probe.enclosure(
                HybridExactBox(source, space, ExactBoxType({
                                                  {-0.25_x, 0.25_x}
            })));

            const auto full = probe._process_working_set(
                {fully_active, HybridEvolutionPhase::PROCESS_STARTING_EVENTS},
                horizon, Semantics::LOWER, {});
            bool full_transition = false;
            for (const auto& successor: full.successors) {
                if (successor.kind == HybridEvolutionSuccessorKind::TRANSITION) {
                    full_transition = true;
                    ARIADNE_TESTING_EXPECT_EQUAL(successor.event, event);
                    ARIADNE_TESTING_EXPECT_EQUAL(successor.event_kind, kind);
                    if (kind == EventKind::URGENT)
                        ARIADNE_TESTING_EXPECT_EQUAL(successor.resolution,
                                                     HybridEvolutionTransitionResolution::FULL_SOURCE);
                }
            }
            ARIADNE_TESTING_PROVE(full_transition);

            const auto lower = probe._process_working_set(
                {partially_active, HybridEvolutionPhase::PROCESS_STARTING_EVENTS},
                horizon, Semantics::LOWER, {});
            for (const auto& successor: lower.successors)
                ARIADNE_TESTING_PROVE(successor.kind != HybridEvolutionSuccessorKind::TRANSITION);
            ARIADNE_TESTING_PROVE(!lower.successors.empty());

            const auto upper = probe._process_working_set(
                {partially_active, HybridEvolutionPhase::PROCESS_STARTING_EVENTS},
                horizon, Semantics::UPPER, {});
            bool possible_transition = false;
            for (const auto& successor: upper.successors) {
                if (successor.kind == HybridEvolutionSuccessorKind::TRANSITION) {
                    possible_transition = true;
                    ARIADNE_TESTING_EXPECT_EQUAL(successor.resolution,
                                                 HybridEvolutionTransitionResolution::UNRESOLVED);
                }
            }
            ARIADNE_TESTING_PROVE(possible_transition);
        }
    }

    void test_step_preserves_transition_and_subdivision_provenance() {
        DiscreteLocation source(1), target(2);
        DiscreteEvent event("jump");
        RealVariable x("x"), y("y");
        HybridAutomaton system;
        system.new_mode(source, {dot(x) = 1, dot(y) = 0});
        system.new_mode(target, {dot(x) = 0, dot(y) = 0});
        system.new_transition(source, event, target, {next(x) = x, next(y) = y}, x >= 0, EventKind::URGENT);
        StepProbe probe(system);
        HybridTime horizon(2, 1);
        auto active = probe.enclosure(HybridExactBox(source, RealSpace({
                                                                 x, y
        }),
                                                     ExactBoxType({{1, 1}, {0, 0}})));
        auto jumped = probe._process_working_set({active, HybridEvolutionPhase::PROCESS_STARTING_EVENTS}, horizon, Semantics::UPPER, {});
        ARIADNE_TESTING_EXPECT_EQUAL(jumped.successors.size(), 1u);
        const auto& successor = jumped.successors.at(0);
        ARIADNE_TESTING_PROVE(successor.kind == HybridEvolutionSuccessorKind::TRANSITION);
        ARIADNE_TESTING_PROVE(successor.work.phase == HybridEvolutionPhase::PROCESS_STARTING_EVENTS);
        ARIADNE_TESTING_PROVE(successor.event == event);
        ARIADNE_TESTING_PROVE(successor.event_kind == EventKind::URGENT);
        ARIADNE_TESTING_PROVE(!successor.crossing.has_value());
        ARIADNE_TESTING_PROVE(successor.resolution == HybridEvolutionTransitionResolution::FULL_SOURCE);
        ARIADNE_TESTING_EXPECT_EQUAL(jumped.contributions.size(), 2u);
        ARIADNE_TESTING_PROVE(jumped.contributions.at(0).kind == HybridEvolutionSetKind::REACH);
        ARIADNE_TESTING_PROVE(jumped.contributions.at(0).enclosure.location() == source);
        ARIADNE_TESTING_PROVE(jumped.contributions.at(1).kind == HybridEvolutionSetKind::INTERMEDIATE);
        ARIADNE_TESTING_PROVE(jumped.contributions.at(1).enclosure.location() == target);

        probe.configuration().set_enable_subdivisions(true);
        probe.configuration().set_maximum_enclosure_radius(0.25_x);
        auto wide = probe.enclosure(HybridExactBox(source, RealSpace({
                                                               x, y
        }),
                                                   ExactBoxType({{-1, 1}, {0, 0}})));
        auto split = probe._process_working_set({wide, HybridEvolutionPhase::EVOLVE_CONTINUOUSLY}, horizon, Semantics::UPPER, {});
        ARIADNE_TESTING_EXPECT_EQUAL(split.continuous_evolution_time, 0.0);
        ARIADNE_TESTING_PROVE(split.successors.size() > 1u);
        ARIADNE_TESTING_PROVE(split.contributions.empty());
        for (const auto& piece: split.successors) {
            ARIADNE_TESTING_PROVE(piece.kind == HybridEvolutionSuccessorKind::SUBDIVISION);
            ARIADNE_TESTING_PROVE(piece.work.phase == HybridEvolutionPhase::EVOLVE_CONTINUOUSLY);
            ARIADNE_TESTING_PROVE(!piece.event.has_value());
        }
        auto partial_jump = probe._process_working_set({wide, HybridEvolutionPhase::PROCESS_STARTING_EVENTS}, horizon, Semantics::UPPER, {});
        ARIADNE_TESTING_PROVE(partial_jump.successors.size() > 1u);
        for (const auto& piece: partial_jump.successors) {
            ARIADNE_TESTING_PROVE(piece.resolution == HybridEvolutionTransitionResolution::UNRESOLVED);
        }
        auto discarded = probe._process_working_set({wide, HybridEvolutionPhase::EVOLVE_CONTINUOUSLY}, horizon, Semantics::LOWER, {});
        ARIADNE_TESTING_PROVE(discarded.successors.empty());
        ARIADNE_TESTING_PROVE(discarded.contributions.empty());
        ARIADNE_TESTING_PROVE(discarded.stop_reason == HybridEvolutionStopReason::LOWER_SIZE_LIMIT);
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_evolution_backend");
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_temporal_progress_uses_the_initial_time_and_terminal_branch_credit);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_invalid_cost_metadata_does_not_reset_valid_temporal_progress);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_progress_retains_the_120_second_prefix_when_permissive_branches_accumulate);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_progress_uses_lower_clock_bounds_for_source_and_replacement_estimates);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_progress_revises_remaining_work_without_recounting_history);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_progress_handles_zero_and_unknown_horizons);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_flat_and_tree_report_the_same_time_work_after_splitting);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_numerical_progress_uses_the_selected_continuous_duration);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_output_and_work_are_distinct);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_failed_step_cannot_release_an_orbit);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_evolver_accepts_an_independent_scheduler);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_lower_initial_activation_requires_the_whole_source);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_step_preserves_transition_and_subdivision_provenance);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
