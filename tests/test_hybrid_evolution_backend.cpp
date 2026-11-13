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
#include <sstream>
#include <stdexcept>

#include "evolution/hybrid_evolution_backend.hpp"
#include "evolution/hybrid_evolver.hpp"
#include "system/hybrid_automata.hpp"
#include "test_logging.hpp"
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

class StepProbe : public GeneralHybridEvolver {
  public:
    using GeneralHybridEvolver::GeneralHybridEvolver;
    using HybridEvolverBase::_process_working_set;
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

    void test_output_and_work_are_distinct() {
        DiscreteLocation q(1);
        RealVariable x("x"), y("y");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0});
        GeneralHybridEvolver evolver(system);
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
        GeneralHybridEvolver evolver(system);
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
            GeneralHybridEvolver evolver(system);
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
        ARIADNE_TESTING_PROVE(!successor.crossing_kind.has_value());
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
    redirect_test_logging_to_null();
    ARIADNE_TESTING_SUITE("test_hybrid_evolution_backend");
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_output_and_work_are_distinct);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_failed_step_cannot_release_an_orbit);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_evolver_accepts_an_independent_scheduler);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionBackend, test_step_preserves_transition_and_subdivision_provenance);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
