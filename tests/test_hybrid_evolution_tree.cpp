/***************************************************************************
 *            test_hybrid_evolution_tree.cpp
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

/** @file tests/test_hybrid_evolution_tree.cpp
 *  @brief Verify trajectory history, parent-first scheduling and Orbit compatibility.
 */

#include <algorithm>
#include "concise_logging.hpp"
#include <memory>
#include <sstream>
#include <stdexcept>

#include "evolution/hybrid_evolution_tree.hpp"
#include "evolution/hybrid_evolver.hpp"
#include "system/hybrid_automata.hpp"
#include "threading/thread_manager.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

HybridEnclosure seed() {
    DiscreteLocation q(1);
    RealVariable x("x"), y("y");
    HybridAutomaton system;
    system.new_mode(q, {dot(x) = 1, dot(y) = 0});
    GeneralHybridEvolver evolver(system);
    return evolver.enclosure(HybridExactBox(q, RealSpace({
                                                   x, y
    }),
                                            ExactBoxType({{0, 0}, {0, 0}})));
}

// Clock values label synthetic work items in scheduler-only tests.
HybridEnclosure marked(const HybridEnclosure& initial, int label) {
    auto result = initial;
    result.set_time_function(initial.configuration().function_factory().create_constant(initial.parameter_domain(), Real(label)));
    return result;
}

int label(const HybridEvolutionWorkItem& item) {
    return static_cast<int>(item.enclosure.time_range().lower_bound().get_d());
}

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

/** @brief Sound no-change numerical strategy, used to verify injected recombiners.
 *  @details One output per input is valid but cannot reduce the number of
 *  trajectories; the backend should preserve the original pending children.
 */
class IdentityEnclosureRecombiner final : public EnclosureRecombinerInterface {
  public:
    IdentityEnclosureRecombiner* clone() const override { return new IdentityEnclosureRecombiner(*this); }

    EnclosureRecombinationResult recombine(
        const List<Enclosure>& enclosures, const EnclosureRecombinationContext&) const override {
        EnclosureRecombinationResult result;
        for (SizeType index = 0; index != enclosures.size(); ++index) {
            EnclosureRecombinationOutput output{enclosures[index], {}};
            output.source_indices.append(index);
            result.outputs.append(output);
        }
        return result;
    }

    OutputStream& _write(OutputStream& os) const override { return os << "IdentityEnclosureRecombiner"; }
};

class ConcurrencyScope {
    std::size_t _previous = ThreadManager::instance().concurrency();

  public:
    ConcurrencyScope() { ThreadManager::instance().set_concurrency(0); }
    ~ConcurrencyScope() { ThreadManager::instance().set_concurrency(_previous); }
};

} // namespace

class TestHybridEvolutionTree : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

  private:
    void compare_backend_orbits(GeneralHybridEvolver& evolver, const HybridEnclosure& initial,
                                const HybridTerminationCriterion& termination, Semantics semantics) {
        std::ostringstream original;
        original << initial;
        evolver.set_evolution_backend_factory(FlatHybridEvolutionBackendFactory());
        const auto expected = evolver.orbit(initial, termination, semantics);
        evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory());
        const auto actual = evolver.orbit(initial, termination, semantics);

        ARIADNE_TESTING_EXPECT_EQUAL(actual.reach().size(), expected.reach().size());
        ARIADNE_TESTING_EXPECT_EQUAL(actual.intermediate().size(), expected.intermediate().size());
        ARIADNE_TESTING_EXPECT_EQUAL(actual.final().size(), expected.final().size());
        ARIADNE_TESTING_PROVE(descriptions(actual.reach()) == descriptions(expected.reach()));
        ARIADNE_TESTING_PROVE(descriptions(actual.intermediate()) == descriptions(expected.intermediate()));
        ARIADNE_TESTING_PROVE(descriptions(actual.final()) == descriptions(expected.final()));
        std::ostringstream expected_initial, actual_initial;
        expected_initial << expected.initial();
        actual_initial << actual.initial();
        ARIADNE_TESTING_EXPECT_EQUAL(expected_initial.str(), original.str());
        ARIADNE_TESTING_EXPECT_EQUAL(actual_initial.str(), original.str());
    }

  public:
    /** @brief Compare complete numerical Orbit collections with the flat oracle.
     *  @details With recombination disabled, the frontier is permitted to change
     *  work order and genealogy, not reachability results or classifications.
     */
    void compare_flat_and_frontier_orbits(
        GeneralHybridEvolver& evolver, const HybridEnclosure& initial,
        const HybridTerminationCriterion& termination, Semantics semantics) {
        evolver.set_evolution_backend_factory(FlatHybridEvolutionBackendFactory());
        const auto flat = evolver.orbit(initial, termination, semantics);
        evolver.set_evolution_backend_factory(
            TreeHybridEvolutionBackendFactory(NoFrontierRecombinationPolicy{}));
        const auto frontier = evolver.orbit(initial, termination, semantics);

        ARIADNE_TESTING_EXPECT_EQUAL(frontier.reach().size(), flat.reach().size());
        ARIADNE_TESTING_EXPECT_EQUAL(frontier.intermediate().size(), flat.intermediate().size());
        ARIADNE_TESTING_EXPECT_EQUAL(frontier.final().size(), flat.final().size());
        ARIADNE_TESTING_PROVE(descriptions(frontier.reach()) == descriptions(flat.reach()));
        ARIADNE_TESTING_PROVE(descriptions(frontier.intermediate()) == descriptions(flat.intermediate()));
        ARIADNE_TESTING_PROVE(descriptions(frontier.final()) == descriptions(flat.final()));
        std::ostringstream expected_initial, actual_initial;
        expected_initial << flat.initial();
        actual_initial << frontier.initial();
        ARIADNE_TESTING_EXPECT_EQUAL(actual_initial.str(), expected_initial.str());
    }

  public:
    void test_flat_and_frontier_orbits_match_on_real_hybrid_transitions() {
        ConcurrencyScope concurrency;
        for (auto event_kind: {EventKind::URGENT, EventKind::PERMISSIVE}) {
            DiscreteLocation source(1), target(2);
            RealVariable x("x"), y("y");
            DiscreteEvent jump("jump");
            HybridAutomaton system;
            system.new_mode(source, {dot(x) = 1, dot(y) = 0});
            system.new_mode(target, {dot(x) = 0, dot(y) = 0});
            system.new_transition(source, jump, target,
                                  {next(x) = x, next(y) = y}, x >= 0, event_kind);
            GeneralHybridEvolver evolver(system);
            evolver.configuration().set_maximum_step_size(0.25_x);
            const auto initial = evolver.enclosure(HybridExactBox(
                source, RealSpace({x, y}),
                ExactBoxType({{-1, -0.75_x}, {0, 0}})));
            for (auto semantics: {Semantics::UPPER, Semantics::LOWER})
                compare_flat_and_frontier_orbits(evolver, initial, HybridTime(2, 1), semantics);
        }
    }

    void test_flat_and_frontier_orbits_match_on_real_splitting_and_invariants() {
        ConcurrencyScope concurrency;
        DiscreteLocation source(1), target(2);
        RealVariable x("x"), y("y");
        DiscreteEvent jump("jump"), invariant("invariant");
        HybridAutomaton system;
        system.new_mode(source, {dot(x) = 1, dot(y) = 0});
        system.new_mode(target, {dot(x) = 0, dot(y) = 0});
        system.new_transition(source, jump, target,
                              {next(x) = x, next(y) = y}, x >= 0, EventKind::URGENT);
        GeneralHybridEvolver evolver(system);
        evolver.configuration().set_enable_subdivisions(true);
        evolver.configuration().set_maximum_enclosure_radius(0.25_x);
        const auto wide = evolver.enclosure(HybridExactBox(
            source, RealSpace({x, y}), ExactBoxType({{-1, 1}, {0, 0}})));
        for (auto semantics: {Semantics::UPPER, Semantics::LOWER}) {
            compare_flat_and_frontier_orbits(evolver, wide, HybridTime(0.03125_x, 1), semantics);
            compare_flat_and_frontier_orbits(evolver, wide, HybridTime(0, 1), semantics);
        }
        system.new_invariant(source, x <= -0.5_dec, invariant);
        GeneralHybridEvolver constrained(system);
        constrained.configuration().set_enable_subdivisions(true);
        constrained.configuration().set_maximum_enclosure_radius(0.25_x);
        const auto restricted = constrained.enclosure(HybridExactBox(
            source, RealSpace({x, y}), ExactBoxType({{-1, 1}, {0, 0}})));
        for (auto semantics: {Semantics::UPPER, Semantics::LOWER})
            compare_flat_and_frontier_orbits(
                constrained, restricted, HybridTime(0.03125_x, 1), semantics);
    }

    void test_parent_finishes_before_distinct_derived_trajectories() {
        auto initial = seed();
        TreeHybridEvolutionBackend backend(initial);
        const auto root = backend.tree().root().value();
        std::vector<int> order;
        DiscreteEvent event("crossing");
        ARIADNE_TESTING_EXPECT_EXCEPTION(backend.release_orbit());
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            const auto marker = label(item);
            order.push_back(marker);
            HybridEvolutionStepResult result;
            if (marker < 2) {
                ARIADNE_TESTING_PROVE(!backend.tree().metadata(root).complete);
                if (marker == 1) {
                    auto first = backend.tree().children(root).front();
                    ARIADNE_TESTING_PROVE(!backend.tree().metadata(first).complete);
                    ARIADNE_TESTING_PROVE(!backend.tree().object(first).entries.front().result.has_value());
                }
                result.transition_to(marked(initial, marker == 0 ? 100 : 200), event, EventKind::URGENT);
                result.continue_with(marked(initial, marker + 1));
            } else if (marker == 2) {
                result.record_final(item.enclosure);
                result.transition_to(marked(initial, 300), event, EventKind::URGENT, std::nullopt, HybridEvolutionTransitionResolution::FULL_SOURCE);
            } else {
                ARIADNE_TESTING_PROVE(backend.tree().metadata(root).complete);
                ARIADNE_TESTING_EXPECT_EQUAL(backend.tree().children(root).size(), 3u);
                if (marker == 200)
                    result.record_final(item.enclosure);
                else
                    result.stop_reason = HybridEvolutionStopReason::EMPTY;
            }
            return result;
        });
        ARIADNE_TESTING_PROVE(order == std::vector<int>({0, 1, 2, 100, 200, 300}));
        const auto& tree = backend.tree();
        ARIADNE_TESTING_EXPECT_EQUAL(tree.size(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(tree.object(root).entries.size(), 3u);
        const auto children = tree.children(root);
        ARIADNE_TESTING_PROVE(children[0] != children[1]);
        for (std::size_t index = 0; index != children.size(); ++index) {
            ARIADNE_TESTING_PROVE(tree.parent(children[index]) == root);
            ARIADNE_TESTING_PROVE(tree.metadata(children[index]).complete);
            const auto origin = tree.metadata(children[index]).origin.value();
            ARIADNE_TESTING_EXPECT_EQUAL(origin.entry, index);
            const auto& source = tree.object(root).entries.at(origin.entry).result.value().successors.at(origin.successor);
            ARIADNE_TESTING_PROVE(source.event == event);
        }
        auto orbit = backend.release_orbit();
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.reach().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.initial().time_range().lower_bound().get_d(), 0.0);
        ARIADNE_TESTING_PROVE(tree.contains(children[0]));
        ARIADNE_TESTING_EXPECT_EXCEPTION(backend.release_orbit());
        ARIADNE_TESTING_EXPECT_EXCEPTION(backend.run({}));
    }

    void test_only_certified_single_urgent_transition_continues_the_node() {
        auto initial = seed();
        DiscreteEvent event("jump");
        for (auto kind: {EventKind::URGENT, EventKind::PERMISSIVE}) {
            for (auto resolution: {HybridEvolutionTransitionResolution::UNRESOLVED, HybridEvolutionTransitionResolution::FULL_SOURCE}) {
                TreeHybridEvolutionBackend backend(initial);
                const auto root = backend.tree().root().value();
                backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
                    HybridEvolutionStepResult result;
                    if (label(item) == 0)
                        result.transition_to(marked(initial, 100), event, kind, std::nullopt, resolution);
                    else
                        result.record_final(item.enclosure);
                    return result;
                });
                const bool continues = kind == EventKind::URGENT && resolution == HybridEvolutionTransitionResolution::FULL_SOURCE;
                ARIADNE_TESTING_EXPECT_EQUAL(backend.tree().size(), continues ? 1u : 2u);
                ARIADNE_TESTING_EXPECT_EQUAL(backend.tree().object(root).entries.size(), continues ? 2u : 1u);
                if (continues) {
                    const auto origin = backend.tree().object(root).entries.at(1).origin.value();
                    ARIADNE_TESTING_EXPECT_EQUAL(origin.entry, 0u);
                    ARIADNE_TESTING_EXPECT_EQUAL(origin.successor, 0u);
                }
                ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 1u);
            }
        }
    }

    void test_noop_and_window_recombination_policies() {
        const auto initial = seed();
        const DiscreteEvent jump("jump");
        const DiscreteEvent other("other");
        const auto advance = [&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            if (label(item) == 0) {
                result.transition_to(marked(initial, 100), jump, EventKind::URGENT);
                result.transition_to(marked(initial, 100), jump, EventKind::URGENT);
                result.transition_to(marked(initial, 200), other, EventKind::URGENT);
            } else {
                result.record_final(item.enclosure);
            }
            return result;
        };

        NoRecombinationSchedulingPolicy noop;
        TreeHybridEvolutionBackend baseline(initial, noop);
        baseline.set_semantics(Semantics::UPPER);
        baseline.run(advance);
        ARIADNE_TESTING_EXPECT_EQUAL(baseline.tree().size(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(baseline.release_orbit().final().size(), 3u);

        SiblingWindowRecombinationSchedulingPolicy window(2u, 0.0);
        TreeHybridEvolutionBackend recombined(initial, window);
        const auto root = recombined.tree().root().value();
        recombined.set_semantics(Semantics::UPPER);
        recombined.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics& stats) {
            if (label(item) != 0)
                ARIADNE_TESTING_PROVE(recombined.tree().metadata(root).complete);
            return advance(item, stats);
        });
        const auto children = recombined.tree().children(root);
        ARIADNE_TESTING_EXPECT_EQUAL(children.size(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(recombined.tree().size(), 5u);
        ARIADNE_TESTING_PROVE(recombined.tree().metadata(children[0]).retired_for_recombination);
        ARIADNE_TESTING_PROVE(recombined.tree().metadata(children[1]).retired_for_recombination);
        ARIADNE_TESTING_PROVE(!recombined.tree().object(children[0]).entries[0].result);
        ARIADNE_TESTING_PROVE(!recombined.tree().object(children[1]).entries[0].result);
        ARIADNE_TESTING_PROVE(!recombined.tree().metadata(children[2]).retired_for_recombination);
        const auto& replacement = recombined.tree().metadata(children[3]);
        ARIADNE_TESTING_PROVE(replacement.complete);
        ARIADNE_TESTING_PROVE(!replacement.origin);
        ARIADNE_TESTING_EXPECT_EQUAL(replacement.recombined_sources.size(), 2u);
        ARIADNE_TESTING_PROVE(replacement.recombined_sources[0] == children[0]);
        ARIADNE_TESTING_PROVE(replacement.recombined_sources[1] == children[1]);
        const auto orbit = recombined.release_orbit();
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.reach().size(), 0u);
    }

    void test_recombination_policy_rejects_lower_semantics_and_invalid_settings() {
        const auto initial = seed();
        ARIADNE_TESTING_EXPECT_EXCEPTION(SiblingWindowRecombinationSchedulingPolicy(1u, 1.0));
        ARIADNE_TESTING_EXPECT_EXCEPTION(SiblingWindowRecombinationSchedulingPolicy(2u, -1.0));
        TreeHybridEvolutionBackend lower(initial, SiblingWindowRecombinationSchedulingPolicy{});
        lower.set_semantics(Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EXCEPTION(lower.run([](const HybridEvolutionWorkItem&, const HybridEvolutionStatistics&) {
            return HybridEvolutionStepResult{};
        }));
        ARIADNE_TESTING_EXPECT_EQUAL(lower.tree().size(), 1u);
        ARIADNE_TESTING_PROVE(!lower.tree().metadata(lower.tree().root().value()).complete);
    }

    void test_recombination_evaluation_triggers() {
        const auto initial = seed();
        const DiscreteEvent jump("jump");
        const auto verify = [&](const HybridRecombinationSchedulingPolicy& grouping,
                                const HybridRecombinationTriggerPolicy& timing,
                                bool after_first, bool after_second, bool after_complete,
                                std::size_t expected_final_count) {
            TreeHybridEvolutionBackend backend(initial, grouping, timing);
            const auto root = backend.tree().root().value();
            bool observed_after_first = false;
            bool observed_after_second = false;
            backend.set_semantics(Semantics::UPPER);
            backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
                const auto marker = label(item);
                HybridEvolutionStepResult result;
                if (marker == 0) {
                    // Two discrete successors remain blocked while the parent continues.
                    result.transition_to(marked(initial, 100), jump, EventKind::URGENT);
                    result.transition_to(marked(initial, 100), jump, EventKind::URGENT);
                    result.continue_with(marked(initial, 1));
                } else if (marker == 1 || marker == 2) {
                    const auto children = backend.tree().children(root);
                    ARIADNE_TESTING_PROVE(!backend.tree().metadata(root).complete);
                    ARIADNE_TESTING_PROVE(children.size() >= 2u);
                    const bool retired = backend.tree().metadata(children[0]).retired_for_recombination;
                    if (marker == 1)
                        observed_after_first = retired;
                    else
                        observed_after_second = retired;
                    if (marker == 1)
                        result.continue_with(marked(initial, 2));
                    else
                        result.record_final(item.enclosure);
                } else {
                    ARIADNE_TESTING_PROVE(backend.tree().metadata(root).complete);
                    result.record_final(item.enclosure);
                }
                return result;
            });
            const auto children = backend.tree().children(root);
            ARIADNE_TESTING_EXPECT_EQUAL(observed_after_first, after_first);
            ARIADNE_TESTING_EXPECT_EQUAL(observed_after_second, after_second);
            ARIADNE_TESTING_EXPECT_EQUAL(backend.tree().metadata(children[0]).retired_for_recombination,
                                         after_complete);
            ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), expected_final_count);
        };
        const SiblingWindowRecombinationSchedulingPolicy window(2u, 0.0);
        verify(NoRecombinationSchedulingPolicy{}, AfterEachStepRecombinationTriggerPolicy{},
               false, false, false, 3u);
        verify(window, OnParentCompletionRecombinationTriggerPolicy{},
               false, false, true, 2u);
        verify(window, AfterEachStepRecombinationTriggerPolicy{},
               true, true, true, 2u);
        verify(window, EveryNStepsRecombinationTriggerPolicy(2u),
               false, true, true, 2u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(EveryNStepsRecombinationTriggerPolicy(0u));
    }

    void test_numerical_recombination_strategy_is_injected() {
        const auto initial = seed();
        const DiscreteEvent jump("jump");
        TreeHybridEvolutionBackendFactory factory(
            SiblingWindowRecombinationSchedulingPolicy(2u, 0.0),
            AfterEachStepRecombinationTriggerPolicy{},
            IdentityEnclosureRecombiner{});
        // The factory owns clones of all three temporary strategies.
        const auto advance = [&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            if (label(item) == 0) {
                result.transition_to(marked(initial, 100), jump, EventKind::URGENT);
                result.transition_to(marked(initial, 100), jump, EventKind::URGENT);
            } else {
                result.record_final(item.enclosure);
            }
            return result;
        };
        auto first = factory.create(initial, Real(1));
        first->set_semantics(Semantics::UPPER);
        first->run(advance);
        // A sound identity strategy supplies two outputs for two sources;
        // the scheduler correctly rejects this ineffective replacement.
        ARIADNE_TESTING_EXPECT_EQUAL(first->release_orbit().final().size(), 2u);
        auto second = factory.create(initial, Real(1));
        second->set_semantics(Semantics::UPPER);
        second->run(advance);
        ARIADNE_TESTING_EXPECT_EQUAL(second->release_orbit().final().size(), 2u);
    }

    void test_frontier_recombines_continuous_split_branches_after_evolution() {
        const auto initial = seed();
        const auto advance = [&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            switch (label(item)) {
                case 0:
                    result.continue_with(marked(initial, 10), HybridEvolutionSuccessorKind::SUBDIVISION);
                    result.continue_with(marked(initial, 20), HybridEvolutionSuccessorKind::SUBDIVISION);
                    break;
                case 10:
                case 20:
                    result.continue_with(marked(initial, 30));
                    break;
                case 30:
                    result.record_final(item.enclosure);
                    break;
                default:
                    throw std::logic_error("Unexpected frontier checkpoint");
            }
            return result;
        };
        // The no-op frontier still interleaves numerical pieces without merging.
        TreeHybridEvolutionBackend baseline(initial, NoFrontierRecombinationPolicy{});
        std::vector<int> baseline_order;
        baseline.run([&](const HybridEvolutionWorkItem& work, const HybridEvolutionStatistics& stats) {
            baseline_order.push_back(label(work));
            return advance(work, stats);
        });
        ARIADNE_TESTING_PROVE(baseline_order == std::vector<int>({0, 10, 20, 30, 30}));
        ARIADNE_TESTING_EXPECT_EQUAL(baseline.release_orbit().final().size(), 2u);

        // A three-step trigger waits until both split branches have evolved
        // independently once; neither branch is complete at the merge point.
        TreeHybridEvolutionBackend merged(
            initial, FrontierWindowRecombinationPolicy(2u, 0.0),
            EveryNStepsRecombinationTriggerPolicy(3u));
        merged.set_semantics(Semantics::UPPER);
        std::vector<int> merged_order;
        merged.run([&](const HybridEvolutionWorkItem& work, const HybridEvolutionStatistics& stats) {
            merged_order.push_back(label(work));
            return advance(work, stats);
        });
        ARIADNE_TESTING_PROVE(merged_order == std::vector<int>({0, 10, 20, 30}));
        const auto root = merged.tree().root().value();
        const auto children = merged.tree().children(root);
        ARIADNE_TESTING_EXPECT_EQUAL(children.size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(merged.tree().object(children[0]).retired_entries.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(merged.tree().object(children[1]).retired_entries.size(), 1u);
        ARIADNE_TESTING_PROVE(merged.tree().metadata(children[0]).complete);
        ARIADNE_TESTING_PROVE(merged.tree().metadata(children[1]).complete);
        const auto& provenance = merged.tree().metadata(children[2]).recombined_checkpoints;
        ARIADNE_TESTING_EXPECT_EQUAL(provenance.size(), 2u);
        const HybridEvolutionCheckpointId first_expected{children[0], 1u};
        const HybridEvolutionCheckpointId second_expected{children[1], 1u};
        ARIADNE_TESTING_PROVE(provenance[0] == first_expected);
        ARIADNE_TESTING_PROVE(provenance[1] == second_expected);
        ARIADNE_TESTING_EXPECT_EQUAL(merged.release_orbit().final().size(), 1u);
    }

    void test_frontier_recombines_cross_trajectory_transitions() {
        const auto initial = seed();
        const DiscreteEvent jump("jump");
        for (auto event_kind: {EventKind::PERMISSIVE, EventKind::URGENT}) {
            const auto advance = [&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
                HybridEvolutionStepResult result;
                switch (label(item)) {
                    case 0:
                        result.continue_with(marked(initial, 10), HybridEvolutionSuccessorKind::SUBDIVISION);
                        result.continue_with(marked(initial, 20), HybridEvolutionSuccessorKind::SUBDIVISION);
                        break;
                    case 10:
                    case 20:
                        result.transition_to(marked(initial, 100), jump, event_kind);
                        break;
                    case 100:
                        result.record_final(item.enclosure);
                        break;
                    default:
                        throw std::logic_error("Unexpected transition checkpoint");
                }
                return result;
            };
            TreeHybridEvolutionBackend baseline(initial, NoFrontierRecombinationPolicy{});
            baseline.run(advance);
            ARIADNE_TESTING_EXPECT_EQUAL(baseline.release_orbit().final().size(), 2u);

            TreeHybridEvolutionBackend merged(
                initial, FrontierWindowRecombinationPolicy(2u, 0.0),
                EveryNStepsRecombinationTriggerPolicy(3u));
            merged.set_semantics(Semantics::UPPER);
            merged.run(advance);
            const auto root = merged.tree().root().value();
            const auto children = merged.tree().children(root);
            ARIADNE_TESTING_EXPECT_EQUAL(children.size(), 3u);
            const auto& sources = merged.tree().metadata(children[2]).recombined_checkpoints;
            ARIADNE_TESTING_EXPECT_EQUAL(sources.size(), 2u);
            ARIADNE_TESTING_PROVE(merged.tree().parent(sources[0].trajectory) == children[0]);
            ARIADNE_TESTING_PROVE(merged.tree().parent(sources[1].trajectory) == children[1]);
            ARIADNE_TESTING_EXPECT_EQUAL(merged.release_orbit().final().size(), 1u);
        }
        ARIADNE_TESTING_EXPECT_EXCEPTION(FrontierWindowRecombinationPolicy(1u, 0.125));
        ARIADNE_TESTING_EXPECT_EXCEPTION(FrontierWindowRecombinationPolicy(2u, -1.0));
        TreeHybridEvolutionBackend lower(initial, FrontierWindowRecombinationPolicy{});
        lower.set_semantics(Semantics::LOWER);
        ARIADNE_TESTING_EXPECT_EXCEPTION(lower.run([](const HybridEvolutionWorkItem&, const HybridEvolutionStatistics&) {
            return HybridEvolutionStepResult{};
        }));
    }

    void test_split_pieces_are_parallel_and_completion_waits_for_all_of_them() {
        auto initial = seed();
        TreeHybridEvolutionBackend backend(initial);
        const auto root = backend.tree().root().value();
        std::vector<int> order;
        int visits_to_first_piece = 0;
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            const auto marker = label(item);
            order.push_back(marker);
            HybridEvolutionStepResult result;
            if (marker == 0) {
                result.continue_with(marked(initial, 1), HybridEvolutionSuccessorKind::SUBDIVISION);
                result.continue_with(marked(initial, 2), HybridEvolutionSuccessorKind::SUBDIVISION);
            } else if (marker == 1 && ++visits_to_first_piece == 1) {
                result.continue_with(item.enclosure, HybridEvolutionSuccessorKind::RECONDITIONING);
            } else if (marker == 2) {
                result.transition_to(marked(initial, 100), DiscreteEvent("jump"), EventKind::URGENT, std::nullopt, HybridEvolutionTransitionResolution::FULL_SOURCE);
            } else if (marker == 100) {
                ARIADNE_TESTING_EXPECT_EQUAL(visits_to_first_piece, 2);
                ARIADNE_TESTING_PROVE(backend.tree().metadata(root).complete);
            }
            return result;
        });
        ARIADNE_TESTING_PROVE(order == std::vector<int>({0, 1, 2, 1, 100}));
        const auto& entries = backend.tree().object(root).entries;
        ARIADNE_TESTING_EXPECT_EQUAL(entries.size(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(entries.at(1).origin->entry, 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(entries.at(2).origin->entry, 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(entries.at(1).origin->successor, 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(entries.at(2).origin->successor, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(entries.at(3).origin->entry, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(backend.tree().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 0u);
    }

    void test_failure_preserves_incomplete_history_without_dispatching_children() {
        auto initial = seed();
        TreeHybridEvolutionBackend backend(initial);
        std::vector<int> order;
        ARIADNE_TESTING_EXPECT_EXCEPTION(backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            order.push_back(label(item));
            if (label(item) != 0)
                throw std::runtime_error("Failed continuation");
            HybridEvolutionStepResult result;
            result.transition_to(marked(initial, 100), DiscreteEvent("jump"), EventKind::URGENT);
            result.continue_with(marked(initial, 1));
            return result;
        }));
        const auto& tree = backend.tree();
        const auto root = tree.root().value();
        ARIADNE_TESTING_PROVE(order == std::vector<int>({0, 1}));
        ARIADNE_TESTING_PROVE(!tree.metadata(root).complete);
        ARIADNE_TESTING_PROVE(tree.object(root).entries.at(0).result.has_value());
        ARIADNE_TESTING_PROVE(!tree.object(root).entries.at(1).result.has_value());
        const auto child = tree.children(root).front();
        ARIADNE_TESTING_PROVE(!tree.metadata(child).complete);
        ARIADNE_TESTING_PROVE(!tree.object(child).entries.front().result.has_value());
        ARIADNE_TESTING_EXPECT_EXCEPTION(backend.release_orbit());
    }

    void test_flat_and_tree_orbits_match_for_urgent_and_permissive_systems() {
        ConcurrencyScope concurrency;
        for (auto kind: {EventKind::URGENT, EventKind::PERMISSIVE}) {
            DiscreteLocation source(1), target(2);
            RealVariable x("x"), y("y");
            HybridAutomaton system;
            system.new_mode(source, {dot(x) = 1, dot(y) = 0});
            system.new_mode(target, {dot(x) = 0, dot(y) = 0});
            system.new_transition(source, DiscreteEvent("jump"), target, {next(x) = x, next(y) = y}, x >= 0, kind);
            GeneralHybridEvolver evolver(system);
            evolver.configuration().set_maximum_step_size(0.25_x);
            auto initial = evolver.enclosure(HybridExactBox(source, RealSpace({
                                                                        x, y
            }),
                                                            ExactBoxType({{-1, -0.75_x}, {0, 0}})));
            for (auto semantics: {Semantics::UPPER, Semantics::LOWER}) {
                compare_backend_orbits(evolver, initial, HybridTime(2, 1), semantics);
            }
        }
    }

    void test_flat_and_tree_orbit_counts_match_at_time_and_discrete_limits() {
        ConcurrencyScope concurrency;
        DiscreteLocation q(1);
        RealVariable x("x"), y("y");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0});
        GeneralHybridEvolver evolver(system);
        evolver.configuration().set_maximum_step_size(0.25_x);
        const auto initial = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                                     x, y
        }),
                                                              ExactBoxType({{-1, -0.75_x}, {0, 0}})));
        for (auto semantics: {Semantics::UPPER, Semantics::LOWER}) {
            compare_backend_orbits(evolver, initial, HybridTime(0, 1), semantics);
            compare_backend_orbits(evolver, initial, HybridTime(2, 0), semantics);
            compare_backend_orbits(evolver, initial, HybridTime(2, 1), semantics);
        }
    }

    void test_flat_and_tree_orbit_counts_match_with_initial_jumps_and_splitting() {
        ConcurrencyScope concurrency;
        DiscreteLocation q(1), target(2);
        DiscreteEvent jump("jump"), invariant("invariant");
        RealVariable x("x"), y("y");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0});
        system.new_mode(target, {dot(x) = 0, dot(y) = 0});
        system.new_transition(q, jump, target, {next(x) = x, next(y) = y}, x >= 0, EventKind::URGENT);
        GeneralHybridEvolver evolver(system);
        evolver.configuration().set_enable_subdivisions(true);
        evolver.configuration().set_maximum_enclosure_radius(0.25_x);
        const auto active = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                                    x, y
        }),
                                                             ExactBoxType({{1, 1}, {0, 0}})));
        for (auto semantics: {Semantics::UPPER, Semantics::LOWER}) {
            compare_backend_orbits(evolver, active, HybridTime(1, 1), semantics);
            compare_backend_orbits(evolver, active, HybridTerminationCriterion(Real(1), Integer(5), Set<DiscreteEvent>{jump}), semantics);
        }
        const auto wide = evolver.enclosure(HybridExactBox(q, RealSpace({
                                                                  x, y
        }),
                                                           ExactBoxType({{-1, 1}, {0, 0}})));
        for (auto semantics: {Semantics::UPPER, Semantics::LOWER}) {
            compare_backend_orbits(evolver, wide, HybridTime(0.03125_x, 1), semantics);
        }

        system.new_invariant(q, x <= -0.5_dec, invariant);
        GeneralHybridEvolver constrained(system);
        constrained.configuration().set_enable_subdivisions(true);
        constrained.configuration().set_maximum_enclosure_radius(0.25_x);
        for (auto semantics: {Semantics::UPPER, Semantics::LOWER}) {
            const auto restricted = constrained.enclosure(HybridExactBox(q, RealSpace({
                                                                                x, y
            }),
                                                                         ExactBoxType({{-1, 1}, {0, 0}})));
            compare_backend_orbits(constrained, restricted, HybridTime(0.03125_x, 1), semantics);
            const auto excluded = constrained.enclosure(HybridExactBox(q, RealSpace({
                                                                              x, y
            }),
                                                                       ExactBoxType({{1, 2}, {0, 0}})));
            compare_backend_orbits(constrained, excluded, HybridTime(0.03125_x, 1), semantics);
        }
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_evolution_tree");
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_parent_finishes_before_distinct_derived_trajectories);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_only_certified_single_urgent_transition_continues_the_node);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_noop_and_window_recombination_policies);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_recombination_policy_rejects_lower_semantics_and_invalid_settings);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_recombination_evaluation_triggers);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_numerical_recombination_strategy_is_injected);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_frontier_recombines_continuous_split_branches_after_evolution);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_frontier_recombines_cross_trajectory_transitions);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_flat_and_frontier_orbits_match_on_real_hybrid_transitions);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_flat_and_frontier_orbits_match_on_real_splitting_and_invariants);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_split_pieces_are_parallel_and_completion_waits_for_all_of_them);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_failure_preserves_incomplete_history_without_dispatching_children);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_flat_and_tree_orbits_match_for_urgent_and_permissive_systems);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_flat_and_tree_orbit_counts_match_at_time_and_discrete_limits);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_flat_and_tree_orbit_counts_match_with_initial_jumps_and_splitting);
    return run_evolution_test_suite(ariadne_test_suite, argc, argv);
}
