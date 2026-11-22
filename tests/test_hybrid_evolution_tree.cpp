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
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>

#include "evolution/hybrid_continuous_time_evolver.hpp"
#include "evolution/hybrid_evolution_tree.hpp"
#include "logging/test_runner.hpp"
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
    HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
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

HybridEvolutionCheckpointId add_pending_checkpoint(
    HybridEvolutionTree& tree, const HybridEnclosure& enclosure,
    std::optional<DerivationTreeNodeId> parent = std::nullopt) {
    HybridEvolutionTrajectory trajectory;
    trajectory.entries.push_back({
        {enclosure, HybridEvolutionPhase::PROCESS_STARTING_EVENTS},
        std::nullopt,
        std::nullopt
    });
    const auto node = tree.empty()
                        ? tree.create_root(std::move(trajectory), HybridEvolutionTrajectoryMetadata{})
                        : tree.add_child(parent.value_or(tree.root().value()), std::move(trajectory), HybridEvolutionTrajectoryMetadata{});
    return {node, 0u};
}

HybridEnclosure jumped(const HybridEnclosure& source, DiscreteLocation target,
                       DiscreteEvent event, int time_label) {
    auto result = marked(source, time_label);
    result.apply_reset(event, target, source.state_space(),
                       EffectiveVectorMultivariateFunction::identity(source.state_dimension()),
                       source.auxiliary_space(), source.continuous_set().auxiliary_mapping());
    return result;
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

/** @brief Inspect a complete real-evolver frontier run without changing the scheduler. */
struct FrontierIntegrationAudit {
    std::size_t selection_checks = 0u;
    std::size_t selected_groups = 0u;
    std::size_t numerical_calls = 0u;
    std::size_t replacement_trajectories = 0u;
    std::size_t retired_checkpoints = 0u;
    std::size_t processed_checkpoints = 0u;
    std::size_t pending_checkpoints = 0u;
    std::size_t split_successors = 0u;
    std::size_t transition_successors = 0u;
    bool valid_provenance = true;
    bool all_trajectories_complete = true;
};

/** @brief Bound test-only recombination to one compatible window.
 *  @details This prevents a merged wide enclosure from immediately being
 *  recombined again if the physical evolver decides to subdivide it.
 */
class SingleWindowFrontierPolicy final : public HybridRecombinationGroupingPolicy {
    std::shared_ptr<FrontierIntegrationAudit> _audit;
    WindowRecombinationPolicy _window;
    bool _require_evolved_checkpoint;
    bool _require_reset_checkpoint;
    mutable bool _selected = false;

  public:
    explicit SingleWindowFrontierPolicy(std::shared_ptr<FrontierIntegrationAudit> audit,
                                        bool require_evolved_checkpoint = false,
                                        bool require_reset_checkpoint = false)
        : _audit(std::move(audit)), _window(2u, 0.5),
          _require_evolved_checkpoint(require_evolved_checkpoint),
          _require_reset_checkpoint(require_reset_checkpoint) {}

    SingleWindowFrontierPolicy* clone() const override {
        return new SingleWindowFrontierPolicy(*this);
    }

    std::vector<std::vector<CheckpointId>> select(
        const HybridEvolutionTree& tree, const std::vector<CheckpointId>& pending) const override {
        ++_audit->selection_checks;
        if (_selected)
            return {};
        std::vector<CheckpointId> eligible;
        for (const auto checkpoint: pending) {
            const auto& work = tree.object(checkpoint.trajectory).entries.at(checkpoint.entry).work;
            const double time_lower = work.enclosure.time_range().lower_bound().raw().get_d();
            if ((!_require_evolved_checkpoint || (checkpoint.entry > 0u && time_lower > 0.0)) &&
                (!_require_reset_checkpoint || !work.enclosure.previous_events().empty()))
                eligible.push_back(checkpoint);
        }
        const auto groups = _window.select(tree, eligible);
        if (groups.empty())
            return {};
        _selected = true;
        ++_audit->selected_groups;
        return {groups.front()};
    }

    bool enables_recombination() const override { return true; }
};

/** @brief Count actual numerical invocations while using the production algorithm. */
class ObservedReferenceEnclosureRecombiner final : public EnclosureRecombinerInterface {
    std::shared_ptr<FrontierIntegrationAudit> _audit;

  public:
    explicit ObservedReferenceEnclosureRecombiner(std::shared_ptr<FrontierIntegrationAudit> audit)
        : _audit(std::move(audit)) {}
    ObservedReferenceEnclosureRecombiner* clone() const override {
        return new ObservedReferenceEnclosureRecombiner(*this);
    }
    EnclosureRecombinationResult recombine(
        const List<Enclosure>& sources, const EnclosureRecombinationContext& context) const override {
        ++_audit->numerical_calls;
        return ReferenceEnclosureRecombiner().recombine(sources, context);
    }
    OutputStream& _write(OutputStream& os) const override {
        return os << "ObservedReferenceEnclosureRecombiner";
    }
};

/** @brief Expose a read-only genealogical audit through an evolver backend factory. */
class AuditedFrontierBackend final : public HybridEvolutionBackendInterface {
    TreeHybridEvolutionBackend _backend;
    std::shared_ptr<FrontierIntegrationAudit> _audit;

  public:
    AuditedFrontierBackend(const HybridEnclosure& initial,
                           const HybridCheckpointAdmissionPolicy& admission,
                           const HybridRecombinationGroupingPolicy& grouping,
                           const EnclosureRecombinationContext& context,
                           std::shared_ptr<FrontierIntegrationAudit> audit)
        : _backend(initial, admission, grouping, BeforeAdmissionRecombinationTriggerPolicy{},
                   ObservedReferenceEnclosureRecombiner(audit), context),
          _audit(std::move(audit)) {}
    AuditedFrontierBackend(const HybridEnclosure& initial,
                           const HybridRecombinationGroupingPolicy& grouping,
                           const EnclosureRecombinationContext& context,
                           std::shared_ptr<FrontierIntegrationAudit> audit)
        : _backend(initial, ImmediateCheckpointAdmissionPolicy{}, grouping, AfterEachStepRecombinationTriggerPolicy{},
                   ObservedReferenceEnclosureRecombiner(audit), context),
          _audit(std::move(audit)) {}

    Void set_semantics(Semantics semantics) override { _backend.set_semantics(semantics); }
    Void run(const StepFunction& step) override { _backend.run(step); }

    Orbit<HybridEnclosure> release_orbit() override {
        auto orbit = _backend.release_orbit();
        const auto& tree = _backend.tree();
        const auto root = tree.root().value();
        std::vector<HybridEvolutionTree::NodeId> nodes{root};
        const auto descendants = tree.descendants(root);
        nodes.insert(nodes.end(), descendants.begin(), descendants.end());
        for (const auto node: nodes) {
            const auto& trajectory = tree.object(node);
            const auto& metadata = tree.metadata(node);
            _audit->all_trajectories_complete =
                _audit->all_trajectories_complete && static_cast<bool>(metadata.complete);
            for (std::size_t entry = 0u; entry != trajectory.entries.size(); ++entry) {
                const bool retired = std::find(trajectory.retired_entries.begin(),
                                               trajectory.retired_entries.end(), entry) !=
                                     trajectory.retired_entries.end();
                if (retired)
                    ++_audit->retired_checkpoints;
                else if (trajectory.entries[entry].result)
                    ++_audit->processed_checkpoints;
                else
                    ++_audit->pending_checkpoints;
                if (trajectory.entries[entry].result) {
                    for (const auto& successor: trajectory.entries[entry].result->successors) {
                        if (successor.kind == HybridEvolutionSuccessorKind::SUBDIVISION)
                            ++_audit->split_successors;
                        if (successor.kind == HybridEvolutionSuccessorKind::TRANSITION)
                            ++_audit->transition_successors;
                    }
                }
            }
            if (!metadata.recombined_checkpoints.empty()) {
                ++_audit->replacement_trajectories;
                for (const auto source: metadata.recombined_checkpoints) {
                    if (!tree.contains(source.trajectory)) {
                        _audit->valid_provenance = false;
                        continue;
                    }
                    const auto& original = tree.object(source.trajectory);
                    if (source.entry >= original.entries.size() ||
                        original.entries[source.entry].result ||
                        std::find(original.retired_entries.begin(), original.retired_entries.end(),
                                  source.entry) == original.retired_entries.end())
                        _audit->valid_provenance = false;
                }
            }
        }
        return orbit;
    }
};

class AuditedFrontierBackendFactory final : public HybridEvolutionBackendFactoryInterface {
    std::shared_ptr<const HybridCheckpointAdmissionPolicy> _admission;
    std::shared_ptr<const HybridRecombinationGroupingPolicy> _grouping;
    EnclosureRecombinationContext _context;
    std::shared_ptr<FrontierIntegrationAudit> _audit;

  public:
    AuditedFrontierBackendFactory(const HybridCheckpointAdmissionPolicy& admission,
                                  const HybridRecombinationGroupingPolicy& grouping,
                                  const EnclosureRecombinationContext& context,
                                  std::shared_ptr<FrontierIntegrationAudit> audit)
        : _admission(admission.clone()), _grouping(grouping.clone()), _context(context), _audit(std::move(audit)) {}
    AuditedFrontierBackendFactory(const HybridRecombinationGroupingPolicy& grouping,
                                  const EnclosureRecombinationContext& context,
                                  std::shared_ptr<FrontierIntegrationAudit> audit)
        : _grouping(grouping.clone()), _context(context), _audit(std::move(audit)) {}
    AuditedFrontierBackendFactory* clone() const override {
        return new AuditedFrontierBackendFactory(*this);
    }
    std::unique_ptr<HybridEvolutionBackendInterface> create(
        const HybridEnclosure& initial, const Real&) const override {
        if (_admission)
            return std::make_unique<AuditedFrontierBackend>(initial, *_admission, *_grouping, _context, _audit);
        return std::make_unique<AuditedFrontierBackend>(initial, *_grouping, _context, _audit);
    }
};

class ConcurrencyScope {
    std::size_t _previous = ThreadManager::instance().concurrency();

  public:
    ConcurrencyScope() { ThreadManager::instance().set_concurrency(0); }
    ~ConcurrencyScope() { ThreadManager::instance().set_concurrency(_previous); }
};

/** @brief Restore the thread pool when a parallel regression finishes. */
class ParallelConcurrencyScope {
    std::size_t _previous = ThreadManager::instance().concurrency();

  public:
    ParallelConcurrencyScope() { ThreadManager::instance().set_concurrency(2u); }
    ~ParallelConcurrencyScope() { ThreadManager::instance().set_concurrency(_previous); }
};

} // namespace

class TestHybridEvolutionTree : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

  private:
    void compare_backend_orbits(HybridContinuousTimeEvolver<ValidatedTag>& evolver, const HybridEnclosure& initial,
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
        HybridContinuousTimeEvolver<ValidatedTag>& evolver, const HybridEnclosure& initial,
        const HybridTerminationCriterion& termination, Semantics semantics) {
        evolver.set_evolution_backend_factory(FlatHybridEvolutionBackendFactory());
        const auto flat = evolver.orbit(initial, termination, semantics);
        evolver.set_evolution_backend_factory(
            TreeHybridEvolutionBackendFactory(ImmediateCheckpointAdmissionPolicy{}, NoRecombinationPolicy{}));
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
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            evolver.configuration().set_maximum_step_size(0.25_x);
            const auto initial = evolver.enclosure(HybridExactBox(
                source, RealSpace({
                            x, y
            }),
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
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        evolver.configuration().set_enable_subdivisions(true);
        evolver.configuration().set_maximum_enclosure_radius(0.25_x);
        const auto wide = evolver.enclosure(HybridExactBox(
            source, RealSpace({
                        x, y
        }),
            ExactBoxType({{-1, 1}, {0, 0}})));
        for (auto semantics: {Semantics::UPPER, Semantics::LOWER}) {
            compare_flat_and_frontier_orbits(evolver, wide, HybridTime(0.03125_x, 1), semantics);
            compare_flat_and_frontier_orbits(evolver, wide, HybridTime(0, 1), semantics);
        }
        system.new_invariant(source, x <= -0.5_dec, invariant);
        HybridContinuousTimeEvolver<ValidatedTag> constrained(system);
        constrained.configuration().set_enable_subdivisions(true);
        constrained.configuration().set_maximum_enclosure_radius(0.25_x);
        const auto restricted = constrained.enclosure(HybridExactBox(
            source, RealSpace({
                        x, y
        }),
            ExactBoxType({{-1, 1}, {0, 0}})));
        for (auto semantics: {Semantics::UPPER, Semantics::LOWER})
            compare_flat_and_frontier_orbits(
                constrained, restricted, HybridTime(0.03125_x, 1), semantics);
    }

    void test_real_frontier_recombination_after_continuous_subdivision() {
        ConcurrencyScope concurrency;
        DiscreteLocation q(1);
        RealVariable x("x"), y("y");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0});
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        evolver.configuration().set_enable_subdivisions(true);
        evolver.configuration().set_maximum_enclosure_radius(0.25_x);
        evolver.configuration().set_maximum_step_size(0.25_x);
        const auto initial = evolver.enclosure(HybridExactBox(
            q, RealSpace({
                   x, y
        }),
            ExactBoxType({{-1, 1}, {0, 0}})));

        auto audit = std::make_shared<FrontierIntegrationAudit>();
        const SingleWindowFrontierPolicy grouping(audit);
        evolver.set_evolution_backend_factory(
            AuditedFrontierBackendFactory(grouping, {}, audit));
        // A zero discrete-step limit stops before continuous subdivision.
        // One step allows continuous evolution; this system has no discrete transitions.
        const auto orbit = evolver.orbit(initial, HybridTime(0.03125_x, 1), Semantics::UPPER);

        // A real numerical split must have produced a compatible candidate,
        // the production recombiner must have been invoked, and original work
        // must have been replaced without erasing processed contributions.
        ARIADNE_TESTING_PROVE(audit->split_successors >= 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(audit->selected_groups, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(audit->numerical_calls, 1u);
        ARIADNE_TESTING_PROVE(audit->replacement_trajectories >= 1u);
        ARIADNE_TESTING_PROVE(audit->retired_checkpoints >= 2u);
        ARIADNE_TESTING_PROVE(audit->processed_checkpoints > 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(audit->pending_checkpoints, 0u);
        ARIADNE_TESTING_PROVE(audit->all_trajectories_complete);
        ARIADNE_TESTING_PROVE(audit->valid_provenance);
        ARIADNE_TESTING_PROVE(orbit.final().size() > 0u);
    }

    void test_real_frontier_recombination_radius_rejection_preserves_work() {
        ConcurrencyScope concurrency;
        DiscreteLocation q(1);
        RealVariable x("x"), y("y");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0});
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        evolver.configuration().set_enable_subdivisions(true);
        evolver.configuration().set_maximum_enclosure_radius(0.25_x);
        evolver.configuration().set_maximum_step_size(0.25_x);
        const auto initial = evolver.enclosure(HybridExactBox(
            q, RealSpace({
                   x, y
        }),
            ExactBoxType({{-1, 1}, {0, 0}})));

        evolver.set_evolution_backend_factory(
            TreeHybridEvolutionBackendFactory(ImmediateCheckpointAdmissionPolicy{}, NoRecombinationPolicy{}));
        // Use the same positive discrete-step limit in both runs: zero stops
        // before continuous subdivision, even without discrete transitions.
        const auto baseline = evolver.orbit(initial, HybridTime(0.03125_x, 1), Semantics::UPPER);

        auto audit = std::make_shared<FrontierIntegrationAudit>();
        const SingleWindowFrontierPolicy grouping(audit);
        evolver.set_evolution_backend_factory(AuditedFrontierBackendFactory(
            grouping, EnclosureRecombinationContext(0.0_x), audit));
        const auto actual = evolver.orbit(initial, HybridTime(0.03125_x, 1), Semantics::UPPER);

        ARIADNE_TESTING_PROVE(audit->split_successors >= 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(audit->selected_groups, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(audit->numerical_calls, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(audit->replacement_trajectories, 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(audit->retired_checkpoints, 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(audit->pending_checkpoints, 0u);
        ARIADNE_TESTING_PROVE(audit->all_trajectories_complete);
        ARIADNE_TESTING_EXPECT_EQUAL(actual.reach().size(), baseline.reach().size());
        ARIADNE_TESTING_EXPECT_EQUAL(actual.intermediate().size(), baseline.intermediate().size());
        ARIADNE_TESTING_EXPECT_EQUAL(actual.final().size(), baseline.final().size());
        ARIADNE_TESTING_PROVE(descriptions(actual.reach()) == descriptions(baseline.reach()));
        ARIADNE_TESTING_PROVE(descriptions(actual.intermediate()) == descriptions(baseline.intermediate()));
        ARIADNE_TESTING_PROVE(descriptions(actual.final()) == descriptions(baseline.final()));
    }

    void test_real_frontier_enabled_across_urgent_and_permissive_resets() {
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
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            // The wide crossing interval is visited over several continuous
            // steps before every piece can reach the guard.
            evolver.configuration().set_maximum_step_size(0.125_x);
            const auto initial = evolver.enclosure(HybridExactBox(
                source, RealSpace({
                            x, y
            }),
                ExactBoxType({{-1, -0.75_x}, {0, 0}})));

            auto audit = std::make_shared<FrontierIntegrationAudit>();
            const SingleWindowFrontierPolicy grouping(audit);
            evolver.set_evolution_backend_factory(
                AuditedFrontierBackendFactory(grouping, {}, audit));
            const auto orbit = evolver.orbit(initial, HybridTime(2, 1), Semantics::UPPER);
            ARIADNE_TESTING_PROVE(audit->processed_checkpoints > 1u);
            ARIADNE_TESTING_PROVE(audit->transition_successors > 0u);
            ARIADNE_TESTING_PROVE(audit->all_trajectories_complete);
            ARIADNE_TESTING_PROVE(audit->valid_provenance);
            ARIADNE_TESTING_EXPECT_EQUAL(audit->pending_checkpoints, 0u);
            std::size_t target_final = 0u;
            std::size_t source_final = 0u;
            for (const auto& final: orbit.final()) {
                if (final.location() == target) {
                    ++target_final;
                    ARIADNE_TESTING_PROVE(!final.previous_events().empty());
                    ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().back(), jump);
                } else {
                    ARIADNE_TESTING_EXPECT_EQUAL(final.location(), source);
                    ++source_final;
                    ARIADNE_TESTING_PROVE(final.previous_events().empty());
                }
            }
            ARIADNE_TESTING_PROVE(target_final > 0u);
            if (event_kind == EventKind::PERMISSIVE)
                ARIADNE_TESTING_PROVE(source_final > 0u);
        }
    }

    void test_real_frontier_recombines_only_after_independent_continuous_steps() {
        ConcurrencyScope concurrency;
        DiscreteLocation q(1);
        RealVariable x("x"), y("y");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0});
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        evolver.configuration().set_enable_subdivisions(true);
        evolver.configuration().set_maximum_enclosure_radius(0.75_x);
        evolver.configuration().set_maximum_step_size(0.0625_x);
        const auto initial = evolver.enclosure(HybridExactBox(
            q, RealSpace({
                   x, y
        }),
            ExactBoxType({{0, 0}, {-1, 1}})));

        auto audit = std::make_shared<FrontierIntegrationAudit>();
        // Require the candidate checkpoint to be an actual continuation,
        // already advanced to strictly positive physical time.
        const SingleWindowFrontierPolicy grouping(audit, true);
        evolver.set_evolution_backend_factory(
            AuditedFrontierBackendFactory(grouping, {}, audit));
        const auto orbit = evolver.orbit(initial, HybridTime(0.25_x, 1), Semantics::UPPER);

        ARIADNE_TESTING_PROVE(audit->split_successors >= 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(audit->selected_groups, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(audit->numerical_calls, 1u);
        ARIADNE_TESTING_PROVE(audit->replacement_trajectories >= 1u);
        ARIADNE_TESTING_PROVE(audit->retired_checkpoints >= 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(audit->pending_checkpoints, 0u);
        ARIADNE_TESTING_PROVE(audit->all_trajectories_complete);
        ARIADNE_TESTING_PROVE(audit->valid_provenance);
        ARIADNE_TESTING_PROVE(orbit.final().size() > 0u);
    }

    void test_real_post_reset_recombination_of_split_trajectories() {
        ConcurrencyScope concurrency;
        // Splitting along y leaves x and its guard-crossing time identical in
        // both branches. The numerical evolver performs the actual crossing;
        // the test policy only considers post-reset checkpoints.
        for (auto event_kind: {EventKind::PERMISSIVE, EventKind::URGENT}) {
            DiscreteLocation source(1), target(2);
            RealVariable x("x"), y("y");
            DiscreteEvent jump("jump");
            HybridAutomaton system;
            system.new_mode(source, {dot(x) = 1, dot(y) = 0});
            system.new_mode(target, {dot(x) = 0, dot(y) = 0});
            system.new_transition(source, jump, target,
                                  {next(x) = x, next(y) = y}, x >= 0, event_kind);
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            evolver.configuration().set_enable_subdivisions(true);
            evolver.configuration().set_maximum_enclosure_radius(0.5_x);
            evolver.configuration().set_maximum_step_size(0.125_x);
            const auto initial = evolver.enclosure(HybridExactBox(
                source, RealSpace({
                            x, y
            }),
                ExactBoxType({{-0.5_x, -0.5_x}, {-1, 1}})));

            auto audit = std::make_shared<FrontierIntegrationAudit>();
            const SingleWindowFrontierPolicy grouping(audit, false, true);
            evolver.set_evolution_backend_factory(
                AuditedFrontierBackendFactory(grouping, {}, audit));
            const auto orbit = evolver.orbit(initial, HybridTime(0.75_x, 1), Semantics::UPPER);

            ARIADNE_TESTING_PROVE(audit->split_successors >= 2u);
            ARIADNE_TESTING_PROVE(audit->transition_successors >= 2u);
            ARIADNE_TESTING_PROVE(audit->selection_checks > 0u);
            ARIADNE_TESTING_EXPECT_EQUAL(audit->selected_groups, 1u);
            ARIADNE_TESTING_EXPECT_EQUAL(audit->numerical_calls, 1u);
            ARIADNE_TESTING_PROVE(audit->replacement_trajectories >= 1u);
            ARIADNE_TESTING_PROVE(audit->retired_checkpoints >= 2u);
            ARIADNE_TESTING_PROVE(audit->processed_checkpoints > 2u);
            ARIADNE_TESTING_EXPECT_EQUAL(audit->pending_checkpoints, 0u);
            ARIADNE_TESTING_PROVE(audit->all_trajectories_complete);
            ARIADNE_TESTING_PROVE(audit->valid_provenance);
            std::size_t target_finals = 0u;
            for (const auto& final: orbit.final()) {
                if (final.location() == target) {
                    ++target_finals;
                    ARIADNE_TESTING_PROVE(!final.previous_events().empty());
                    ARIADNE_TESTING_EXPECT_EQUAL(final.previous_events().back(), jump);
                }
            }
            ARIADNE_TESTING_PROVE(target_finals > 0u);
        }
    }

    void test_frontier_recombination_requires_progress_after_resplitting() {
        const auto initial = seed();
        TreeHybridEvolutionBackend backend(
            initial, ImmediateCheckpointAdmissionPolicy{}, WindowRecombinationPolicy(2u, 0.0),
            AfterEachStepRecombinationTriggerPolicy{});
        backend.set_semantics(Semantics::UPPER);
        std::size_t calls = 0u;
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            if (calls++ < 2u) {
                result.continue_with(item.enclosure, HybridEvolutionSuccessorKind::SUBDIVISION);
                result.continue_with(item.enclosure, HybridEvolutionSuccessorKind::SUBDIVISION);
            } else {
                result.record_final(item.enclosure);
            }
            return result;
        });
        const auto top = backend.tree().root().value();
        std::size_t replacements = 0u, retired = 0u;
        const auto nodes = backend.tree().descendants(top);
        for (const auto node: nodes) {
            if (!backend.tree().metadata(node).recombined_checkpoints.empty())
                ++replacements;
            retired += backend.tree().object(node).retired_entries.size();
            ARIADNE_TESTING_PROVE(backend.tree().metadata(node).complete);
        }
        ARIADNE_TESTING_EXPECT_EQUAL(replacements, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(retired, 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(calls, 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 2u);
    }

    void test_bounded_frontier_wait_recovers_late_merge_partners() {
        const auto initial = seed();
        const auto advance = [&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            const int value = label(item);
            HybridEvolutionStepResult result;
            if (value == 0) {
                // The more advanced branch enters the queue first.
                result.continue_with(marked(initial, 2), HybridEvolutionSuccessorKind::SUBDIVISION);
                result.continue_with(marked(initial, 1), HybridEvolutionSuccessorKind::SUBDIVISION);
            } else if (value == 1) {
                result.continue_with(marked(initial, 2));
            } else if (value == 2) {
                result.record_final(item.enclosure);
            } else {
                throw std::logic_error("Unexpected checkpoint");
            }
            return result;
        };
        TreeHybridEvolutionBackend immediate(initial, ImmediateCheckpointAdmissionPolicy{}, WindowRecombinationPolicy(2u, 0.0));
        immediate.set_semantics(Semantics::UPPER);
        std::vector<int> direct_order;
        immediate.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics& stats) {
            direct_order.push_back(label(item));
            return advance(item, stats);
        });
        ARIADNE_TESTING_PROVE(direct_order == std::vector<int>({0, 2, 1, 2}));
        ARIADNE_TESTING_EXPECT_EQUAL(immediate.release_orbit().final().size(), 2u);

        TreeHybridEvolutionBackend waiting(
            initial, ImmediateCheckpointAdmissionPolicy{}, WindowRecombinationPolicy(2u, 0.0),
            AfterEachStepRecombinationTriggerPolicy{}, ReferenceEnclosureRecombiner{}, {},
            BoundedTimeLagCheckpointDispatchPolicy(2u, 2.0));
        waiting.set_semantics(Semantics::UPPER);
        std::vector<int> delayed_order;
        waiting.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics& stats) {
            delayed_order.push_back(label(item));
            return advance(item, stats);
        });
        ARIADNE_TESTING_PROVE(delayed_order == std::vector<int>({0, 1, 2}));
        const auto root = waiting.tree().root().value();
        std::size_t replacements = 0u, retired = 0u;
        for (const auto node: waiting.tree().descendants(root)) {
            if (!waiting.tree().metadata(node).recombined_checkpoints.empty())
                ++replacements;
            retired += waiting.tree().object(node).retired_entries.size();
        }
        ARIADNE_TESTING_EXPECT_EQUAL(replacements, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(retired, 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(waiting.release_orbit().final().size(), 1u);
    }

    void test_bounded_frontier_wait_does_not_starve_unmatched_work() {
        const auto initial = seed();
        TreeHybridEvolutionBackend backend(
            initial, ImmediateCheckpointAdmissionPolicy{}, NoRecombinationPolicy{},
            AfterEachStepRecombinationTriggerPolicy{}, ReferenceEnclosureRecombiner{}, {},
            BoundedTimeLagCheckpointDispatchPolicy(1u, 10.0));
        std::vector<int> order;
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            const int value = label(item);
            order.push_back(value);
            HybridEvolutionStepResult result;
            if (value == 0) {
                result.continue_with(marked(initial, 3), HybridEvolutionSuccessorKind::SUBDIVISION);
                result.continue_with(marked(initial, 1), HybridEvolutionSuccessorKind::SUBDIVISION);
            } else {
                result.record_final(item.enclosure);
            }
            return result;
        });
        ARIADNE_TESTING_PROVE(order == std::vector<int>({0, 1, 3}));
        ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 2u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(BoundedTimeLagCheckpointDispatchPolicy(1u, -1.0));
    }

    void test_radius_trend_gate_accepts_only_profitable_continuations() {
        DiscreteLocation q(1);
        RealVariable x("x"), y("y");
        HybridAutomaton system;
        system.new_mode(q, {dot(x) = 1, dot(y) = 0});
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
        const auto initial = evolver.enclosure(HybridExactBox(
            q, RealSpace({
                   x, y
        }),
            ExactBoxType({{-1, 1}, {0, 0}})));
        const auto advance = [&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            if (label(item) == 0) {
                result.continue_with(marked(initial, 1), HybridEvolutionSuccessorKind::SUBDIVISION);
                result.continue_with(marked(initial, 1), HybridEvolutionSuccessorKind::SUBDIVISION);
            } else if (label(item) == 1) {
                result.continue_with(marked(initial, 2));
            } else {
                result.record_final(item.enclosure);
            }
            return result;
        };
        const auto verify = [&](const HybridRecombinationAcceptancePolicy& acceptance,
                                std::size_t expected_replacements, std::size_t expected_finals) {
            TreeHybridEvolutionBackend backend(
                initial, ImmediateCheckpointAdmissionPolicy{}, WindowRecombinationPolicy(2u, 0.0),
                EveryNStepsRecombinationTriggerPolicy(3u), ReferenceEnclosureRecombiner{}, {},
                ImmediateCheckpointDispatchPolicy{}, acceptance);
            backend.set_semantics(Semantics::UPPER);
            backend.run(advance);
            const auto root = backend.tree().root().value();
            std::size_t replacements = 0u;
            for (const auto node: backend.tree().descendants(root))
                if (!backend.tree().metadata(node).recombined_checkpoints.empty())
                    ++replacements;
            ARIADNE_TESTING_EXPECT_EQUAL(replacements, expected_replacements);
            ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), expected_finals);
        };
        verify(RadiusTrendRecombinationAcceptancePolicy(100.0, 1.0), 1u, 1u);
        verify(RadiusTrendRecombinationAcceptancePolicy(0.1, 1.0), 0u, 2u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(RadiusTrendRecombinationAcceptancePolicy(0.0, 1.0));
        ARIADNE_TESTING_EXPECT_EXCEPTION(RadiusTrendRecombinationAcceptancePolicy(2.0, -1.0));
    }

    void test_parallel_parent_first_last_publisher_merges_held_children() {
        if (ThreadManager::instance().maximum_concurrency() < 2u)
            return;
        ParallelConcurrencyScope concurrency;
        const auto initial = seed();
        const DiscreteEvent jump("jump");
        std::mutex gate_mutex;
        std::condition_variable gate;
        unsigned running = 0u;
        std::atomic<unsigned> numerical_finished{0u};
        std::atomic<unsigned> successors_executed{0u};
        std::atomic<bool> premature_release{false};
        TreeHybridEvolutionBackend backend(
            initial, ParentCheckpointAdmissionPolicy{}, AllCompatibleRecombinationPolicy{},
            BeforeAdmissionRecombinationTriggerPolicy{});
        backend.set_semantics(Semantics::UPPER);
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            switch (label(item)) {
                case 0:
                    result.continue_with(marked(initial, 1), HybridEvolutionSuccessorKind::SUBDIVISION);
                    result.continue_with(marked(initial, 2), HybridEvolutionSuccessorKind::SUBDIVISION);
                    break;
                case 1:
                case 2: {
                    // Both numerical tasks must be simultaneously active:
                    // a serial implementation fails the bounded wait.
                    std::unique_lock<std::mutex> lock(gate_mutex);
                    ++running;
                    gate.notify_all();
                    if (!gate.wait_for(lock, std::chrono::seconds(15),
                                       [&] { return running == 2u; }))
                        throw std::runtime_error("Parent split pieces were not executed in parallel");
                    lock.unlock();
                    result.transition_to(marked(initial, 100), jump, EventKind::URGENT);
                    ++numerical_finished;
                    break;
                }
                case 100:
                    if (numerical_finished.load() != 2u)
                        premature_release = true;
                    ++successors_executed;
                    result.record_final(item.enclosure);
                    break;
                default:
                    throw std::logic_error("Unexpected parallel parent-first checkpoint");
            }
            return result;
        });
        const auto root = backend.tree().root().value();
        const auto children = backend.tree().children(root);
        ARIADNE_TESTING_EXPECT_EQUAL(running, 2u);
        ARIADNE_TESTING_PROVE(!premature_release.load());
        ARIADNE_TESTING_EXPECT_EQUAL(successors_executed.load(), 1u);
        ARIADNE_TESTING_PROVE(backend.tree().metadata(root).complete);
        ARIADNE_TESTING_EXPECT_EQUAL(children.size(), 3u);
        ARIADNE_TESTING_PROVE(!backend.tree().object(backend.tree().children(children[0]).front()).retired_entries.empty());
        ARIADNE_TESTING_PROVE(!backend.tree().object(backend.tree().children(children[1]).front()).retired_entries.empty());
        ARIADNE_TESTING_EXPECT_EQUAL(
            backend.tree().metadata(children[2]).recombined_checkpoints.size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 1u);
    }

    void test_parallel_frontier_executes_independent_split_checkpoints() {
        if (ThreadManager::instance().maximum_concurrency() < 2u)
            return;
        ParallelConcurrencyScope concurrency;
        const auto initial = seed();
        std::mutex gate_mutex;
        std::condition_variable gate;
        unsigned running = 0u;
        TreeHybridEvolutionBackend backend(initial, ImmediateCheckpointAdmissionPolicy{}, NoRecombinationPolicy{});
        backend.set_semantics(Semantics::UPPER);
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            switch (label(item)) {
                case 0:
                    result.continue_with(marked(initial, 1), HybridEvolutionSuccessorKind::SUBDIVISION);
                    result.continue_with(marked(initial, 2), HybridEvolutionSuccessorKind::SUBDIVISION);
                    break;
                case 1:
                case 2: {
                    std::unique_lock<std::mutex> lock(gate_mutex);
                    ++running;
                    gate.notify_all();
                    if (!gate.wait_for(lock, std::chrono::seconds(15),
                                       [&] { return running == 2u; }))
                        throw std::runtime_error("Frontier split pieces were not executed in parallel");
                    lock.unlock();
                    result.record_final(item.enclosure);
                    break;
                }
                default:
                    throw std::logic_error("Unexpected parallel frontier checkpoint");
            }
            return result;
        });
        const auto root = backend.tree().root().value();
        ARIADNE_TESTING_EXPECT_EQUAL(running, 2u);
        for (const auto node: backend.tree().descendants(root))
            ARIADNE_TESTING_PROVE(backend.tree().metadata(node).complete);
        ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 2u);
    }

    void test_location_batch_merges_different_parents_after_all_producers_finish() {
        ConcurrencyScope concurrency;
        const auto initial = seed();
        const DiscreteLocation target(2);
        const DiscreteEvent jump("jump");
        TreeHybridEvolutionBackend backend(initial, LocationCheckpointAdmissionPolicy{}, AllCompatibleRecombinationPolicy{});
        backend.set_semantics(Semantics::UPPER);
        unsigned producers_finished = 0u, target_executions = 0u;
        std::vector<int> source_order;
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            if (item.enclosure.location() == target) {
                ARIADNE_TESTING_EXPECT_EQUAL(producers_finished, 2u);
                ++target_executions;
                result.record_final(item.enclosure);
            } else {
                source_order.push_back(label(item));
                switch (label(item)) {
                    case 0:
                        result.continue_with(marked(initial, 1), HybridEvolutionSuccessorKind::SUBDIVISION);
                        result.continue_with(marked(initial, 2), HybridEvolutionSuccessorKind::SUBDIVISION);
                        break;
                    case 1:
                        result.transition_to(jumped(item.enclosure, target, jump, 10), jump, EventKind::URGENT);
                        ++producers_finished;
                        break;
                    case 2:
                        result.continue_with(marked(item.enclosure, 3));
                        break;
                    case 3:
                        result.transition_to(jumped(item.enclosure, target, jump, 20), jump, EventKind::URGENT);
                        ++producers_finished;
                        break;
                    default:
                        throw std::logic_error("Unexpected location-batch source");
                }
            }
            return result;
        });
        ARIADNE_TESTING_PROVE(source_order == std::vector<int>({0, 1, 2, 3}));
        ARIADNE_TESTING_EXPECT_EQUAL(target_executions, 1u);
        std::size_t replacements = 0u;
        for (const auto node: backend.tree().descendants(backend.tree().root().value())) {
            const auto& sources = backend.tree().metadata(node).recombined_checkpoints;
            if (sources.empty())
                continue;
            ++replacements;
            ARIADNE_TESTING_EXPECT_EQUAL(sources.size(), 2u);
            ARIADNE_TESTING_PROVE(backend.tree().parent(sources[0].trajectory) != backend.tree().parent(sources[1].trajectory));
        }
        ARIADNE_TESTING_EXPECT_EQUAL(replacements, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 1u);
    }

    void test_location_early_replacements_stay_held_and_absorb_later_arrivals() {
        ConcurrencyScope concurrency;
        const auto initial = seed();
        const DiscreteLocation target(2);
        const DiscreteEvent jump("jump");
        TreeHybridEvolutionBackend backend(initial, LocationCheckpointAdmissionPolicy{},
                                           AllCompatibleRecombinationPolicy{}, AfterEachStepRecombinationTriggerPolicy{});
        backend.set_semantics(Semantics::UPPER);
        unsigned finished = 0u, executions = 0u;
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            if (item.enclosure.location() == target) {
                ARIADNE_TESTING_EXPECT_EQUAL(finished, 3u);
                ++executions;
                result.record_final(item.enclosure);
            } else if (label(item) == 0) {
                for (int index = 1; index <= 3; ++index)
                    result.continue_with(marked(initial, index), HybridEvolutionSuccessorKind::SUBDIVISION);
            } else if (label(item) == 3) {
                result.continue_with(marked(item.enclosure, 4));
            } else {
                result.transition_to(jumped(item.enclosure, target, jump, label(item) * 10), jump, EventKind::URGENT);
                ++finished;
            }
            return result;
        });
        ARIADNE_TESTING_EXPECT_EQUAL(executions, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 1u);
    }

    void test_location_batches_keep_destinations_separate_and_reenter_later() {
        ConcurrencyScope concurrency;
        const auto initial = seed();
        const DiscreteLocation second(2), third(3);
        const DiscreteEvent to_second("to-second"), to_third("to-third"), again("again");
        TreeHybridEvolutionBackend backend(initial, LocationCheckpointAdmissionPolicy{}, AllCompatibleRecombinationPolicy{});
        backend.set_semantics(Semantics::UPPER);
        unsigned source_finished = 0u;
        std::vector<DiscreteLocation> batches;
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            if (!item.enclosure.previous_events().empty()) {
                ARIADNE_TESTING_EXPECT_EQUAL(source_finished, 2u);
                batches.push_back(item.enclosure.location());
                result.record_final(item.enclosure);
            } else if (label(item) == 0) {
                result.continue_with(marked(initial, 1), HybridEvolutionSuccessorKind::SUBDIVISION);
                result.continue_with(marked(initial, 2), HybridEvolutionSuccessorKind::SUBDIVISION);
            } else {
                const auto time_label = label(item) * 10;
                result.transition_to(jumped(item.enclosure, second, to_second, time_label), to_second, EventKind::PERMISSIVE);
                result.transition_to(jumped(item.enclosure, third, to_third, time_label), to_third, EventKind::PERMISSIVE);
                result.transition_to(jumped(item.enclosure, initial.location(), again, time_label), again, EventKind::PERMISSIVE);
                ++source_finished;
            }
            return result;
        });
        ARIADNE_TESTING_PROVE(batches == std::vector<DiscreteLocation>({second, third, initial.location()}));
        ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 3u);
    }

    void test_parallel_location_batch_waits_for_in_flight_producers() {
        if (ThreadManager::instance().maximum_concurrency() < 2u)
            return;
        ParallelConcurrencyScope concurrency;
        const auto initial = seed();
        const DiscreteLocation target(2);
        const DiscreteEvent jump("jump");
        std::mutex gate_mutex;
        std::condition_variable gate;
        unsigned running = 0u;
        std::atomic<unsigned> finished{0u}, dispatched{0u};
        std::atomic<bool> premature{false};
        TreeHybridEvolutionBackend backend(initial, LocationCheckpointAdmissionPolicy{}, AllCompatibleRecombinationPolicy{});
        backend.set_semantics(Semantics::UPPER);
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            if (item.enclosure.location() == target) {
                if (finished.load() != 2u)
                    premature = true;
                ++dispatched;
                result.record_final(item.enclosure);
            } else if (label(item) == 0) {
                result.continue_with(marked(initial, 1), HybridEvolutionSuccessorKind::SUBDIVISION);
                result.continue_with(marked(initial, 2), HybridEvolutionSuccessorKind::SUBDIVISION);
            } else {
                std::unique_lock<std::mutex> lock(gate_mutex);
                ++running;
                gate.notify_all();
                if (!gate.wait_for(lock, std::chrono::seconds(15), [&] { return running == 2u; }))
                    throw std::runtime_error("Location pieces did not execute concurrently");
                lock.unlock();
                result.transition_to(jumped(item.enclosure, target, jump, label(item) * 10), jump, EventKind::URGENT);
                ++finished;
            }
            return result;
        });
        ARIADNE_TESTING_PROVE(!premature.load());
        ARIADNE_TESTING_EXPECT_EQUAL(dispatched.load(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 1u);
    }

    void test_location_factory_preserves_real_orbits_without_recombination() {
        ConcurrencyScope concurrency;
        const DiscreteLocation source(1), target(2);
        const DiscreteEvent jump("jump");
        RealVariable x("x"), y("y");
        for (const auto kind: {EventKind::URGENT, EventKind::PERMISSIVE}) {
            HybridAutomaton system;
            system.new_mode(source, {dot(x) = 1, dot(y) = sqr(x)});
            system.new_mode(target, {dot(x) = 1 + sqr(y), dot(y) = 1});
            system.new_transition(source, jump, target, {next(x) = x, next(y) = y}, x >= 0, kind);
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
            evolver.configuration().set_maximum_step_size(0.03125_x);
            evolver.configuration().set_maximum_enclosure_radius(0.03125_x);
            const auto initial = evolver.enclosure(HybridExactBox(source, RealSpace({
                                                                              x, y
            }),
                                                                  ExactBoxType({{-0.0625_x, 0.0625_x}, {0, 0}})));
            const auto expected = evolver.orbit(initial, HybridTime(0.125_x, 2), Semantics::UPPER);
            evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory(
                LocationCheckpointAdmissionPolicy{}, NoRecombinationPolicy{}));
            const auto actual = evolver.orbit(initial, HybridTime(0.125_x, 2), Semantics::UPPER);
            ARIADNE_TESTING_PROVE(descriptions(expected.reach()) == descriptions(actual.reach()));
            ARIADNE_TESTING_PROVE(descriptions(expected.intermediate()) == descriptions(actual.intermediate()));
            ARIADNE_TESTING_PROVE(descriptions(expected.final()) == descriptions(actual.final()));
            const auto audit = std::make_shared<FrontierIntegrationAudit>();
            evolver.set_evolution_backend_factory(AuditedFrontierBackendFactory(
                LocationCheckpointAdmissionPolicy{}, AllCompatibleRecombinationPolicy{}, {}, audit));
            const auto recombined = evolver.orbit(initial, HybridTime(0.125_x, 2), Semantics::UPPER);
            ARIADNE_TESTING_PROVE(audit->numerical_calls > 0u);
            ARIADNE_TESTING_PROVE(audit->replacement_trajectories > 0u);
            ARIADNE_TESTING_PROVE(audit->retired_checkpoints >= 2u);
            ARIADNE_TESTING_EXPECT_EQUAL(audit->pending_checkpoints, 0u);
            ARIADNE_TESTING_PROVE(audit->all_trajectories_complete);
            ARIADNE_TESTING_PROVE(audit->valid_provenance);
            ARIADNE_TESTING_PROVE(recombined.final().size() > 0u);
            const auto rejected_audit = std::make_shared<FrontierIntegrationAudit>();
            evolver.set_evolution_backend_factory(AuditedFrontierBackendFactory(
                LocationCheckpointAdmissionPolicy{}, AllCompatibleRecombinationPolicy{},
                EnclosureRecombinationContext(0.0_x), rejected_audit));
            const auto rejected = evolver.orbit(initial, HybridTime(0.125_x, 2), Semantics::UPPER);
            ARIADNE_TESTING_PROVE(rejected_audit->numerical_calls > 0u);
            ARIADNE_TESTING_EXPECT_EQUAL(rejected_audit->replacement_trajectories, 0u);
            ARIADNE_TESTING_EXPECT_EQUAL(rejected_audit->retired_checkpoints, 0u);
            ARIADNE_TESTING_PROVE(descriptions(expected.reach()) == descriptions(rejected.reach()));
            ARIADNE_TESTING_PROVE(descriptions(expected.intermediate()) == descriptions(rejected.intermediate()));
            ARIADNE_TESTING_PROVE(descriptions(expected.final()) == descriptions(rejected.final()));
        }
    }

    void test_progress_recombines_waiting_time_without_inventing_evolution() {
        ConcurrencyScope concurrency;
        const auto initial = seed();
        const DiscreteLocation target(2);
        const DiscreteEvent jump("jump");
        TreeHybridEvolutionBackend backend(initial, LocationCheckpointAdmissionPolicy{}, AllCompatibleRecombinationPolicy{});
        backend.set_time_horizon(Real(10));
        backend.set_semantics(Semantics::UPPER);
        bool split = false;
        unsigned targets = 0u;
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            const auto progress = backend.progress();
            if (item.enclosure.location() == target) {
                ARIADNE_TESTING_EXPECT_EQUAL(progress.processed_time, 4.0);
                ARIADNE_TESTING_EXPECT_EQUAL(progress.pending, 1u);
                ARIADNE_TESTING_EXPECT_EQUAL(progress.remaining_time.value(), 7.0);
                if (item.phase == HybridEvolutionPhase::PROCESS_STARTING_EVENTS) {
                    ++targets;
                    result.continue_with(item.enclosure);
                } else {
                    result.continuous_evolution_time = 7.0;
                    result.record_final(marked(item.enclosure, 10));
                }
            } else if (item.phase == HybridEvolutionPhase::PROCESS_STARTING_EVENTS) {
                result.continue_with(initial);
            } else if (label(item) == 0) {
                result.continuous_evolution_time = 2.0;
                result.continue_with(marked(initial, 2));
            } else if (!split) {
                split = true;
                result.continue_with(item.enclosure, HybridEvolutionSuccessorKind::SUBDIVISION);
                result.continue_with(item.enclosure, HybridEvolutionSuccessorKind::SUBDIVISION);
            } else {
                // One pending set is still in the source; the earlier reset
                // image already waits in the target and must remain counted.
                ARIADNE_TESTING_EXPECT_EQUAL(progress.pending, 2u);
                ARIADNE_TESTING_EXPECT_EQUAL(progress.remaining_time.value(), progress.processed_time == 2.0 ? 16.0 : 15.0);
                result.continuous_evolution_time = 1.0;
                result.transition_to(jumped(item.enclosure, target, jump, 3), jump, EventKind::URGENT);
            }
            return result;
        });
        const auto progress = backend.progress();
        ARIADNE_TESTING_EXPECT_EQUAL(targets, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(progress.processed_time, 11.0);
        ARIADNE_TESTING_EXPECT_EQUAL(progress.remaining_time.value(), 0.0);
        ARIADNE_TESTING_PROVE(progress.complete);
        ARIADNE_TESTING_EXPECT_EXCEPTION(backend.set_time_horizon(Real(20)));
    }

    void test_progress_counts_parallel_in_flight_work() {
        if (ThreadManager::instance().maximum_concurrency() < 2u)
            return;
        ParallelConcurrencyScope concurrency;
        const auto initial = seed();
        std::vector<std::shared_ptr<const HybridEvolutionBackendFactoryInterface>> factories{
            std::make_shared<FlatHybridEvolutionBackendFactory>(),
            std::make_shared<TreeHybridEvolutionBackendFactory>(),
            std::make_shared<TreeHybridEvolutionBackendFactory>(ImmediateCheckpointAdmissionPolicy{}),
            std::make_shared<TreeHybridEvolutionBackendFactory>(LocationCheckpointAdmissionPolicy{})};
        for (const auto& factory: factories) {
            auto backend = factory->create(initial, Real(10));
            std::mutex gate_mutex;
            std::condition_variable gate;
            unsigned entered = 0u;
            HybridEvolutionProgress simultaneous;
            backend->run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
                HybridEvolutionStepResult result;
                if (label(item) == 0) {
                    result.continuous_evolution_time = 2.0;
                    result.continue_with(marked(initial, 2), HybridEvolutionSuccessorKind::SUBDIVISION);
                    result.continue_with(marked(initial, 2), HybridEvolutionSuccessorKind::SUBDIVISION);
                } else {
                    std::unique_lock<std::mutex> lock(gate_mutex);
                    ++entered;
                    if (entered == 2u)
                        simultaneous = backend->progress();
                    gate.notify_all();
                    if (!gate.wait_for(lock, std::chrono::seconds(15), [&] { return entered == 2u; }))
                        throw std::runtime_error("Parallel progress regression did not claim both pieces");
                    result.continuous_evolution_time = 8.0;
                    result.record_final(marked(item.enclosure, 10));
                }
                return result;
            });
            ARIADNE_TESTING_EXPECT_EQUAL(simultaneous.pending, 2u);
            ARIADNE_TESTING_EXPECT_EQUAL(simultaneous.processed_time, 2.0);
            ARIADNE_TESTING_EXPECT_EQUAL(simultaneous.remaining_time.value(), 16.0);
            ARIADNE_TESTING_EXPECT_EQUAL(backend->progress().processed_time, 18.0);
            ARIADNE_TESTING_PROVE(backend->progress().complete == true);
        }
    }

    void test_parent_admission_releases_completed_families_without_global_barrier() {
        ConcurrencyScope concurrency;
        const auto initial = seed();
        TreeHybridEvolutionBackend backend(initial);
        bool slow_finished = false;
        bool fast_executed = false;
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            switch (label(item)) {
                case 0:
                    result.transition_to(marked(initial, 10), DiscreteEvent("fast"), EventKind::URGENT);
                    result.transition_to(marked(initial, 20), DiscreteEvent("slow"), EventKind::URGENT);
                    break;
                case 10:
                    result.transition_to(marked(initial, 100), DiscreteEvent("jump"), EventKind::URGENT);
                    break;
                case 20:
                    result.continue_with(marked(initial, 21));
                    break;
                case 21:
                    slow_finished = true;
                    break;
                case 100:
                    ARIADNE_TESTING_PROVE(!slow_finished);
                    fast_executed = true;
                    break;
                default:
                    throw std::logic_error("Unexpected parent-admission checkpoint");
            }
            return result;
        });
        ARIADNE_TESTING_PROVE(fast_executed && slow_finished);
    }

    void test_parent_admission_preserves_all_dependencies_of_cross_parent_replacements() {
        ConcurrencyScope concurrency;
        const auto initial = seed();
        TreeHybridEvolutionBackend backend(initial, ParentCheckpointAdmissionPolicy{},
                                           WindowRecombinationPolicy(8u, 0.0),
                                           AfterEachStepRecombinationTriggerPolicy{});
        backend.set_semantics(Semantics::UPPER);
        bool first_finished = false, second_finished = false;
        std::size_t targets = 0u;
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            switch (label(item)) {
                case 0:
                    result.transition_to(marked(initial, 10), DiscreteEvent("first"), EventKind::URGENT);
                    result.transition_to(marked(initial, 20), DiscreteEvent("second"), EventKind::URGENT);
                    break;
                case 10:
                case 20:
                    result.transition_to(marked(initial, 100), DiscreteEvent("jump"), EventKind::URGENT);
                    result.continue_with(marked(initial, label(item) + 1));
                    break;
                case 11:
                    first_finished = true;
                    break;
                case 21:
                    result.continue_with(marked(initial, 22));
                    break;
                case 22:
                    second_finished = true;
                    break;
                case 100:
                    ARIADNE_TESTING_PROVE(first_finished && second_finished);
                    ++targets;
                    result.record_final(item.enclosure);
                    break;
                default:
                    throw std::logic_error("Unexpected parent-admission checkpoint");
            }
            return result;
        });
        ARIADNE_TESTING_EXPECT_EQUAL(targets, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 1u);
    }

    void test_periodic_ready_recombination_does_not_flush_on_each_admission() {
        ConcurrencyScope concurrency;
        const auto initial = seed();
        const auto execute = [&](const HybridRecombinationTriggerPolicy& trigger) {
            TreeHybridEvolutionBackend backend(initial, ImmediateCheckpointAdmissionPolicy{},
                                               WindowRecombinationPolicy(2u, 0.0), trigger);
            backend.set_semantics(Semantics::UPPER);
            backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
                HybridEvolutionStepResult result;
                if (label(item) == 0) {
                    result.continue_with(marked(initial, 1), HybridEvolutionSuccessorKind::SUBDIVISION);
                    result.continue_with(marked(initial, 1), HybridEvolutionSuccessorKind::SUBDIVISION);
                } else {
                    result.record_final(item.enclosure);
                }
                return result;
            });
            return backend.release_orbit().final().size();
        };
        ARIADNE_TESTING_EXPECT_EQUAL(execute(BeforeAdmissionRecombinationTriggerPolicy{}), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(execute(EveryNStepsRecombinationTriggerPolicy(3u)), 2u);
    }

    void test_sibling_window_ignores_interleaved_unrelated_parents() {
        ConcurrencyScope concurrency;
        const auto initial = seed();
        TreeHybridEvolutionBackend backend(initial, ParentCheckpointAdmissionPolicy{},
                                           WindowRecombinationPolicy(2u, 0.0, HybridRecombinationScope::SIBLINGS));
        backend.set_semantics(Semantics::UPPER);
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            switch (label(item)) {
                case 0:
                    result.transition_to(marked(initial, 10), DiscreteEvent("first"), EventKind::URGENT);
                    result.transition_to(marked(initial, 20), DiscreteEvent("second"), EventKind::URGENT);
                    break;
                case 10:
                case 20:
                    result.transition_to(marked(initial, 100), DiscreteEvent("jump"), EventKind::URGENT);
                    result.continue_with(marked(initial, label(item) + 1));
                    break;
                case 11:
                case 21:
                    result.transition_to(marked(initial, 100), DiscreteEvent("jump"), EventKind::URGENT);
                    break;
                case 100:
                    result.record_final(item.enclosure);
                    break;
                default:
                    throw std::logic_error("Unexpected sibling-window checkpoint");
            }
            return result;
        });
        std::size_t merged = 0u;
        const auto root = backend.tree().root().value();
        for (const auto node: backend.tree().descendants(root)) {
            const auto& sources = backend.tree().metadata(node).recombined_checkpoints;
            if (!sources.empty()) {
                ++merged;
                ARIADNE_TESTING_EXPECT_EQUAL(sources.size(), 2u);
                ARIADNE_TESTING_PROVE(backend.tree().parent(sources[0].trajectory) == backend.tree().parent(sources[1].trajectory));
            }
        }
        ARIADNE_TESTING_EXPECT_EQUAL(merged, 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 2u);
    }

    void test_waiting_admission_rejects_continuous_only_quality_gate() {
        const auto initial = seed();
        ARIADNE_TESTING_EXPECT_EXCEPTION(TreeHybridEvolutionBackend(
            initial, LocationCheckpointAdmissionPolicy{}, AllCompatibleRecombinationPolicy{},
            BeforeAdmissionRecombinationTriggerPolicy{}, ReferenceEnclosureRecombiner{}, {},
            ImmediateCheckpointDispatchPolicy{}, RadiusTrendRecombinationAcceptancePolicy{}));
        ARIADNE_TESTING_EXPECT_EXCEPTION(TreeHybridEvolutionBackendFactory(
            ParentCheckpointAdmissionPolicy{}, AllCompatibleRecombinationPolicy{},
            BeforeAdmissionRecombinationTriggerPolicy{}, ReferenceEnclosureRecombiner{}, {},
            ImmediateCheckpointDispatchPolicy{}, RadiusTrendRecombinationAcceptancePolicy{}));
        // Disabled recombination never invokes the unused acceptance gate.
        TreeHybridEvolutionBackend backend(initial, LocationCheckpointAdmissionPolicy{}, NoRecombinationPolicy{},
                                           BeforeAdmissionRecombinationTriggerPolicy{}, ReferenceEnclosureRecombiner{}, {},
                                           ImmediateCheckpointDispatchPolicy{}, RadiusTrendRecombinationAcceptancePolicy{});
        backend.run([](const HybridEvolutionWorkItem&, const HybridEvolutionStatistics&) { return HybridEvolutionStepResult{}; });
        ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 0u);
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

        NoRecombinationPolicy noop;
        TreeHybridEvolutionBackend baseline(initial, ParentCheckpointAdmissionPolicy{}, noop);
        baseline.set_semantics(Semantics::UPPER);
        baseline.run(advance);
        ARIADNE_TESTING_EXPECT_EQUAL(baseline.tree().size(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(baseline.release_orbit().final().size(), 3u);

        WindowRecombinationPolicy window(2u, 0.0, HybridRecombinationScope::SIBLINGS);
        TreeHybridEvolutionBackend recombined(initial, ParentCheckpointAdmissionPolicy{}, window);
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
        ARIADNE_TESTING_PROVE(!recombined.tree().object(children[0]).retired_entries.empty());
        ARIADNE_TESTING_PROVE(!recombined.tree().object(children[1]).retired_entries.empty());
        ARIADNE_TESTING_PROVE(!recombined.tree().object(children[0]).entries[0].result);
        ARIADNE_TESTING_PROVE(!recombined.tree().object(children[1]).entries[0].result);
        ARIADNE_TESTING_PROVE(recombined.tree().object(children[2]).retired_entries.empty());
        const auto& replacement = recombined.tree().metadata(children[3]);
        ARIADNE_TESTING_PROVE(replacement.complete);
        ARIADNE_TESTING_PROVE(!replacement.origin);
        ARIADNE_TESTING_EXPECT_EQUAL(replacement.recombined_checkpoints.size(), 2u);
        ARIADNE_TESTING_PROVE(replacement.recombined_checkpoints[0].trajectory == children[0]);
        ARIADNE_TESTING_PROVE(replacement.recombined_checkpoints[1].trajectory == children[1]);
        const auto orbit = recombined.release_orbit();
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.reach().size(), 0u);
    }

    void test_proximity_bounds_normalized_diameter_radius_and_cardinality() {
        const auto initial = seed();
        const auto boxed = [&](ExactBoxType box) {
            return HybridEnclosure(HybridExactBox(initial.location(), initial.state_space(), box), initial.configuration());
        };
        HybridEvolutionTree tree;
        const auto a = add_pending_checkpoint(tree, boxed({
                                                        {0, 0},
                                                        {0, 0}
        }));
        const auto b = add_pending_checkpoint(tree, boxed({
                                                        {75, 75},
                                                        { 0,  0}
        }));
        const auto c = add_pending_checkpoint(tree, boxed({
                                                        {150, 150},
                                                        {  0,   0}
        }),
                                              b.trajectory);
        const Scalings scales(1.0_x, {
                                         {RealVariable("x"), 100.0_x}
        });
        const std::vector<HybridEvolutionCheckpointId> pending{a, b, c};
        ARIADNE_TESTING_PROVE(ProximityRecombinationPolicy(1.0).select(tree, pending).empty());
        const ProximityRecombinationPolicy local(1.0, 8u, 40.0, scales);
        const auto groups = local.select(tree, pending);
        ARIADNE_TESTING_EXPECT_EQUAL(groups.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(groups.front().size(), 2u);
        ARIADNE_TESTING_PROVE(groups.front()[0] == a);
        ARIADNE_TESTING_PROVE(groups.front()[1] == b);
        // 0 -> .75 -> 1.5 is a chain, but its full normalized diameter exceeds one.
        ARIADNE_TESTING_PROVE(ProximityRecombinationPolicy(1.0, 8u, 30.0, scales).select(tree, pending).empty());
        // Geometric grouping is independent of immediate genealogy parents.
        const auto unrelated = local.select(tree, {b, c});
        ARIADNE_TESTING_EXPECT_EQUAL(unrelated.size(), 1u);
        const auto far_y = add_pending_checkpoint(tree, boxed({
                                                            {0, 0},
                                                            {2, 2}
        }));
        ARIADNE_TESTING_PROVE(local.select(tree, {a, far_y}).empty());
        const auto boundary_left = add_pending_checkpoint(tree, boxed({
                                                                    {99, 99},
                                                                    { 0,  0}
        }));
        const auto boundary_right = add_pending_checkpoint(tree, boxed({
                                                                     {101, 101},
                                                                     {  0,   0}
        }));
        ARIADNE_TESTING_EXPECT_EQUAL(local.select(tree, {boundary_left, boundary_right}).size(), 1u);

        HybridEvolutionTree wide;
        const auto first = add_pending_checkpoint(wide, boxed({
                                                            {-50, 50},
                                                            {  0,  0}
        }));
        const auto second = add_pending_checkpoint(wide, boxed({
                                                             {25, 125},
                                                             { 0,   0}
        }));
        ARIADNE_TESTING_PROVE(ProximityRecombinationPolicy(1.0, 8u, 60.0, scales).select(wide, {first, second}).empty());
        ARIADNE_TESTING_EXPECT_EQUAL(ProximityRecombinationPolicy(1.0, 8u, 90.0, scales).select(wide, {first, second}).size(), 1u);

        HybridEvolutionTree dense;
        std::vector<HybridEvolutionCheckpointId> coincident;
        for (unsigned i = 0; i != 10u; ++i)
            coincident.push_back(add_pending_checkpoint(dense, initial));
        const auto bounded = ProximityRecombinationPolicy(1.0, 3u).select(dense, coincident);
        ARIADNE_TESTING_EXPECT_EQUAL(bounded.size(), 3u);
        for (const auto& group: bounded)
            ARIADNE_TESTING_EXPECT_EQUAL(group.size(), 3u);
    }

    void test_proximity_uses_reduced_domains_and_preserves_compatibility() {
        const auto initial = seed();
        auto constrained = HybridEnclosure(
            HybridExactBox(initial.location(), initial.state_space(), ExactBoxType({
                                                                          {0, 100},
                                                                          {0,   0}
        })),
            initial.configuration());
        const auto coordinate = EffectiveVectorMultivariateFunction::identity(2u);
        constrained.new_state_constraint(DiscreteEvent("bound"), ValidatedConstraint(0, coordinate[0], 1));
        constrained.reduce();
        HybridEvolutionTree tree;
        const auto a = add_pending_checkpoint(tree, constrained);
        // Disjoint clock intervals are allowed; spatial distance is measured on the reduced set.
        const auto b = add_pending_checkpoint(tree, marked(initial, 30));
        const ProximityRecombinationPolicy policy(1.0, 8u, 1.0);
        ARIADNE_TESTING_EXPECT_EQUAL(policy.select(tree, {a, b}).size(), 1u);
        const auto other = add_pending_checkpoint(tree, jumped(initial, DiscreteLocation(2), DiscreteEvent("jump"), 30));
        ARIADNE_TESTING_PROVE(policy.select(tree, {a, other}).empty());
        const auto different_event = add_pending_checkpoint(tree, jumped(initial, initial.location(), DiscreteEvent("jump"), 30));
        ARIADNE_TESTING_PROVE(policy.select(tree, {a, different_event}).empty());
        tree.object(b.trajectory).entries[0].work.phase = HybridEvolutionPhase::EVOLVE_CONTINUOUSLY;
        ARIADNE_TESTING_PROVE(policy.select(tree, {a, b}).empty());
    }

    void test_proximity_waiting_replacements_evolve_before_remerging() {
        ConcurrencyScope concurrency;
        const auto initial = seed();
        const DiscreteLocation target(2);
        const DiscreteEvent jump("jump");
        TreeHybridEvolutionBackend backend(
            initial, LocationCheckpointAdmissionPolicy{}, ProximityRecombinationPolicy(1.0, 2u, 1.0),
            EveryNStepsRecombinationTriggerPolicy(2u), ReferenceEnclosureRecombiner{}, EnclosureRecombinationContext(1.0_x));
        backend.set_semantics(Semantics::UPPER);
        unsigned emitted = 0u, executed = 0u;
        bool merged_before_completion = false;
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            if (item.enclosure.location() == target) {
                ARIADNE_TESTING_EXPECT_EQUAL(emitted, 8u);
                ++executed;
                result.record_final(item.enclosure);
            } else {
                if (label(item) == 2) {
                    for (const auto node: backend.tree().descendants(backend.tree().root().value()))
                        merged_before_completion = merged_before_completion || !backend.tree().metadata(node).recombined_checkpoints.empty();
                }
                ++emitted;
                result.transition_to(jumped(initial, target, jump, 100), jump, EventKind::PERMISSIVE);
                if (emitted < 8u)
                    result.continue_with(marked(initial, static_cast<int>(emitted)));
            }
            return result;
        });
        ARIADNE_TESTING_PROVE(merged_before_completion);
        ARIADNE_TESTING_EXPECT_EQUAL(executed, 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 4u);
        // A replacement seed cannot be merged again while waiting, but a continuation can.
        HybridEvolutionTree tree;
        const auto a = add_pending_checkpoint(tree, initial);
        const auto b = add_pending_checkpoint(tree, initial);
        tree.metadata(a.trajectory).recombined_checkpoints.push_back(b);
        const ProximityRecombinationPolicy policy(1.0);
        ARIADNE_TESTING_PROVE(policy.select(tree, {a, b}).empty());
        tree.object(a.trajectory).entries.push_back(tree.object(a.trajectory).entries[0]);
        ARIADNE_TESTING_EXPECT_EQUAL(policy.select(tree, {
                                                             {a.trajectory, 1u},
                                                             b
        })
                                         .size(),
                                     1u);
    }

    void test_proximity_selects_disjoint_local_pairs_in_a_large_frontier() {
        const auto initial = seed();
        HybridEvolutionTree tree;
        std::vector<HybridEvolutionCheckpointId> pending;
        for (int i = 0; i != 1024; ++i) {
            const ExactDouble coordinate(4 * i);
            const auto enclosure = HybridEnclosure(
                HybridExactBox(initial.location(), initial.state_space(), ExactBoxType({
                                                                              {coordinate, coordinate},
                                                                              {         0,          0}
            })),
                initial.configuration());
            pending.push_back(add_pending_checkpoint(tree, enclosure));
            pending.push_back(add_pending_checkpoint(tree, enclosure));
        }
        const auto groups = ProximityRecombinationPolicy(1.0, 8u, 1.0).select(tree, pending);
        ARIADNE_TESTING_EXPECT_EQUAL(groups.size(), 1024u);
        for (const auto& group: groups)
            ARIADNE_TESTING_EXPECT_EQUAL(group.size(), 2u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(ProximityRecombinationPolicy(0.0));
        ARIADNE_TESTING_EXPECT_EXCEPTION(ProximityRecombinationPolicy(-1.0));
        ARIADNE_TESTING_EXPECT_EXCEPTION(ProximityRecombinationPolicy(1.0, 1u));
        ARIADNE_TESTING_EXPECT_EXCEPTION(ProximityRecombinationPolicy(1.0, 8u, -1.0));
        ARIADNE_TESTING_EXPECT_EXCEPTION(ProximityRecombinationPolicy{std::numeric_limits<double>::infinity()});
    }

    void test_all_compatible_siblings_merge_nonadjacent_unbounded_groups() {
        const auto initial = seed();
        const DiscreteEvent jump("jump"), other("other"), isolated("isolated");
        const auto advance = [&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            if (label(item) == 0) {
                // A and B successors are interleaved intentionally. The
                // last source is incompatible and must stay independent.
                result.transition_to(marked(initial, 10), jump, EventKind::URGENT);
                result.transition_to(marked(initial, 10), other, EventKind::PERMISSIVE);
                result.transition_to(marked(initial, 10), jump, EventKind::URGENT);
                result.transition_to(marked(initial, 10), isolated, EventKind::URGENT);
                result.transition_to(marked(initial, 10), other, EventKind::PERMISSIVE);
                result.transition_to(marked(initial, 10), jump, EventKind::URGENT);
            } else if (label(item) == 10) {
                result.record_final(item.enclosure);
            } else {
                throw std::logic_error("Unexpected work item");
            }
            return result;
        };

        TreeHybridEvolutionBackend baseline(initial);
        baseline.run(advance);
        ARIADNE_TESTING_EXPECT_EQUAL(baseline.release_orbit().final().size(), 6u);

        TreeHybridEvolutionBackend merged(
            initial, ParentCheckpointAdmissionPolicy{}, AllCompatibleRecombinationPolicy(HybridRecombinationScope::SIBLINGS),
            BeforeAdmissionRecombinationTriggerPolicy{});
        merged.set_semantics(Semantics::UPPER);
        merged.run(advance);
        const auto root = merged.tree().root().value();
        const auto children = merged.tree().children(root);
        ARIADNE_TESTING_EXPECT_EQUAL(children.size(), 8u);
        ARIADNE_TESTING_EXPECT_EQUAL(merged.tree().size(), 9u);
        for (std::size_t i = 0u; i != 6u; ++i) {
            const bool should_retire = i != 3u;
            ARIADNE_TESTING_PROVE(
                static_cast<bool>(!merged.tree().object(children[i]).retired_entries.empty()) == should_retire);
            if (should_retire)
                ARIADNE_TESTING_PROVE(!merged.tree().object(children[i]).entries.front().result);
        }
        const auto& first = merged.tree().metadata(children[6]);
        const auto& second = merged.tree().metadata(children[7]);
        ARIADNE_TESTING_EXPECT_EQUAL(first.recombined_checkpoints.size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(second.recombined_checkpoints.size(), 2u);
        ARIADNE_TESTING_PROVE(first.recombined_checkpoints[0].trajectory == children[0]);
        ARIADNE_TESTING_PROVE(first.recombined_checkpoints[1].trajectory == children[2]);
        ARIADNE_TESTING_PROVE(first.recombined_checkpoints[2].trajectory == children[5]);
        ARIADNE_TESTING_PROVE(second.recombined_checkpoints[0].trajectory == children[1]);
        ARIADNE_TESTING_PROVE(second.recombined_checkpoints[1].trajectory == children[4]);
        ARIADNE_TESTING_PROVE(first.complete && second.complete);
        ARIADNE_TESTING_EXPECT_EQUAL(merged.release_orbit().final().size(), 3u);
    }

    void test_all_compatible_siblings_can_resplit_after_unbounded_merge() {
        const auto initial = seed();
        const DiscreteEvent jump("jump");
        std::size_t replacements_executed = 0u;
        TreeHybridEvolutionBackend merged(
            initial, ParentCheckpointAdmissionPolicy{}, AllCompatibleRecombinationPolicy(HybridRecombinationScope::SIBLINGS));
        merged.set_semantics(Semantics::UPPER);
        merged.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            HybridEvolutionStepResult result;
            switch (label(item)) {
                case 0:
                    result.transition_to(marked(initial, 10), jump, EventKind::URGENT);
                    result.transition_to(marked(initial, 10), jump, EventKind::URGENT);
                    result.transition_to(marked(initial, 10), jump, EventKind::URGENT);
                    break;
                case 10:
                    ++replacements_executed;
                    // The numerical evolver is permitted to split a merged
                    // enclosure again; all policies now record explicit split trajectories.
                    result.continue_with(marked(initial, 11),
                                         HybridEvolutionSuccessorKind::SUBDIVISION);
                    result.continue_with(marked(initial, 12),
                                         HybridEvolutionSuccessorKind::SUBDIVISION);
                    break;
                case 11:
                case 12:
                    result.record_final(item.enclosure);
                    break;
                default:
                    throw std::logic_error("Unexpected work item after recombination");
            }
            return result;
        });
        const auto root = merged.tree().root().value();
        const auto children = merged.tree().children(root);
        ARIADNE_TESTING_EXPECT_EQUAL(children.size(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            merged.tree().metadata(children[3]).recombined_checkpoints.size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(replacements_executed, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(merged.tree().object(children[3]).entries.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(merged.tree().children(children[3]).size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(merged.release_orbit().final().size(), 2u);
    }

    void test_waiting_keeps_early_replacements_held_until_parent_completion() {
        const auto initial = seed();
        const DiscreteEvent first("first"), second("second");
        TreeHybridEvolutionBackend backend(
            initial, ParentCheckpointAdmissionPolicy{}, AllCompatibleRecombinationPolicy(HybridRecombinationScope::SIBLINGS),
            AfterEachStepRecombinationTriggerPolicy{});
        backend.set_semantics(Semantics::UPPER);
        const auto root = backend.tree().root().value();
        std::vector<int> observed;
        backend.run([&](const HybridEvolutionWorkItem& item, const HybridEvolutionStatistics&) {
            const auto marker = label(item);
            observed.push_back(marker);
            HybridEvolutionStepResult result;
            if (marker == 0) {
                result.transition_to(marked(initial, 100), first, EventKind::URGENT);
                result.transition_to(marked(initial, 100), first, EventKind::URGENT);
                result.continue_with(marked(initial, 1));
            } else if (marker == 1) {
                // The first pair has already been merged at AFTER_STEP, but
                // its replacement must still be held while the parent runs.
                const auto children = backend.tree().children(root);
                ARIADNE_TESTING_EXPECT_EQUAL(children.size(), 3u);
                ARIADNE_TESTING_PROVE(!backend.tree().metadata(root).complete);
                ARIADNE_TESTING_PROVE(!backend.tree().object(children[0]).retired_entries.empty());
                ARIADNE_TESTING_PROVE(!backend.tree().object(children[1]).retired_entries.empty());
                ARIADNE_TESTING_PROVE(!backend.tree().object(children[2]).entries.front().result);
                result.transition_to(marked(initial, 200), second, EventKind::PERMISSIVE);
                result.transition_to(marked(initial, 200), second, EventKind::PERMISSIVE);
                result.continue_with(marked(initial, 2));
            } else if (marker == 2) {
                // Both accepted replacements are still withheld.
                const auto children = backend.tree().children(root);
                ARIADNE_TESTING_EXPECT_EQUAL(children.size(), 6u);
                ARIADNE_TESTING_PROVE(!backend.tree().metadata(root).complete);
                ARIADNE_TESTING_PROVE(!backend.tree().object(children[2]).entries.front().result);
                ARIADNE_TESTING_PROVE(!backend.tree().object(children[5]).entries.front().result);
                result.record_final(item.enclosure);
            } else if (marker == 100 || marker == 200) {
                ARIADNE_TESTING_PROVE(backend.tree().metadata(root).complete);
                result.record_final(item.enclosure);
            } else {
                throw std::logic_error("Unexpected waiting-queue checkpoint");
            }
            return result;
        });
        const auto children = backend.tree().children(root);
        ARIADNE_TESTING_PROVE(observed == std::vector<int>({0, 1, 2, 100, 200}));
        ARIADNE_TESTING_EXPECT_EQUAL(children.size(), 6u);
        for (auto i: {0u, 1u, 3u, 4u})
            ARIADNE_TESTING_PROVE(!backend.tree().object(children[i]).retired_entries.empty());
        for (auto i: {2u, 5u}) {
            ARIADNE_TESTING_PROVE(backend.tree().metadata(children[i]).complete);
            ARIADNE_TESTING_EXPECT_EQUAL(
                backend.tree().metadata(children[i]).recombined_checkpoints.size(), 2u);
        }
        ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), 3u);
    }

    void test_recombination_policy_rejects_lower_semantics_and_invalid_settings() {
        const auto initial = seed();
        ARIADNE_TESTING_EXPECT_EXCEPTION(WindowRecombinationPolicy(1u, 1.0, HybridRecombinationScope::SIBLINGS));
        ARIADNE_TESTING_EXPECT_EXCEPTION(WindowRecombinationPolicy(2u, -1.0, HybridRecombinationScope::SIBLINGS));
        TreeHybridEvolutionBackend lower(initial, ParentCheckpointAdmissionPolicy{}, WindowRecombinationPolicy(8u, 0.125, HybridRecombinationScope::SIBLINGS));
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
        const auto verify = [&](const HybridRecombinationGroupingPolicy& grouping,
                                const HybridRecombinationTriggerPolicy& timing,
                                bool after_first, bool after_second, bool after_complete,
                                std::size_t expected_final_count) {
            TreeHybridEvolutionBackend backend(initial, ParentCheckpointAdmissionPolicy{}, grouping, timing);
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
                    const bool retired = !backend.tree().object(children[0]).retired_entries.empty();
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
            ARIADNE_TESTING_EXPECT_EQUAL(!backend.tree().object(children[0]).retired_entries.empty(),
                                         after_complete);
            ARIADNE_TESTING_EXPECT_EQUAL(backend.release_orbit().final().size(), expected_final_count);
        };
        const WindowRecombinationPolicy window(2u, 0.0, HybridRecombinationScope::SIBLINGS);
        verify(NoRecombinationPolicy{}, AfterEachStepRecombinationTriggerPolicy{},
               false, false, false, 3u);
        verify(window, BeforeAdmissionRecombinationTriggerPolicy{},
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
        TreeHybridEvolutionBackendFactory factory(ParentCheckpointAdmissionPolicy{},
                                                  WindowRecombinationPolicy(2u, 0.0, HybridRecombinationScope::SIBLINGS),
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
        TreeHybridEvolutionBackend baseline(initial, ImmediateCheckpointAdmissionPolicy{}, NoRecombinationPolicy{});
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
            initial, ImmediateCheckpointAdmissionPolicy{}, WindowRecombinationPolicy(2u, 0.0),
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
            TreeHybridEvolutionBackend baseline(initial, ImmediateCheckpointAdmissionPolicy{}, NoRecombinationPolicy{});
            baseline.run(advance);
            ARIADNE_TESTING_EXPECT_EQUAL(baseline.release_orbit().final().size(), 2u);

            TreeHybridEvolutionBackend merged(
                initial, ImmediateCheckpointAdmissionPolicy{}, WindowRecombinationPolicy(2u, 0.0),
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
        ARIADNE_TESTING_EXPECT_EXCEPTION(WindowRecombinationPolicy(1u, 0.125));
        ARIADNE_TESTING_EXPECT_EXCEPTION(WindowRecombinationPolicy(2u, -1.0));
        TreeHybridEvolutionBackend lower(initial, ImmediateCheckpointAdmissionPolicy{}, WindowRecombinationPolicy{});
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
        const auto pieces = backend.tree().children(root);
        ARIADNE_TESTING_EXPECT_EQUAL(pieces.size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(backend.tree().object(root).entries.size(), 1u);
        for (std::size_t i = 0; i != pieces.size(); ++i) {
            const auto origin = backend.tree().metadata(pieces[i]).origin.value();
            ARIADNE_TESTING_EXPECT_EQUAL(origin.entry, 0u);
            ARIADNE_TESTING_EXPECT_EQUAL(origin.successor, i);
            ARIADNE_TESTING_PROVE(backend.tree().metadata(pieces[i]).complete);
        }
        const auto& first_entries = backend.tree().object(pieces[0]).entries;
        ARIADNE_TESTING_EXPECT_EQUAL(first_entries.size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(first_entries.at(1).origin->entry, 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(backend.tree().children(pieces[1]).size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(backend.tree().size(), 4u);
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
            HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
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
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
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
        HybridContinuousTimeEvolver<ValidatedTag> evolver(system);
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
        HybridContinuousTimeEvolver<ValidatedTag> constrained(system);
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
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_progress_recombines_waiting_time_without_inventing_evolution);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_progress_counts_parallel_in_flight_work);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_periodic_ready_recombination_does_not_flush_on_each_admission);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_sibling_window_ignores_interleaved_unrelated_parents);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_parent_admission_releases_completed_families_without_global_barrier);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_parent_admission_preserves_all_dependencies_of_cross_parent_replacements);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_waiting_admission_rejects_continuous_only_quality_gate);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_location_batch_merges_different_parents_after_all_producers_finish);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_location_early_replacements_stay_held_and_absorb_later_arrivals);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_location_batches_keep_destinations_separate_and_reenter_later);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_parallel_location_batch_waits_for_in_flight_producers);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_location_factory_preserves_real_orbits_without_recombination);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_parent_finishes_before_distinct_derived_trajectories);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_parallel_parent_first_last_publisher_merges_held_children);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_parallel_frontier_executes_independent_split_checkpoints);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_only_certified_single_urgent_transition_continues_the_node);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_noop_and_window_recombination_policies);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_proximity_bounds_normalized_diameter_radius_and_cardinality);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_proximity_uses_reduced_domains_and_preserves_compatibility);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_proximity_waiting_replacements_evolve_before_remerging);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_proximity_selects_disjoint_local_pairs_in_a_large_frontier);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_all_compatible_siblings_merge_nonadjacent_unbounded_groups);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_all_compatible_siblings_can_resplit_after_unbounded_merge);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_waiting_keeps_early_replacements_held_until_parent_completion);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_recombination_policy_rejects_lower_semantics_and_invalid_settings);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_recombination_evaluation_triggers);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_numerical_recombination_strategy_is_injected);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_frontier_recombines_continuous_split_branches_after_evolution);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_frontier_recombines_cross_trajectory_transitions);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_flat_and_frontier_orbits_match_on_real_hybrid_transitions);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_flat_and_frontier_orbits_match_on_real_splitting_and_invariants);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_real_frontier_recombination_after_continuous_subdivision);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_real_frontier_recombination_radius_rejection_preserves_work);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_real_frontier_enabled_across_urgent_and_permissive_resets);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_real_frontier_recombines_only_after_independent_continuous_steps);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_real_post_reset_recombination_of_split_trajectories);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_frontier_recombination_requires_progress_after_resplitting);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_bounded_frontier_wait_recovers_late_merge_partners);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_bounded_frontier_wait_does_not_starve_unmatched_work);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_radius_trend_gate_accepts_only_profitable_continuations);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_split_pieces_are_parallel_and_completion_waits_for_all_of_them);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_failure_preserves_incomplete_history_without_dispatching_children);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_flat_and_tree_orbits_match_for_urgent_and_permissive_systems);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_flat_and_tree_orbit_counts_match_at_time_and_discrete_limits);
    ARIADNE_TESTING_ADD_TEST(TestHybridEvolutionTree, test_flat_and_tree_orbit_counts_match_with_initial_jumps_and_splitting);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
