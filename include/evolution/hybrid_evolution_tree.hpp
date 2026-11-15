/***************************************************************************
 *            evolution/hybrid_evolution_tree.hpp
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

/** @file evolution/hybrid_evolution_tree.hpp
 *  @brief Trajectory histories and a serial parent-first hybrid evolution backend.
 */

#ifndef ARIADNE_HYBRID_EVOLUTION_TREE_HPP
#define ARIADNE_HYBRID_EVOLUTION_TREE_HPP

#include "evolution/enclosure_recombiner.hpp"
#include "evolution/hybrid_evolution_backend.hpp"
#include "tree/derivation_tree.hpp"

namespace Ariadne {

/** @brief Source entry and successor ordinal identifying an append-only derivation.
 *  @details For an entry this refers within its own trajectory; for node metadata
 *  it refers to the parent trajectory. No enclosure is duplicated in this link.
 */
struct HybridEvolutionEntryOrigin {
    std::size_t entry;     //!< Index of the computation that generated this checkpoint.
    std::size_t successor; //!< Index in that computation's StepResult::successors.
};

/** @brief One checkpoint and its local computation, distinct from its Orbit outputs.
 *  @details Entries are ordered by creation. Equal predecessors identify parallel
 *  pieces; adjacency in the vector does not imply sequential physical evolution.
 */
struct HybridEvolutionEntry {
    HybridEvolutionWorkItem work;                     //!< Computational checkpoint, including its entry-processing phase.
    std::optional<HybridEvolutionEntryOrigin> origin; //!< Same-trajectory predecessor; absent for the seed.
    std::optional<HybridEvolutionStepResult> result;  //!< Absent until this checkpoint is processed successfully.
};

/** @brief Stable checkpoint identity: one entry of a genealogy trajectory. */
struct HybridEvolutionCheckpointId {
    DerivationTreeNodeId trajectory;
    std::size_t entry;

    friend bool operator==(const HybridEvolutionCheckpointId&,
                           const HybridEvolutionCheckpointId&) = default;
};

/** @brief One trajectory's ordered history; a trajectory is one generic tree node. */
struct HybridEvolutionTrajectory {
    std::vector<HybridEvolutionEntry> entries; //!< Seed followed by continuations and numerical pieces.
    std::vector<std::size_t> retired_entries;   //!< Unprocessed checkpoints replaced by certified frontier recombination.
};

/** @brief Local trajectory completion and provenance of its birth. */
struct HybridEvolutionTrajectoryMetadata {
    Bool complete = false;                            //!< All local pieces processed and all derived children published; independent of descendants.
    std::optional<HybridEvolutionEntryOrigin> origin; //!< Parent entry/successor producing the seed; absent at the root and at a recombined replacement.
    bool retired_for_recombination = false;            //!< Original child kept as provenance, never dispatched or projected.
    std::vector<DerivationTreeNodeId> recombined_sources; //!< Original sibling IDs represented by a replacement.
    std::vector<HybridEvolutionCheckpointId> recombined_checkpoints; //!< Exact entry-level provenance of a frontier replacement.
    std::size_t pending_entries = 0u; //!< Number of still-pending work items in frontier mode.
    bool whole_trajectory = true; //!< Whole-source urgent continuation eligibility in frontier mode.
    std::optional<double> recombination_time_barrier; //!< Upper time bound of the most recent frontier merge along this lineage; block re-merging until time has advanced strictly beyond it.
};

/** @brief Generic two-parameter genealogy with Evolution-owned objects and metadata. */
using HybridEvolutionTree = DerivationTree<HybridEvolutionTrajectory, HybridEvolutionTrajectoryMetadata>;

/** @brief Points at which a scheduler may reconsider pending recombination. */
enum class HybridRecombinationEvaluationPoint : std::uint8_t {
    AFTER_STEP,       //!< A numerical work item has completed and its results have been committed.
    PARENT_COMPLETE   //!< All work of this trajectory has completed.
};

/** @brief Immutable policy specifying when pending children should be inspected.
 *  @details Called with the completed step count of the parent trajectory.
 *  Scheduling decisions must be stateless and must not mutate the history.
 */
class HybridRecombinationTriggerPolicy {
  public:
    using NodeId = HybridEvolutionTree::NodeId;
    virtual ~HybridRecombinationTriggerPolicy() = default;
    /** @brief Clone an independently owned policy with the same configuration. */
    virtual HybridRecombinationTriggerPolicy* clone() const = 0;
    /** @brief Decide whether to invoke grouping at this committed work boundary. */
    virtual bool should_evaluate(const HybridEvolutionTree& tree, NodeId parent,
                                 HybridRecombinationEvaluationPoint point,
                                 std::size_t completed_steps) const = 0;
};

/** @brief Inspect children only when their parent trajectory is complete. */
class OnParentCompletionRecombinationTriggerPolicy final : public HybridRecombinationTriggerPolicy {
  public:
    OnParentCompletionRecombinationTriggerPolicy* clone() const override;
    bool should_evaluate(const HybridEvolutionTree&, NodeId,
                         HybridRecombinationEvaluationPoint point, std::size_t) const override;
};

/** @brief Inspect children after every committed numerical step, continuous or discrete. */
class AfterEachStepRecombinationTriggerPolicy final : public HybridRecombinationTriggerPolicy {
  public:
    AfterEachStepRecombinationTriggerPolicy* clone() const override;
    bool should_evaluate(const HybridEvolutionTree&, NodeId,
                         HybridRecombinationEvaluationPoint point, std::size_t) const override;
};

/** @brief Inspect children after each N steps and flush any remainder at parent completion. */
class EveryNStepsRecombinationTriggerPolicy final : public HybridRecombinationTriggerPolicy {
    std::size_t _period;

  public:
    /** @brief Configure the positive number of completed steps between inspections. */
    explicit EveryNStepsRecombinationTriggerPolicy(std::size_t period);
    EveryNStepsRecombinationTriggerPolicy* clone() const override;
    bool should_evaluate(const HybridEvolutionTree&, NodeId,
                         HybridRecombinationEvaluationPoint point,
                         std::size_t completed_steps) const override;
};

/** @brief Policy for selecting disjoint, unevolved sibling groups at evaluation points.
 *  @details The trigger chooses when this selection runs; this policy decides
 *  which pending children are candidates. Numerical construction is delegated
 *  to a separately injected enclosure recombiner. Policies must not mutate the tree.
 */
class HybridRecombinationSchedulingPolicy {
  public:
    using NodeId = HybridEvolutionTree::NodeId;
    virtual ~HybridRecombinationSchedulingPolicy() = default;
    /** @brief Clone the immutable selection settings. */
    virtual HybridRecombinationSchedulingPolicy* clone() const = 0;
    /** @brief Select disjoint groups of unevolved children at the requested boundary. */
    virtual std::vector<std::vector<NodeId>> select(const HybridEvolutionTree& tree, NodeId parent) const = 0;
    /** @brief Whether this policy can request upper-enclosure recombination. */
    virtual bool enables_recombination() const = 0;
};

/** @brief Preserve the existing parent-first schedule without recombination. */
class NoRecombinationSchedulingPolicy final : public HybridRecombinationSchedulingPolicy {
  public:
    NoRecombinationSchedulingPolicy* clone() const override;
    std::vector<std::vector<NodeId>> select(const HybridEvolutionTree& tree, NodeId parent) const override;
    bool enables_recombination() const override { return false; }
};

/** @brief Greedily group consecutive compatible siblings within a birth-time window.
 *  @details Groups have at least two and at most max_group_size sources. Different
 *  events, target locations, state/auxiliary spaces, phases, event histories and
 *  symbolic auxiliary maps are not combined. A group is never split across
 *  unrelated siblings. Upper semantics is required for numerical recombination.
 */
class SiblingWindowRecombinationSchedulingPolicy final : public HybridRecombinationSchedulingPolicy {
    std::size_t _max_group_size;
    double _max_time_spread;

  public:
    /** @brief Configure a bounded batch size and birth-time interval width. */
    explicit SiblingWindowRecombinationSchedulingPolicy(std::size_t max_group_size = 8u,
                                                         double max_time_spread = 0.125);
    SiblingWindowRecombinationSchedulingPolicy* clone() const override;
    std::vector<std::vector<NodeId>> select(const HybridEvolutionTree& tree, NodeId parent) const override;
    bool enables_recombination() const override { return true; }
};

/** @brief Select arbitrary compatible, unevolved checkpoints from the global frontier.
 *  @details The policy sees all queued checkpoints, not merely siblings; each
 *  selected checkpoint must be present, unprocessed and unique. The backend
 *  validates numerical/discrete compatibility independently of the policy.
 */
class HybridFrontierRecombinationPolicy {
  public:
    using CheckpointId = HybridEvolutionCheckpointId;
    virtual ~HybridFrontierRecombinationPolicy() = default;
    /** @brief Clone immutable grouping configuration. */
    virtual HybridFrontierRecombinationPolicy* clone() const = 0;
    /** @brief Propose disjoint groups of pending checkpoint identities. */
    virtual std::vector<std::vector<CheckpointId>> select(
        const HybridEvolutionTree& history, const std::vector<CheckpointId>& pending) const = 0;
    /** @brief True when numerical upper-covering recombination may be attempted. */
    virtual bool enables_recombination() const = 0;
};

/** @brief Frontier scheduler without numerical recombination. */
class NoFrontierRecombinationPolicy final : public HybridFrontierRecombinationPolicy {
  public:
    NoFrontierRecombinationPolicy* clone() const override;
    std::vector<std::vector<CheckpointId>> select(
        const HybridEvolutionTree&, const std::vector<CheckpointId>&) const override;
    bool enables_recombination() const override { return false; }
};

/** @brief Greedy global time-window grouping across split branches and reset trajectories.
 *  @details Enforces a maximum time span and batch size; nonadjacent pending
 *  checkpoints can be grouped. Backend rechecks compatibility and coverage.
 *  The default permits only overlapping time intervals (including identical
 *  point times); no arbitrary time shift is introduced.
 */
class FrontierWindowRecombinationPolicy final : public HybridFrontierRecombinationPolicy {
    std::size_t _max_group_size;
    double _max_time_spread;

  public:
    /** @brief Configure number of sources and maximum common time-window span. */
    explicit FrontierWindowRecombinationPolicy(std::size_t max_group_size = 8u,
                                               double max_time_spread = 0.125);
    FrontierWindowRecombinationPolicy* clone() const override;
    std::vector<std::vector<CheckpointId>> select(
        const HybridEvolutionTree& history, const std::vector<CheckpointId>& pending) const override;
    bool enables_recombination() const override { return true; }
};

/** @brief Serial backend with parent-first and opt-in checkpoint-frontier schedules.
 *  @details The default/legacy grouping path finishes a parent before its derived
 *  trajectories and preserves the pre-existing order. The opt-in frontier path
 *  schedules individual work items across trajectories; subdivision creates
 *  explicit child trajectories, and completed branches can meet again through
 *  time-compatible pending checkpoints, including after continuous evolution.
 *  No processed result is removed or recombined retroactively. Original pending
 *  work is retired only after certified full coverage of the selected sources.
 *  Numerical recombination requires upper semantics in both schedules.
 */
class TreeHybridEvolutionBackend final : public HybridEvolutionBackendInterface {
    class Impl;
    std::unique_ptr<Impl> _impl;

  public:
    /** @brief Create the root trajectory with its unprocessed initial checkpoint. */
    explicit TreeHybridEvolutionBackend(
        const HybridEnclosure& initial,
        const HybridRecombinationSchedulingPolicy& grouping = NoRecombinationSchedulingPolicy{},
        const HybridRecombinationTriggerPolicy& trigger = OnParentCompletionRecombinationTriggerPolicy{},
        const EnclosureRecombinerInterface& recombiner = ReferenceEnclosureRecombiner{},
        const EnclosureRecombinationContext& context = {});
    /** @brief Use a global checkpoint frontier, with configurable time and grouping policies.
     *  @details Unlike the parent-first constructor, this permits continuous
     *  branches from earlier subdivisions to be recombined before they terminate.
     */
    explicit TreeHybridEvolutionBackend(
        const HybridEnclosure& initial,
        const HybridFrontierRecombinationPolicy& grouping,
        const HybridRecombinationTriggerPolicy& trigger = AfterEachStepRecombinationTriggerPolicy{},
        const EnclosureRecombinerInterface& recombiner = ReferenceEnclosureRecombiner{},
        const EnclosureRecombinationContext& context = {});
    /** @brief Release the owned trajectory history. */
    ~TreeHybridEvolutionBackend() override;
    /** @brief Disallow copying the owner and its trajectory identities. */
    TreeHybridEvolutionBackend(const TreeHybridEvolutionBackend&) = delete;
    /** @brief Disallow copying the owner and its trajectory identities. */
    TreeHybridEvolutionBackend& operator=(const TreeHybridEvolutionBackend&) = delete;
    /** @brief Run once, processing each trajectory's pieces before releasing its children.
     *  @note Callback failures leave an inspectable incomplete history and prohibit Orbit release.
     */
    Void set_semantics(Semantics semantics) override;
    Void run(const StepFunction& advance) override;
    /** @brief Materialize the conventional Orbit once, retaining history for const inspection. */
    Orbit<HybridEnclosure> release_orbit() override;
    /** @brief Read genealogy, checkpoints, step results and local completion.
     *  @note Access is sequential. A failed run can contain processed and unprocessed entries.
     *  The reference is valid until this backend is destroyed; growth can invalidate entry references.
     */
    const HybridEvolutionTree& tree() const noexcept;
};

/** @brief Factory selecting the serial trajectory-tree backend for ordinary orbit() calls.
 *  @details The independent grouping, timing, numerical recombiner and radius
 *  context settings are cloned or copied for each run.
 */
class TreeHybridEvolutionBackendFactory : public HybridEvolutionBackendFactoryInterface {
    std::shared_ptr<const HybridRecombinationSchedulingPolicy> _grouping;
    std::shared_ptr<const HybridFrontierRecombinationPolicy> _frontier_grouping;
    std::shared_ptr<const HybridRecombinationTriggerPolicy> _trigger;
    std::shared_ptr<const EnclosureRecombinerInterface> _recombiner;
    EnclosureRecombinationContext _context;

  public:
    /** @brief Use default no-op grouping and parent-completion evaluation. */
    TreeHybridEvolutionBackendFactory();
    /** @brief Independently choose grouping, evaluation timing, numerical method and radius context. */
    explicit TreeHybridEvolutionBackendFactory(
        const HybridRecombinationSchedulingPolicy& grouping,
        const HybridRecombinationTriggerPolicy& trigger = OnParentCompletionRecombinationTriggerPolicy{},
        const EnclosureRecombinerInterface& recombiner = ReferenceEnclosureRecombiner{},
        const EnclosureRecombinationContext& context = {});
    /** @brief Opt into interleaved checkpoint evolution and cross-trajectory recombination. */
    explicit TreeHybridEvolutionBackendFactory(
        const HybridFrontierRecombinationPolicy& grouping,
        const HybridRecombinationTriggerPolicy& trigger = AfterEachStepRecombinationTriggerPolicy{},
        const EnclosureRecombinerInterface& recombiner = ReferenceEnclosureRecombiner{},
        const EnclosureRecombinationContext& context = {});
    /** @brief Clone the immutable factory settings without copying any run history. */
    TreeHybridEvolutionBackendFactory* clone() const override;
    /** @brief Create independent trajectory history; the numerical callback enforces the horizon. */
    std::unique_ptr<HybridEvolutionBackendInterface> create(const HybridEnclosure& initial, const Real& final_time) const override;
};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLUTION_TREE_HPP
