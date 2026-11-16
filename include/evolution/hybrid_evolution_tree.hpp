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
 *  @brief Trajectory histories and policy-driven concurrent hybrid evolution.
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
    DerivationTreeNodeId trajectory; //!< Stable trajectory node identity.
    std::size_t entry;               //!< Index in the trajectory history.

    /** @brief Compare stable trajectory and entry identities. */
    friend bool operator==(const HybridEvolutionCheckpointId&,
                           const HybridEvolutionCheckpointId&) = default;
};

/** @brief One trajectory's ordered history; a trajectory is one generic tree node. */
struct HybridEvolutionTrajectory {
    std::vector<HybridEvolutionEntry> entries; //!< Seed followed by continuations and numerical pieces.
    std::vector<std::size_t> retired_entries;  //!< Unprocessed checkpoints replaced by certified recombination.
};

/** @brief Local trajectory completion and provenance of its birth. */
struct HybridEvolutionTrajectoryMetadata {
    Bool complete = false;                                           //!< All local pieces processed and all derived children published; independent of descendants.
    std::optional<HybridEvolutionEntryOrigin> origin;                //!< Parent entry/successor producing the seed; absent at the root and at a recombined replacement.
    std::vector<HybridEvolutionCheckpointId> recombined_checkpoints; //!< Exact entry-level provenance of a replacement.
    std::size_t pending_entries = 0u;                                //!< Number of still-pending work items in this trajectory.
    bool whole_trajectory = true;                                    //!< Whole-source urgent continuation eligibility.
    std::optional<double> recombination_time_barrier;                //!< Upper time bound of the last merge; ready-queue re-merging requires strict time progress.
};

/** @brief Generic two-parameter genealogy with Evolution-owned objects and metadata. */
using HybridEvolutionTree = DerivationTree<HybridEvolutionTrajectory, HybridEvolutionTrajectoryMetadata>;

/** @brief Pending collection on which a phased admission policy evaluates merges. */
enum class HybridRecombinationQueue : std::uint8_t { READY,
                                                     WAITING };

/** @brief Boundaries at which pending recombination can be reconsidered. */
enum class HybridRecombinationEvaluationPoint : std::uint8_t {
    AFTER_STEP,      //!< A numerical result and all its successors have been published.
    BEFORE_ADMISSION //!< Pending work is about to become eligible for execution.
};

/** @brief Immutable evaluation timing, independent of admission and grouping.
 *  @details completed_steps always counts committed steps globally in this run.
 */
class HybridRecombinationTriggerPolicy {
  public:
    using NodeId = HybridEvolutionTree::NodeId;
    /** @brief Destroy a policy through its interface. */
    virtual ~HybridRecombinationTriggerPolicy() = default;
    /** @brief Clone immutable evaluation settings. */
    virtual HybridRecombinationTriggerPolicy* clone() const = 0;
    /** @brief Decide whether to inspect the configured pending collection. */
    virtual bool should_evaluate(const HybridEvolutionTree&, NodeId producer,
                                 HybridRecombinationEvaluationPoint point,
                                 std::size_t completed_steps, HybridRecombinationQueue queue) const = 0;
};

/** @brief Inspect pending work just before admission.
 *  @details Immediate admission evaluates ready work after each publication;
 *  parent admission evaluates waiting work before releasing a completed
 *  continuous family; location admission evaluates waiting work at its barrier.
 */
class BeforeAdmissionRecombinationTriggerPolicy final : public HybridRecombinationTriggerPolicy {
  public:
    /** @brief Clone immutable policy settings. */
    BeforeAdmissionRecombinationTriggerPolicy* clone() const override;
    /** @brief Evaluate the configured publication or admission boundary. */
    bool should_evaluate(const HybridEvolutionTree&, NodeId,
                         HybridRecombinationEvaluationPoint, std::size_t, HybridRecombinationQueue) const override;
};

/** @brief Inspect pending work after every committed numerical step. */
class AfterEachStepRecombinationTriggerPolicy final : public HybridRecombinationTriggerPolicy {
  public:
    /** @brief Clone immutable policy settings. */
    AfterEachStepRecombinationTriggerPolicy* clone() const override;
    /** @brief Evaluate the configured publication or admission boundary. */
    bool should_evaluate(const HybridEvolutionTree&, NodeId,
                         HybridRecombinationEvaluationPoint, std::size_t, HybridRecombinationQueue) const override;
};

/** @brief Inspect every N global steps and flush before admission. */
class EveryNStepsRecombinationTriggerPolicy final : public HybridRecombinationTriggerPolicy {
    std::size_t _period;

  public:
    /** @brief Configure a strictly positive evaluation period. */
    explicit EveryNStepsRecombinationTriggerPolicy(std::size_t period);
    /** @brief Clone immutable policy settings. */
    EveryNStepsRecombinationTriggerPolicy* clone() const override;
    /** @brief Evaluate the configured publication or admission boundary. */
    bool should_evaluate(const HybridEvolutionTree&, NodeId,
                         HybridRecombinationEvaluationPoint, std::size_t, HybridRecombinationQueue) const override;
};

/** @brief Optional genealogy restriction on otherwise compatible checkpoints. */
enum class HybridRecombinationScope : std::uint8_t {
    ANY,     //!< Allow checkpoints from unrelated branches of this run.
    SIBLINGS //!< Require seed checkpoints with the same immediate parent and successor meaning.
};

/** @brief Immutable grouping of pending checkpoint identities.
 *  @details Only the supplied pending IDs can be selected, in disjoint groups.
 *  The backend independently validates compatibility, provenance and coverage.
 */
class HybridRecombinationGroupingPolicy {
  public:
    using CheckpointId = HybridEvolutionCheckpointId;
    /** @brief Destroy a policy through its interface. */
    virtual ~HybridRecombinationGroupingPolicy() = default;
    /** @brief Clone immutable grouping settings. */
    virtual HybridRecombinationGroupingPolicy* clone() const = 0;
    /** @brief Propose disjoint groups of at least two pending checkpoints. */
    virtual std::vector<std::vector<CheckpointId>> select(
        const HybridEvolutionTree&, const std::vector<CheckpointId>& pending) const = 0;
    /** @brief Whether this policy may request numerical recombination. */
    virtual bool enables_recombination() const = 0;
    /** @brief Whether the backend must additionally require time overlap with the group anchor. */
    virtual bool requires_time_overlap() const { return true; }
};

/** @brief Keep every original checkpoint without numerical recombination. */
class NoRecombinationPolicy final : public HybridRecombinationGroupingPolicy {
  public:
    /** @brief Clone immutable policy settings. */
    NoRecombinationPolicy* clone() const override;
    /** @brief Select disjoint compatible groups from the supplied pending IDs. */
    std::vector<std::vector<CheckpointId>> select(
        const HybridEvolutionTree&, const std::vector<CheckpointId>&) const override;
    /** @brief Report whether numerical recombination is enabled. */
    bool enables_recombination() const override { return false; }
};

/** @brief Group all compatible pending checkpoints without time or cardinality limits. */
class AllCompatibleRecombinationPolicy final : public HybridRecombinationGroupingPolicy {
    HybridRecombinationScope _scope;

  public:
    /** @brief Optionally restrict grouping to sibling seeds. */
    explicit AllCompatibleRecombinationPolicy(HybridRecombinationScope scope = HybridRecombinationScope::ANY);
    /** @brief Clone immutable policy settings. */
    AllCompatibleRecombinationPolicy* clone() const override;
    /** @brief Select disjoint compatible groups from the supplied pending IDs. */
    std::vector<std::vector<CheckpointId>> select(
        const HybridEvolutionTree&, const std::vector<CheckpointId>&) const override;
    /** @brief Report whether numerical recombination is enabled. */
    bool enables_recombination() const override { return true; }
    /** @brief Report whether grouping additionally requires time overlap. */
    bool requires_time_overlap() const override { return false; }
};

/** @brief Greedy grouping with a maximum group size and time span.
 *  @details ANY can select nonadjacent checkpoints but requires time overlap
 *  with the anchor. SIBLINGS groups consecutive compatible sibling seeds and
 *  permits disjoint times within the configured span, preserving sibling-window selection.
 */
class WindowRecombinationPolicy final : public HybridRecombinationGroupingPolicy {
    std::size_t _max_group_size;
    double _max_time_spread;
    HybridRecombinationScope _scope;

  public:
    /** @brief Configure group size >= 2, finite nonnegative span and genealogy scope. */
    explicit WindowRecombinationPolicy(std::size_t max_group_size = 8u,
                                       double max_time_spread = 0.125,
                                       HybridRecombinationScope scope = HybridRecombinationScope::ANY);
    /** @brief Clone immutable policy settings. */
    WindowRecombinationPolicy* clone() const override;
    /** @brief Select disjoint compatible groups from the supplied pending IDs. */
    std::vector<std::vector<CheckpointId>> select(
        const HybridEvolutionTree&, const std::vector<CheckpointId>&) const override;
    /** @brief Report whether numerical recombination is enabled. */
    bool enables_recombination() const override { return true; }
    /** @brief Report whether grouping additionally requires time overlap. */
    bool requires_time_overlap() const override { return _scope == HybridRecombinationScope::ANY; }
};

/** @brief Immutable policy separating checkpoint admission from genealogy and dispatch.
 *  @details The scheduler closes a batch only when its ready queue and in-flight
 *  count are both empty. select_next() is also called at intermediate
 *  publications, allowing parent admission to release independent completed
 *  families. Empty selections are allowed until the batch completes.
 */
class HybridCheckpointAdmissionPolicy {
  public:
    using CheckpointId = HybridEvolutionCheckpointId;
    /** @brief Destroy an admission strategy through its interface. */
    virtual ~HybridCheckpointAdmissionPolicy() = default;
    /** @brief Copy immutable admission settings without copying scheduler state. */
    virtual HybridCheckpointAdmissionPolicy* clone() const = 0;
    /** @brief Whether this successor joins the ready queue immediately.
     *  @param derived_trajectory True if this successor created a new genealogy node.
     */
    virtual bool admit_successor(const HybridEvolutionWorkItem& producer,
                                 const HybridEvolutionSuccessor& successor, bool derived_trajectory) const = 0;
    /** @brief Select unique eligible waiting IDs; completed batches with waiting work require progress.
     *  @param batch_complete True only when ready and in-flight work are both empty.
     */
    virtual std::vector<CheckpointId> select_next(
        const HybridEvolutionTree&, const std::vector<CheckpointId>& waiting, bool batch_complete) const = 0;
    /** @brief Choose the queue inspected at trigger boundaries. */
    virtual HybridRecombinationQueue recombination_queue() const = 0;
};

/** @brief Finish a continuous family before admitting its discrete descendants.
 *  @details Subdivisions remain runnable and form explicit child trajectories.
 *  Their discrete descendants wait until all pieces of that continuous family
 *  finish. Independent families can release children without a global barrier.
 */
class ParentCheckpointAdmissionPolicy final : public HybridCheckpointAdmissionPolicy {
  public:
    /** @brief Clone immutable policy settings. */
    ParentCheckpointAdmissionPolicy* clone() const override;
    /** @brief Decide whether the successor becomes immediately ready. */
    bool admit_successor(const HybridEvolutionWorkItem&, const HybridEvolutionSuccessor&, bool) const override;
    /** @brief Select eligible waiting checkpoints according to this admission policy. */
    std::vector<CheckpointId> select_next(
        const HybridEvolutionTree&, const std::vector<CheckpointId>&, bool) const override;
    /** @brief Return the pending collection used for grouping. */
    HybridRecombinationQueue recombination_queue() const override { return HybridRecombinationQueue::WAITING; }
};

/** @brief Preserve immediate global frontier admission. */
class ImmediateCheckpointAdmissionPolicy final : public HybridCheckpointAdmissionPolicy {
  public:
    /** @brief Clone immutable policy settings. */
    ImmediateCheckpointAdmissionPolicy* clone() const override;
    /** @brief Decide whether the successor becomes immediately ready. */
    bool admit_successor(const HybridEvolutionWorkItem&, const HybridEvolutionSuccessor&, bool) const override;
    /** @brief Select eligible waiting checkpoints according to this admission policy. */
    std::vector<CheckpointId> select_next(
        const HybridEvolutionTree&, const std::vector<CheckpointId>&, bool) const override;
    /** @brief Return the pending collection used for grouping. */
    HybridRecombinationQueue recombination_queue() const override { return HybridRecombinationQueue::READY; }
};

/** @brief Run one location at a time, retaining all discrete successors until batch completion.
 *  @details Continuous continuations, reconditioning and subdivisions remain
 *  ready in the current location. Discrete successors, including self-loops,
 *  wait for a later batch. The first waiting location is selected next, with
 *  all its checkpoints admitted together, independently of their parents.
 */
class LocationCheckpointAdmissionPolicy final : public HybridCheckpointAdmissionPolicy {
  public:
    /** @brief Clone immutable policy settings. */
    LocationCheckpointAdmissionPolicy* clone() const override;
    /** @brief Decide whether the successor becomes immediately ready. */
    bool admit_successor(const HybridEvolutionWorkItem&, const HybridEvolutionSuccessor&, bool) const override;
    /** @brief Select eligible waiting checkpoints according to this admission policy. */
    std::vector<CheckpointId> select_next(
        const HybridEvolutionTree&, const std::vector<CheckpointId>&, bool) const override;
    /** @brief Return the pending collection used for grouping. */
    HybridRecombinationQueue recombination_queue() const override { return HybridRecombinationQueue::WAITING; }
};

/** @brief Choose whether to delay execution of the first ready checkpoint.
 *  @details Called only when another runnable checkpoint is available. Deferral
 *  preserves work and does not change the physical time of any enclosure.
 *  Decisions are allowed to depend on bounded previous deferral counts.
 */
class HybridCheckpointDispatchPolicy {
  public:
    using CheckpointId = HybridEvolutionCheckpointId;
    /** @brief Destroy a policy through its interface. */
    virtual ~HybridCheckpointDispatchPolicy() = default;
    /** @brief Clone immutable policy settings. */
    virtual HybridCheckpointDispatchPolicy* clone() const = 0;
    /** @brief Whether dispatch needs snapshots of alternative ready checkpoints. */
    virtual bool enables_deferral() const { return true; }
    /** @brief Request that a pending checkpoint be moved to the queue tail. */
    virtual bool should_defer(const HybridEvolutionTree& history, CheckpointId candidate,
                              const std::vector<CheckpointId>& other_ready,
                              std::size_t previous_deferrals) const = 0;
};

/** @brief Preserve original FIFO work dispatch without waiting. */
class ImmediateCheckpointDispatchPolicy final : public HybridCheckpointDispatchPolicy {
  public:
    /** @brief Clone immutable policy settings. */
    ImmediateCheckpointDispatchPolicy* clone() const override;
    /** @brief FIFO never needs alternative-ready snapshots. */
    bool enables_deferral() const override { return false; }
    /** @brief Decide whether to postpone this ready checkpoint. */
    bool should_defer(const HybridEvolutionTree&, CheckpointId,
                      const std::vector<CheckpointId>&, std::size_t) const override;
};

/** @brief Give earlier-time work a bounded chance to produce merge partners.
 *  @details A ready checkpoint is delayed only while a different trajectory
 *  has pending work at an earlier physical time within max_time_lag. Every
 *  checkpoint is deferred at most max_deferrals times; progress is guaranteed
 *  even when no compatible partner ever arrives.
 */
class BoundedTimeLagCheckpointDispatchPolicy final : public HybridCheckpointDispatchPolicy {
    std::size_t _max_deferrals;
    double _max_time_lag;

  public:
    explicit BoundedTimeLagCheckpointDispatchPolicy(std::size_t max_deferrals = 8u,
                                                    double max_time_lag = 1.0);
    /** @brief Clone immutable policy settings. */
    BoundedTimeLagCheckpointDispatchPolicy* clone() const override;
    /** @brief Decide whether to postpone this ready checkpoint. */
    bool should_defer(const HybridEvolutionTree& history, CheckpointId candidate,
                      const std::vector<CheckpointId>& other_ready,
                      std::size_t previous_deferrals) const override;
};

/** @brief Optional quality/cost gate after a sound numerical recombination.
 *  @details Rejection keeps all original pending sources intact. This decision
 *  does not replace the recombiner's sound covering contract.
 */
class HybridRecombinationAcceptancePolicy {
  public:
    using CheckpointId = HybridEvolutionCheckpointId;
    /** @brief Destroy a policy through its interface. */
    virtual ~HybridRecombinationAcceptancePolicy() = default;
    /** @brief Clone immutable policy settings. */
    virtual HybridRecombinationAcceptancePolicy* clone() const = 0;
    /** @brief Whether this gate can assess sources from the selected pending collection. */
    virtual bool supports_queue(HybridRecombinationQueue) const { return true; }
    /** @brief Assess the quality of a numerically valid recombination. */
    virtual bool accept(const HybridEvolutionTree& history,
                        const std::vector<CheckpointId>& sources,
                        const EnclosureRecombinationResult& result) const = 0;
};

/** @brief Accept every numerically valid reduction, preserving existing behavior. */
class AcceptAllRecombinationsPolicy final : public HybridRecombinationAcceptancePolicy {
  public:
    /** @brief Clone immutable policy settings. */
    AcceptAllRecombinationsPolicy* clone() const override;
    /** @brief Assess the quality of a numerically valid recombination. */
    bool accept(const HybridEvolutionTree&, const std::vector<CheckpointId>&,
                const EnclosureRecombinationResult&) const override;
};

/** @brief Require controlled output-radius amplification and nonexpansive history.
 *  @details Uses observed consecutive continuous checkpoint radii as a local
 *  heuristic, not a formal contraction certificate for the vector field.
 *  Rejects sources without a previously executed continuous continuation
 *  step, to avoid treating a discrete reset or subdivision as contraction.
 */
class RadiusTrendRecombinationAcceptancePolicy final : public HybridRecombinationAcceptancePolicy {
    double _max_output_amplification;
    double _max_continuation_growth;

  public:
    /** @brief Both ratios must be finite and strictly positive. */
    explicit RadiusTrendRecombinationAcceptancePolicy(double max_output_amplification = 2.0,
                                                      double max_continuation_growth = 1.0);
    /** @brief Clone immutable policy settings. */
    RadiusTrendRecombinationAcceptancePolicy* clone() const override;
    /** @brief Report whether this quality gate supports the pending collection. */
    bool supports_queue(HybridRecombinationQueue queue) const override { return queue == HybridRecombinationQueue::READY; }
    /** @brief Assess the quality of a numerically valid recombination. */
    bool accept(const HybridEvolutionTree& history,
                const std::vector<CheckpointId>& sources,
                const EnclosureRecombinationResult& result) const override;
};

/** @brief Trajectory-tree backend with one checkpoint scheduler for every admission policy.
 *  @details Only unprocessed pending checkpoints can be recombined. Numerical
 *  evolution runs outside the publication lock; admission, history and merges
 *  are committed together. Recombination requires upper semantics.
 */
class TreeHybridEvolutionBackend final : public HybridEvolutionBackendInterface {
    class Impl;
    std::unique_ptr<Impl> _impl;

  public:
    /** @brief Configure admission, grouping, evaluation, numerical construction and optional gates.
     *  @details Defaults preserve parent admission without recombination.
     *  A merge gate unsupported by the selected queue is rejected when recombination is enabled.
     */
    explicit TreeHybridEvolutionBackend(
        const HybridEnclosure& initial,
        const HybridCheckpointAdmissionPolicy& admission = ParentCheckpointAdmissionPolicy{},
        const HybridRecombinationGroupingPolicy& grouping = NoRecombinationPolicy{},
        const HybridRecombinationTriggerPolicy& trigger = BeforeAdmissionRecombinationTriggerPolicy{},
        const EnclosureRecombinerInterface& recombiner = ReferenceEnclosureRecombiner{},
        const EnclosureRecombinationContext& context = {},
        const HybridCheckpointDispatchPolicy& dispatch = ImmediateCheckpointDispatchPolicy{},
        const HybridRecombinationAcceptancePolicy& acceptance = AcceptAllRecombinationsPolicy{});
    /** @brief Release the owned trajectory history. */
    ~TreeHybridEvolutionBackend() override;
    /** @brief Disallow copying the owner and its trajectory identities. */
    TreeHybridEvolutionBackend(const TreeHybridEvolutionBackend&) = delete;
    /** @brief Disallow copying the owner and its trajectory identities. */
    TreeHybridEvolutionBackend& operator=(const TreeHybridEvolutionBackend&) = delete;
    /** @brief Supply the approximation semantics before running. */
    Void set_semantics(Semantics semantics) override;
    /** @brief Supply a finite-time progress horizon before running; omission leaves the estimate unavailable. */
    Void set_time_horizon(const Real& final_time);
    /** @brief Read set-time progress consistently with checkpoint publication and recombination. */
    HybridEvolutionProgress progress() const override;
    /** @brief Run once with the selected admission and recombination policies.
     *  @note Callback failures leave an inspectable incomplete history and prohibit Orbit release.
     */
    Void run(const StepFunction& advance) override;
    /** @brief Materialize the conventional Orbit once, retaining history for const inspection. */
    Orbit<HybridEnclosure> release_orbit() override;
    /** @brief Read genealogy, checkpoints, step results and local completion.
     *  @note Access is sequential. A failed run can contain processed and unprocessed entries.
     *  The reference is valid until this backend is destroyed; growth can invalidate entry references.
     */
    const HybridEvolutionTree& tree() const noexcept;
};

/** @brief Factory selecting policy-driven trajectory-tree execution for ordinary orbit() calls.
 *  @details The independent grouping, timing, numerical recombiner and radius
 *  context settings are cloned or copied for each run.
 */
class TreeHybridEvolutionBackendFactory : public HybridEvolutionBackendFactoryInterface {
    std::shared_ptr<const HybridCheckpointAdmissionPolicy> _admission;
    std::shared_ptr<const HybridRecombinationGroupingPolicy> _grouping;
    std::shared_ptr<const HybridCheckpointDispatchPolicy> _dispatch;
    std::shared_ptr<const HybridRecombinationAcceptancePolicy> _acceptance;
    std::shared_ptr<const HybridRecombinationTriggerPolicy> _trigger;
    std::shared_ptr<const EnclosureRecombinerInterface> _recombiner;
    EnclosureRecombinationContext _context;

  public:
    /** @brief Clone independent admission, grouping and numerical settings for each run. */
    explicit TreeHybridEvolutionBackendFactory(
        const HybridCheckpointAdmissionPolicy& admission = ParentCheckpointAdmissionPolicy{},
        const HybridRecombinationGroupingPolicy& grouping = NoRecombinationPolicy{},
        const HybridRecombinationTriggerPolicy& trigger = BeforeAdmissionRecombinationTriggerPolicy{},
        const EnclosureRecombinerInterface& recombiner = ReferenceEnclosureRecombiner{},
        const EnclosureRecombinationContext& context = {},
        const HybridCheckpointDispatchPolicy& dispatch = ImmediateCheckpointDispatchPolicy{},
        const HybridRecombinationAcceptancePolicy& acceptance = AcceptAllRecombinationsPolicy{});
    /** @brief Clone the immutable factory settings without copying any run history. */
    TreeHybridEvolutionBackendFactory* clone() const override;
    /** @brief Create independent trajectory history; the numerical callback enforces the horizon. */
    std::unique_ptr<HybridEvolutionBackendInterface> create(const HybridEnclosure& initial, const Real& final_time) const override;
};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLUTION_TREE_HPP
