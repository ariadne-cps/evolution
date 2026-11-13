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

/** @brief One trajectory's ordered history; a trajectory is one generic tree node. */
struct HybridEvolutionTrajectory {
    std::vector<HybridEvolutionEntry> entries; //!< Seed followed by continuations and numerical pieces.
};

/** @brief Local trajectory completion and provenance of its birth. */
struct HybridEvolutionTrajectoryMetadata {
    Bool complete = false;                            //!< All local pieces processed and all derived children published; independent of descendants.
    std::optional<HybridEvolutionEntryOrigin> origin; //!< Parent entry/successor producing the seed; absent at the root.
};

/** @brief Generic two-parameter genealogy with Evolution-owned objects and metadata. */
using HybridEvolutionTree = DerivationTree<HybridEvolutionTrajectory, HybridEvolutionTrajectoryMetadata>;

/** @brief Serial backend completing a trajectory before scheduling its derived children.
 *  @details Continuous steps, reconditioning and subdivision extend a node's
 *  history. Unresolved/permissive resets create child trajectories. A sole urgent
 *  successor with a whole-source certificate continues an undivided node. Once
 *  subdivision or an urgent derived branch separates portions, subsequent resets
 *  derive children even if the final remaining piece has a coverage certificate. Orbit is
 *  projected from explicit history contributions after successful completion.
 *  No recombination, history editing or concurrent tree mutation is performed.
 */
class TreeHybridEvolutionBackend final : public HybridEvolutionBackendInterface {
    class Impl;
    std::unique_ptr<Impl> _impl;

  public:
    /** @brief Create the root trajectory with its unprocessed initial checkpoint. */
    explicit TreeHybridEvolutionBackend(const HybridEnclosure& initial);
    /** @brief Release the owned trajectory history. */
    ~TreeHybridEvolutionBackend() override;
    /** @brief Disallow copying the owner and its trajectory identities. */
    TreeHybridEvolutionBackend(const TreeHybridEvolutionBackend&) = delete;
    /** @brief Disallow copying the owner and its trajectory identities. */
    TreeHybridEvolutionBackend& operator=(const TreeHybridEvolutionBackend&) = delete;
    /** @brief Run once, processing each trajectory's pieces before releasing its children.
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

/** @brief Factory selecting the serial trajectory-tree backend for ordinary orbit() calls. */
class TreeHybridEvolutionBackendFactory : public HybridEvolutionBackendFactoryInterface {
  public:
    /** @brief Clone the stateless factory without copying any run history. */
    TreeHybridEvolutionBackendFactory* clone() const override;
    /** @brief Create independent trajectory history; the numerical callback enforces the horizon. */
    std::unique_ptr<HybridEvolutionBackendInterface> create(const HybridEnclosure& initial, const Real& final_time) const override;
};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLUTION_TREE_HPP
