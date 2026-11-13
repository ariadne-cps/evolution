/***************************************************************************
 *            evolution/hybrid_evolution_backend.hpp
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

/** @file evolution/hybrid_evolution_backend.hpp
 *  @brief Step results and replaceable execution backends for hybrid evolution.
 */

#ifndef ARIADNE_HYBRID_EVOLUTION_BACKEND_HPP
#define ARIADNE_HYBRID_EVOLUTION_BACKEND_HPP

#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <vector>

#include "evolution/hybrid_orbit.hpp"
#include "system/discrete_event.hpp"
#include "system/hybrid_automaton_interface.hpp"

namespace Ariadne {

enum class CrossingKind : std::uint8_t;

/** @brief Processing required by an enclosure, independent of its provenance. */
enum class HybridEvolutionPhase : std::uint8_t {
    PROCESS_STARTING_EVENTS, //!< Check termination and initially active events, including at the root.
    EVOLVE_CONTINUOUSLY      //!< Starting events have already been processed.
};

/** @brief One computational checkpoint; it is not an intermediate Orbit contribution. */
struct HybridEvolutionWorkItem {
    HybridEnclosure enclosure;  //!< Enclosure to process, retaining its numerical enclosure configuration.
    HybridEvolutionPhase phase; //!< Entry-processing phase.
};

/** @brief Meaning of a successor, independent of the backend's storage topology. */
enum class HybridEvolutionSuccessorKind : std::uint8_t {
    CONTINUATION,  //!< Continue local processing without deriving a trajectory.
    TRANSITION,    //!< A reset image; resolution and provenance require event information.
    SUBDIVISION,   //!< A numerical piece, not a sequential time step.
    RECONDITIONING //!< A new numerical representation without time progress.
};

/** @brief Numerical evidence about the source covered by an urgent transition.
 *  @note A transverse crossing or a sole successor is not by itself a coverage proof.
 */
enum class HybridEvolutionTransitionResolution : std::uint8_t {
    UNRESOLVED, //!< Coverage of the whole source is not established, or the transition is permissive.
    FULL_SOURCE //!< The urgent transition is certified to apply to the whole invariant-valid source.
};

/** @brief A successor and the information needed to preserve its origin. */
struct HybridEvolutionSuccessor {
    HybridEvolutionWorkItem work;                                                                     //!< Checkpoint to schedule.
    HybridEvolutionSuccessorKind kind;                                                                //!< Reason for producing the checkpoint.
    std::optional<DiscreteEvent> event;                                                               //!< Transition event, when applicable.
    std::optional<EventKind> event_kind;                                                              //!< Urgency/permissiveness, when applicable.
    std::optional<CrossingKind> crossing_kind;                                                        //!< Crossing classification; absent for starting events.
    HybridEvolutionTransitionResolution resolution = HybridEvolutionTransitionResolution::UNRESOLVED; //!< Source coverage evidence.
};

/** @brief Conventional Orbit collection to which an enclosure contributes. */
enum class HybridEvolutionSetKind : std::uint8_t {
    REACH,        //!< A swept set, including explicitly recorded pre-reset sets.
    INTERMEDIATE, //!< An intermediate output, not necessarily the next work item.
    FINAL         //!< An explicitly selected final set.
};

/** @brief One explicitly classified contribution, separate from future work. */
struct HybridEvolutionContribution {
    HybridEnclosure enclosure;   //!< Validated output enclosure.
    HybridEvolutionSetKind kind; //!< Collection in the returned Orbit.
};

/** @brief Reason why this work item has no ordinary continuation; not an Orbit role. */
enum class HybridEvolutionStopReason : std::uint8_t {
    NONE,
    NO_CONTINUATION,
    EMPTY,
    PAST_HORIZON,
    TIME_LIMIT,
    DISCRETE_LIMIT,
    TERMINATING_EVENT,
    LOWER_SIZE_LIMIT
};

/** @brief Complete local result, published by a backend after numerical processing succeeds.
 *  @details Final contributions and successors may coexist. Exhausting work does
 *  not implicitly add a final set or complete all work of a trajectory.
 */
struct HybridEvolutionStepResult {
    std::vector<HybridEvolutionContribution> contributions;                  //!< Output contributions in production order.
    std::vector<HybridEvolutionSuccessor> successors;                        //!< Future work in production order.
    HybridEvolutionStopReason stop_reason = HybridEvolutionStopReason::NONE; //!< Local stopping information.

    /** @brief Record a swept reachable set. */
    Void record_reach(const HybridEnclosure& set);
    /** @brief Record an intermediate output without scheduling it implicitly. */
    Void record_intermediate(const HybridEnclosure& set);
    /** @brief Record a terminal output without assuming there are no successors. */
    Void record_final(const HybridEnclosure& set);
    /** @brief Record a continuation, subdivision or reconditioned checkpoint. */
    Void continue_with(const HybridEnclosure& set, HybridEvolutionSuccessorKind kind = HybridEvolutionSuccessorKind::CONTINUATION);
    /** @brief Record a post-reset checkpoint with explicit event provenance. */
    Void transition_to(const HybridEnclosure& set, DiscreteEvent event, EventKind event_kind,
                       std::optional<CrossingKind> crossing_kind = std::nullopt,
                       HybridEvolutionTransitionResolution resolution = HybridEvolutionTransitionResolution::UNRESOLVED);
};

/** @brief Consistent result counters for logging; not a numerical stopping criterion. */
struct HybridEvolutionStatistics {
    SizeType reach = 0;        //!< Published reach contributions.
    SizeType intermediate = 0; //!< Published intermediate contributions.
    SizeType final = 0;        //!< Published final contributions.
};

/** @brief One run's execution/history owner, independent of the numerical step implementation.
 *  @details A backend schedules checkpoints and consumes complete step results.
 *  Implementations own synchronization and join workers before run() returns.
 */
class HybridEvolutionBackendInterface {
  public:
    using StepFunction = std::function<HybridEvolutionStepResult(const HybridEvolutionWorkItem&, const HybridEvolutionStatistics&)>;
    /** @brief Destroy the backend after all its worker callbacks have joined. */
    virtual ~HybridEvolutionBackendInterface() = default;
    /** @brief Execute the run once, propagating errors after outstanding workers join. */
    virtual Void run(const StepFunction& advance) = 0;
    /** @brief Transfer the conventional Orbit after successful run completion, once. */
    virtual Orbit<HybridEnclosure> release_orbit() = 0;
};

/** @brief Immutable selection policy creating independent state for each orbit() call. */
class HybridEvolutionBackendFactoryInterface {
  public:
    /** @brief Destroy the factory. */
    virtual ~HybridEvolutionBackendFactoryInterface() = default;
    /** @brief Copy factory settings without copying an active run. */
    virtual HybridEvolutionBackendFactoryInterface* clone() const = 0;
    /** @brief Create one run, preserving the original input as Orbit::initial(). */
    virtual std::unique_ptr<HybridEvolutionBackendInterface> create(const HybridEnclosure& initial, const Real& final_time) const = 0;
};

/** @brief Factory for the existing flat workload and conventional Orbit accumulation. */
class FlatHybridEvolutionBackendFactory : public HybridEvolutionBackendFactoryInterface {
  public:
    /** @brief Copy the stateless flat-backend factory. */
    FlatHybridEvolutionBackendFactory* clone() const override;
    /** @brief Create fresh workload/result state; the factory itself retains no Orbit. */
    std::unique_ptr<HybridEvolutionBackendInterface> create(const HybridEnclosure& initial, const Real& final_time) const override;
};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLUTION_BACKEND_HPP
