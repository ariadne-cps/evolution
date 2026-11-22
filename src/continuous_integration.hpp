/***************************************************************************
 *            continuous_integration.hpp
 *  Copyright  2026  Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#ifndef ARIADNE_EVOLUTION_CONTINUOUS_INTEGRATION_HPP
#define ARIADNE_EVOLUTION_CONTINUOUS_INTEGRATION_HPP

#include "evolution/enclosure.hpp"
#include "evolution/evolver_interface.hpp"

namespace Ariadne {
class ValidatedIntegrationStep;
class OrdinaryDifferentialAlgebraicInclusion;
class SelectedParametricMultifunctionIntegrator;
class FirstOrderMultifunctionIntegrator;

namespace EvolutionDetail {
/** @brief Shared inclusion request/bound construction; no workload or Orbit is created. */
ValidatedIntegrationStep compute_inclusion_step(
    const OrdinaryDifferentialAlgebraicInclusion&, const Enclosure&,
    const ValidatedVectorMultivariateFunctionPatch& constant_parameters,
    const Dyadic& time, const Dyadic& maximum_step,
    const SelectedParametricMultifunctionIntegrator& solver);
ValidatedIntegrationStep compute_inclusion_step(
    const OrdinaryDifferentialAlgebraicInclusion&, const Enclosure&,
    const ValidatedVectorMultivariateFunctionPatch& constant_parameters,
    const Dyadic& time, const Dyadic& maximum_step,
    const FirstOrderMultifunctionIntegrator& solver);

/** @brief Whether the maps take physical states or existing enclosure parameters. */
enum class IntegrationStepCoordinates { STATE,
                                        PARAMETERS };

/** @brief Whether uncorrelated outer boxes may replace unavailable certificates. */
enum class IntegrationStepRepresentation { REQUIRE_STRUCTURED,
                                           ALLOW_BOX_ENCLOSURES };

/** @brief Each reach patch retains its time interval and uncertainty domain.
 *  The endpoint is the independent evolve certificate. A box component is an
 *  explicit outer Cartesian enclosure, never a correlated solver model.
 */
struct AppliedIntegrationStep {
    List<Enclosure> reach;
    Enclosure evolved;
    bool used_box_enclosure = false;
};

/** @brief Apply a validated step without mutating the source enclosure.
 *  Preserves available structure, constraints, outputs, time and fixed parameters.
 *  New input coordinates are appended independently for reach and endpoint.
 *  Unavailable maps require explicit ALLOW_BOX_ENCLOSURES and UPPER semantics.
 *  Structured LOWER application still requires a solver guarantee from the caller.
 */
AppliedIntegrationStep apply_integration_step(
    const Enclosure& initial, const ValidatedIntegrationStep& step,
    IntegrationStepCoordinates coordinates, Semantics semantics,
    IntegrationStepRepresentation representation = IntegrationStepRepresentation::REQUIRE_STRUCTURED);

/** @brief Apply a step using a map from existing enclosure parameters to
 *  the numerical [state, constant parameters] arguments, excluding time.
 *  Fixed parameters are composed from their original coordinates on every step.
 */
AppliedIntegrationStep apply_integration_step(
    const Enclosure& initial, const ValidatedIntegrationStep& step,
    const ValidatedVectorMultivariateFunctionPatch& arguments, Semantics semantics,
    IntegrationStepRepresentation representation = IntegrationStepRepresentation::REQUIRE_STRUCTURED);
} // namespace EvolutionDetail
} // namespace Ariadne
#endif
