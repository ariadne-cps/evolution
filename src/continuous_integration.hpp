/***************************************************************************
 *            continuous_integration.hpp
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

/** @file src/continuous_integration.hpp
 *  @brief Shared validated integration step preparation and application.
 */

#ifndef ARIADNE_EVOLUTION_CONTINUOUS_INTEGRATION_HPP
#define ARIADNE_EVOLUTION_CONTINUOUS_INTEGRATION_HPP

#include "evolution/enclosure.hpp"
#include "evolution/evolver_interface.hpp"
#include "integration/validated_integration_step.hpp"

namespace Ariadne {
class ValidatedIntegrationStep;
class OrdinaryDifferentialAlgebraicInclusion;
class SelectedParametricMultifunctionIntegrator;
class FirstOrderMultifunctionIntegrator;

namespace EvolutionDetail {
/** @brief Match the ordinary/hybrid equation field's explicit state/time/parameter layout.
 * Autonomous parameterized fields acquire an ignored numerical clock argument.
 */
EffectiveVectorMultivariateFunction prepare_equation_field(const EffectiveVectorMultivariateFunction&,
                                                           SizeType state_dimension, bool timed, SizeType parameter_count);
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

/** @brief A numerical step conditioned on an enclosure's external clock map.
 *  The first numerical conditioning argument is the clock origin for this step;
 *  it is composed from the source, never adjoined as a physical state or a
 *  trajectory-constant model parameter. Local duration advances that clock.
 */
struct AnchoredIntegrationStep {
    ValidatedIntegrationStep step;
    ValidatedVectorMultivariateFunctionPatch arguments;
};
AnchoredIntegrationStep compute_clocked_equation_step(
    const EffectiveVectorMultivariateFunction& field, const Enclosure&,
    const ValidatedVectorMultivariateFunctionPatch& model_parameters, const Dyadic& maximum_step, const IntegratorInterface&);
AnchoredIntegrationStep compute_clocked_inclusion_step(
    const OrdinaryDifferentialAlgebraicInclusion&, const Enclosure&,
    const ValidatedVectorMultivariateFunctionPatch& model_parameters, const Dyadic& maximum_step, const SelectedParametricMultifunctionIntegrator&);
AnchoredIntegrationStep compute_clocked_inclusion_step(
    const OrdinaryDifferentialAlgebraicInclusion&, const Enclosure&,
    const ValidatedVectorMultivariateFunctionPatch& model_parameters, const Dyadic& maximum_step, const FirstOrderMultifunctionIntegrator&);

/** @brief Whether the maps take physical states or existing enclosure parameters. */
enum class IntegrationStepCoordinates { STATE,
                                        PARAMETERS };

/** @brief Whether uncorrelated outer boxes may replace unavailable certificates. */
enum class IntegrationStepRepresentation { REQUIRE_STRUCTURED,
                                           ALLOW_BOX_UPPER_ENCLOSURES,
                                           ALLOW_BOX_LOWER_ENCLOSURES };

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
 *  Unavailable maps require an explicit permission matching the requested
 *  semantics. LOWER permission means that the caller has proved the certified
 *  boxes enclose the continuation of a known nonempty reached subset; replacing
 *  its correlations by independent ERROR axes only enlarges that enclosure.
 */
AppliedIntegrationStep apply_integration_step(
    const Enclosure& initial, const ValidatedIntegrationStep& step,
    IntegrationStepCoordinates coordinates, Semantics semantics,
    IntegrationStepRepresentation representation = IntegrationStepRepresentation::REQUIRE_STRUCTURED);

/** @brief Apply a step using a map from existing enclosure parameters to
 *  the numerical [state, constant parameters] arguments, excluding time.
 *  Model parameters are composed from their original coordinates on every step.
 *  A clocked step additionally composes its explicit source-clock conditioning
 *  coordinate; this is not persisted as a model parameter.
 */
AppliedIntegrationStep apply_integration_step(
    const Enclosure& initial, const ValidatedIntegrationStep& step,
    const ValidatedVectorMultivariateFunctionPatch& arguments, Semantics semantics,
    IntegrationStepRepresentation representation = IntegrationStepRepresentation::REQUIRE_STRUCTURED);
} // namespace EvolutionDetail
} // namespace Ariadne
#endif
