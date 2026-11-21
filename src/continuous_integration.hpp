/***************************************************************************
 *            continuous_integration.hpp
 *  Copyright  2026  Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#ifndef ARIADNE_EVOLUTION_CONTINUOUS_INTEGRATION_HPP
#define ARIADNE_EVOLUTION_CONTINUOUS_INTEGRATION_HPP

#include "function/function_patch.hpp"

namespace Ariadne {
class Enclosure;
class ValidatedIntegrationStep;

namespace EvolutionDetail {
/** @brief Whether the maps take physical states or existing enclosure parameters. */
enum class IntegrationStepCoordinates { STATE,
                                        PARAMETERS };

/** @brief Apply a structured validated step to separate reach and endpoint copies.
 *  Preserves constraints, auxiliary outputs, time and parameter dependence.
 *  This bridge requires one reach patch and an endpoint map. An absolute-time
 *  reach argument is rebased to local time without discarding fixed-parameter
 *  dependencies. Unsupported representations fail without mutating inputs.
 */
void apply_integration_step(Enclosure& reach, Enclosure& evolved,
                            const ValidatedIntegrationStep& step, IntegrationStepCoordinates coordinates);

/** @brief Apply a step using a map from existing enclosure parameters to
 *  the numerical [state, constant parameters] arguments, excluding time.
 *  Fixed parameters are composed from their original coordinates on every step.
 */
void apply_integration_step(Enclosure& reach, Enclosure& evolved,
                            const ValidatedIntegrationStep& step,
                            const ValidatedVectorMultivariateFunctionPatch& arguments);
} // namespace EvolutionDetail
} // namespace Ariadne
#endif
