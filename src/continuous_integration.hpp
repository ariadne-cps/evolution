/***************************************************************************
 *            continuous_integration.hpp
 *  Copyright  2026  Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#ifndef ARIADNE_EVOLUTION_CONTINUOUS_INTEGRATION_HPP
#define ARIADNE_EVOLUTION_CONTINUOUS_INTEGRATION_HPP

namespace Ariadne {
class Enclosure;
class ValidatedIntegrationStep;

namespace EvolutionDetail {
/** @brief Whether the maps take physical states or existing enclosure parameters. */
enum class IntegrationStepCoordinates { STATE,
                                        PARAMETERS };

/** @brief Apply a common autonomous step to separate reach and endpoint copies.
 *  Preserves constraints, auxiliary outputs, time and parameter dependence.
 *  This bridge currently requires one local-time reach patch and an endpoint
 *  map. Unsupported representations fail before either enclosure is changed.
 */
void apply_integration_step(Enclosure& reach, Enclosure& evolved,
                            const ValidatedIntegrationStep& step, IntegrationStepCoordinates coordinates);
} // namespace EvolutionDetail
} // namespace Ariadne
#endif
