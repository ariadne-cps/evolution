/***************************************************************************
 *            simulation_parameters.cpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#include "simulation_parameters.hpp"

#include <cmath>
#include <stdexcept>

namespace Ariadne::EvolutionDetail {
Void validate_simulation_parameters(const RealSpace& space, const ExactBoxType& domain,
                                    const RealValuation& values, bool exact_names) {
    ARIADNE_PRECONDITION(space.dimension() == domain.dimension());
    if (exact_names && values.values().size() != space.dimension())
        throw std::invalid_argument("Approximate evolution requires exactly the declared constant parameter names");
    for (SizeType i = 0u; i != space.dimension(); ++i) {
        const auto variable = space[i];
        if (!values.values().has_key(variable.name()))
            throw std::invalid_argument("Approximate evolution requires an explicit value for every model parameter");
        const FloatDPBounds bounds(values.get(variable), dp);
        if (!std::isfinite(bounds.lower_raw().get_d()) || !std::isfinite(bounds.upper_raw().get_d()))
            throw std::invalid_argument("Approximate evolution requires finite constant parameter values");
        const ExactIntervalType enclosure(bounds.lower_raw(), bounds.upper_raw());
        if (!definitely(subset(enclosure, domain[i])))
            throw std::invalid_argument("Cannot certify the selected constant parameter value within its declared domain");
    }
}
} // namespace Ariadne::EvolutionDetail
