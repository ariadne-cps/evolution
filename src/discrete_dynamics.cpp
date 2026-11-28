/***************************************************************************
 *            src/discrete_dynamics.cpp
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

/** @file src/discrete_dynamics.cpp
 *  @brief Implement owning native discrete equation/inclusion metadata.
 */

#include "evolution/discrete_dynamics.hpp"

namespace Ariadne {

DiscreteDynamics::DiscreteDynamics(const DifferenceAlgebraicEquation& equation)
    : _declaration(std::shared_ptr<const DifferenceAlgebraicEquation>(equation.clone())) {}

DiscreteDynamics::DiscreteDynamics(const DifferenceAlgebraicInclusion& inclusion)
    : _declaration(std::shared_ptr<const DifferenceAlgebraicInclusion>(inclusion.clone())) {}

DiscreteDynamics* DiscreteDynamics::clone() const { return new DiscreteDynamics(*this); }

SizeType DiscreteDynamics::dimension() const {
    return std::visit([](const auto& system) -> SizeType {
        return system->dimension();
    },
                      _declaration);
}

RealSpace DiscreteDynamics::state_space() const {
    return std::visit([](const auto& system) -> RealSpace {
        return system->state_space();
    },
                      _declaration);
}

RealSpace DiscreteDynamics::auxiliary_space() const {
    return std::visit([](const auto& system) -> RealSpace {
        return system->auxiliary_space();
    },
                      _declaration);
}

RealSpace DiscreteDynamics::state_auxiliary_space() const {
    return std::visit([](const auto& system) -> RealSpace {
        return system->state_auxiliary_space();
    },
                      _declaration);
}

Bool DiscreteDynamics::has_equation() const { return equation() != nullptr; }
Bool DiscreteDynamics::has_inclusion() const { return inclusion() != nullptr; }

const DifferenceAlgebraicEquation* DiscreteDynamics::equation() const {
    const auto* handle = std::get_if<std::shared_ptr<const DifferenceAlgebraicEquation>>(&_declaration);
    return handle ? handle->get() : nullptr;
}

const DifferenceAlgebraicInclusion* DiscreteDynamics::inclusion() const {
    const auto* handle = std::get_if<std::shared_ptr<const DifferenceAlgebraicInclusion>>(&_declaration);
    return handle ? handle->get() : nullptr;
}

} // namespace Ariadne
