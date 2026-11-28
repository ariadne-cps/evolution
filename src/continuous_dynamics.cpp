/***************************************************************************
 *            continuous_dynamics.cpp
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

/** @file src/continuous_dynamics.cpp
 *  @brief Implement shared continuous-system metadata.
 */

#include "evolution/continuous_dynamics.hpp"

namespace Ariadne {
ContinuousDynamics::ContinuousDynamics(const OrdinaryDifferentialAlgebraicEquation& system)
    : _declaration(SharedPointer<const OrdinaryDifferentialAlgebraicEquation>(system.clone())) {}
ContinuousDynamics::ContinuousDynamics(const OrdinaryDifferentialAlgebraicInclusion& system)
    : _declaration(SharedPointer<const OrdinaryDifferentialAlgebraicInclusion>(system.clone())) {}
ContinuousDynamics* ContinuousDynamics::clone() const { return new ContinuousDynamics(*this); }
SizeType ContinuousDynamics::dimension() const {
    return std::visit([](const auto& system) { return system->dimension(); }, _declaration);
}
RealSpace ContinuousDynamics::state_space() const {
    return std::visit([](const auto& system) -> RealSpace { return system->state_space(); }, _declaration);
}
RealSpace ContinuousDynamics::auxiliary_space() const {
    const auto* value = equation();
    return value != nullptr ? value->auxiliary_space() : RealSpace{};
}
RealSpace ContinuousDynamics::state_auxiliary_space() const {
    const auto* value = equation();
    return value != nullptr ? value->state_auxiliary_space() : state_space();
}
const RealSpace& ContinuousDynamics::parameter_space() const {
    return std::visit([](const auto& system) -> const RealSpace& { return system->parameter_space(); }, _declaration);
}
const BoxDomainType& ContinuousDynamics::constant_parameter_domain() const {
    return std::visit([](const auto& system) -> const BoxDomainType& { return system->constant_parameter_domain(); }, _declaration);
}
bool ContinuousDynamics::has_explicit_time() const {
    return std::visit([](const auto& system) { return system->has_explicit_time(); }, _declaration);
}
const OrdinaryDifferentialAlgebraicEquation* ContinuousDynamics::equation() const {
    const auto* value = std::get_if<SharedPointer<const OrdinaryDifferentialAlgebraicEquation>>(&_declaration);
    return value != nullptr ? value->get() : nullptr;
}
const OrdinaryDifferentialAlgebraicInclusion* ContinuousDynamics::inclusion() const {
    const auto* value = std::get_if<SharedPointer<const OrdinaryDifferentialAlgebraicInclusion>>(&_declaration);
    return value != nullptr ? value->get() : nullptr;
}

} // namespace Ariadne
