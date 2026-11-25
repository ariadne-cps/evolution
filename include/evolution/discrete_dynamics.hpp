/***************************************************************************
 *            evolution/discrete_dynamics.hpp
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

/** @file evolution/discrete_dynamics.hpp
 *  @brief Owning native discrete equation/inclusion metadata for approximate evolution.
 */

#ifndef ARIADNE_DISCRETE_DYNAMICS_HPP
#define ARIADNE_DISCRETE_DYNAMICS_HPP

#include <memory>
#include <variant>

#include "evolution/evolver.decl.hpp"
#include "system/difference_algebraic_equation.hpp"
#include "system/difference_algebraic_inclusion.hpp"

namespace Ariadne {

/** @brief An owning view of native autonomous difference equations and bounded-input inclusions.
 *  @details The original System declarations are retained without modifying
 *  their semantics. Input sampling belongs to the approximate evolver;
 *  this wrapper never turns an inclusion into a deterministic equation.
 */
class DiscreteDynamics {
  public:
    using TimeType = Integer;
    explicit DiscreteDynamics(const DifferenceAlgebraicEquation&);
    explicit DiscreteDynamics(const DifferenceAlgebraicInclusion&);
    DiscreteDynamics* clone() const;
    SizeType dimension() const;
    RealSpace state_space() const;
    RealSpace auxiliary_space() const;
    RealSpace state_auxiliary_space() const;
    Bool has_equation() const;
    Bool has_inclusion() const;

  private:
    std::variant<std::shared_ptr<const DifferenceAlgebraicEquation>,
                 std::shared_ptr<const DifferenceAlgebraicInclusion>> _declaration;

    const DifferenceAlgebraicEquation* equation() const;
    const DifferenceAlgebraicInclusion* inclusion() const;
    friend class DiscreteTimeEvolver<ApproximateTag>;
};

} // namespace Ariadne

#endif // ARIADNE_DISCRETE_DYNAMICS_HPP
