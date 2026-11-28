/***************************************************************************
 *            evolution/continuous_dynamics.hpp
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

/** @file evolution/continuous_dynamics.hpp
 *  @brief Owning continuous-system metadata shared by evolution paradigms.
 */

#ifndef ARIADNE_CONTINUOUS_DYNAMICS_HPP
#define ARIADNE_CONTINUOUS_DYNAMICS_HPP
#include <variant>

#include "evolution/evolver.decl.hpp"
#include "system/odae.hpp"
#include "system/odai.hpp"
namespace Ariadne {
/** @brief Owning metadata view of an ODAE or ODAI for the common evolution loop.
 *  The original declaration is retained unchanged; this view defines no solver
 *  or additional physical state. Numerical dispatch remains internal.
 */
class ContinuousDynamics {
  public:
    using TimeType = Real;
    explicit ContinuousDynamics(const OrdinaryDifferentialAlgebraicEquation&);
    explicit ContinuousDynamics(const OrdinaryDifferentialAlgebraicInclusion&);
    ContinuousDynamics* clone() const;
    SizeType dimension() const;
    RealSpace state_space() const;
    RealSpace auxiliary_space() const;
    RealSpace state_auxiliary_space() const;
    const RealSpace& parameter_space() const;
    const BoxDomainType& constant_parameter_domain() const;
    bool has_explicit_time() const;

  private:
    std::variant<SharedPointer<const OrdinaryDifferentialAlgebraicEquation>,
                 SharedPointer<const OrdinaryDifferentialAlgebraicInclusion>>
        _declaration;
    const OrdinaryDifferentialAlgebraicEquation* equation() const;
    const OrdinaryDifferentialAlgebraicInclusion* inclusion() const;
    friend class ContinuousTimeEvolver<ValidatedTag>;
    friend class ContinuousTimeEvolver<ApproximateTag>;
};

} // namespace Ariadne
#endif
