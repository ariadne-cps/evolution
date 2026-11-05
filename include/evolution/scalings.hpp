/***************************************************************************
 *            dynamics/scalings.hpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file dynamics/scalings.hpp
 *  \brief Scalings for real variables.
 */

#ifndef ARIADNE_SCALINGS_HPP
#define ARIADNE_SCALINGS_HPP

#include <iostream>
#include <map>

#include "utility/container.hpp"
#include "utility/stlio.hpp"

#include "numeric/builtin.hpp"
#include "paving/grid.hpp"
#include "symbolic/variable.hpp"
#include "symbolic/space.hpp"

namespace Ariadne {

typedef void Void;
typedef std::ostream OutputStream;

//! \ingroup DynamicsModule
//! \brief A class which defines the state space grid to use given the continuous state variables \a spc.
class Scalings {
    ExactDouble _default_scaling;
    Map<Identifier,ExactDouble> _scalings;
  public:
    Scalings(ApproximateDouble default_scaling)
        : _default_scaling(cast_exact(default_scaling)), _scalings() { ARIADNE_PRECONDITION(decide(default_scaling>0)); }
    Scalings(ExactDouble default_scaling, Map<RealVariable,ExactDouble> const& scalings);
    Void set_scaling(RealVariable const& var, ApproximateDouble scal) {
        ARIADNE_PRECONDITION(decide(scal>0)); _scalings[var.name()]=cast_exact(scal); }
    ExactDouble scaling(const RealVariable& var) const {
        return (this->_scalings.has_key(var.name())) ? this->_scalings[var.name()] : this->_default_scaling; }
    Grid grid(RealSpace const& spc) const;
    friend OutputStream& operator<<(OutputStream& os, Scalings const& s) {
        return os << "Scalings( " << s._scalings << " )"; }
};

inline Scalings::Scalings(ExactDouble default_scaling, Map<RealVariable,ExactDouble> const& scalings)
    : _default_scaling(default_scaling), _scalings()
{
    ARIADNE_PRECONDITION(default_scaling>0);
    for(auto var_scal : scalings) {
        ARIADNE_PRECONDITION(var_scal.second>0);
        this->_scalings.insert(var_scal.first.name(),var_scal.second);
    }
}

inline Grid Scalings::grid(RealSpace const& spc) const {
    Vector<ExactDouble> lengths(spc.dimension(),1.0_x);
    for(SizeType i=0; i<spc.dimension(); ++i) { lengths[i]=this->scaling(spc[i]); }
    return Grid(lengths);
}

} // namespace Ariadne

#endif // ARIADNE_SCALINGS_HPP

