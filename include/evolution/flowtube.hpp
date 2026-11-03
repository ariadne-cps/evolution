/***************************************************************************
 *            dynamics/flowtube.hpp
 *
 *  Copyright  2007-20  Pieter Collins
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

/*! \file dynamics/flowtube.hpp
 *  \brief Flow tubes for dynamic systems
 */

#ifndef ARIADNE_FLOWTUBE_HPP
#define ARIADNE_FLOWTUBE_HPP


#include "function/function_model.hpp"
#include "set/function_set.hpp"

namespace Ariadne {

template<class ES> class Orbit;

using ValidatedImageSet = ValidatedConstrainedImageSet;

class FlowTube {
    ValidatedVectorMultivariateFunctionModelDP _phi;
  public:
    //! Construct a flow map phi(x,t) over its spatial and temporal domain.
    explicit FlowTube(ValidatedVectorMultivariateFunctionModelDP const& phi) : _phi(phi) {
        ARIADNE_PRECONDITION(phi.argument_size()==phi.result_size()+1u);
    }
    SizeType number_of_state_variables() const { return _phi.result_size(); }
    template<class X, class T> requires Same<X,T> decltype(auto) operator()(Vector<X> const& x, Scalar<T> const& t) const {
        return _phi(join(x,t)); }
    decltype(auto) space_domain() const { SizeType n=this->number_of_state_variables(); return _phi.domain()[range(0,n)]; }
    auto time_range() const { SizeType n=this->number_of_state_variables(); return _phi.domain()[n]; }
    auto initial_time() const { return time_range().lower_bound(); }
    auto final_time() const { return time_range().upper_bound(); }
    ValidatedVectorMultivariateFunctionModelDP reach_function() const { return _phi; }
    ValidatedVectorMultivariateFunctionModelDP final_function() const {
        return partial_evaluate(_phi,number_of_state_variables(),ValidatedNumber(final_time())); }
    ValidatedImageSet initial_set() const {
        auto initial=partial_evaluate(_phi,number_of_state_variables(),ValidatedNumber(initial_time()));
        return ValidatedImageSet(space_domain(),initial); }
    ValidatedImageSet reach_set() const { return ValidatedImageSet(_phi.domain(),_phi); }
    ValidatedImageSet final_set() const { return ValidatedImageSet(space_domain(),final_function()); }
};

template<class E, class R=E> struct OrbitStep {
    OrbitStep* _previous;
    E _initial; R _intermediate;
    List<OrbitStep> _next;
    List<E> _final;
};

template<> class Orbit<FlowTube> {
    using ES=FlowTube;
    OrbitStep<ES> _starting;
    List<OrbitStep<ES>*> _working;
  public:
    Orbit(const ES& set);
};

} // namespace Ariadne

#endif // ARIADNE_FLOWTUBE_HPP
