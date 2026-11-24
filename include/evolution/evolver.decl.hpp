/***************************************************************************
 *            evolver.decl.hpp
 *
 *  Copyright  2007-20  Alberto Casagrande, Pieter Collins
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

/** @file evolution/evolver.decl.hpp
 *  @brief Explicit numerical paradigm declarations for evolution algorithms.
 */

#ifndef ARIADNE_EVOLVER_DECL_HPP
#define ARIADNE_EVOLVER_DECL_HPP

#include "foundation/paradigm.hpp"

namespace Ariadne {
template<class P>
class ContinuousTimeEvolver;
template<class P>
class ContinuousTimeEvolverConfiguration;
template<class P>
class DiscreteTimeEvolver;
template<class P>
class DiscreteTimeEvolverConfiguration;
template<class P>
class HybridDiscreteTimeEvolver;
template<class P>
class HybridDiscreteTimeEvolverConfiguration;
template<class P>
class HybridContinuousTimeEvolver;
template<class P>
class HybridContinuousTimeEvolverBase;
template<class P>
class HybridContinuousTimeEvolverBaseConfiguration;
template<class P>
class HybridContinuousTimeEvolverConfiguration;
template<class P>
class HybridContinuousTimeEvolverFactory;
template<class P>
class HybridContinuousTimeEvolverFactoryInterface;
template<class P>
class HybridContinuousTimeEvolverInterface;
} // namespace Ariadne
#endif
