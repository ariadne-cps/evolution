/***************************************************************************
 *            hybrid_orbit.hpp
 *
 *  Copyright  2007-20  Pieter Collins
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

/** @file evolution/hybrid_orbit.hpp
 *  @brief Discrete and continuous orbit collections for hybrid systems.
 */

#ifndef ARIADNE_HYBRID_ORBIT_HPP
#define ARIADNE_HYBRID_ORBIT_HPP

#include <iostream>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "evolution/hybrid_enclosure.hpp"
#include "evolution/orbit.hpp"
#include "hybrid/hybrid_set.hpp"
#include "labelled/graphics_interface.hpp"
#include "numeric/numeric.hpp"
#include "set/function_set.hpp"

namespace Ariadne {

typedef FloatDP Time;

template<class ES>
class Orbit;

template<class BS>
class ListSet;
template<class X>
class Point;

/** @brief Grid represents discrete and continuous orbit collections for hybrid systems. */
class Grid;
/** @brief GridCell represents discrete and continuous orbit collections for hybrid systems. */
class GridCell;
/** @brief BinaryTreePaving represents discrete and continuous orbit collections for hybrid systems. */
class BinaryTreePaving;
/** @brief HybridGrid represents discrete and continuous orbit collections for hybrid systems. */
class HybridGrid;
/** @brief HybridGridCell represents discrete and continuous orbit collections for hybrid systems. */
class HybridGridCell;
/** @brief HybridPaving represents discrete and continuous orbit collections for hybrid systems. */
class HybridPaving;
/** @brief HybridTime represents discrete and continuous orbit collections for hybrid systems. */
class HybridTime;

/** @brief HybridAutomatonInterface represents discrete and continuous orbit collections for hybrid systems. */
class HybridAutomatonInterface;

/** @brief DiscreteLocation represents discrete and continuous orbit collections for hybrid systems. */
class DiscreteLocation;

typedef HybridBasicSet<InterpolatedCurve> HybridInterpolatedCurve;

/** @brief Apply the requested operation. */
template<class ES>
OutputStream& operator<<(OutputStream&, const Orbit<ES>&);

template<>
/** @brief Orbit represents discrete and continuous orbit collections for hybrid systems. */
class Orbit<HybridApproximatePoint>
    : public HybridDrawableInterface {
  public:
    /** @brief Construct the object. */
    Orbit(const HybridApproximatePoint& hpt);
    /** @brief Construct the object. */
    Orbit(List<HybridInterpolatedCurve> crvs);
    /** @brief Perform the specified evolution operation. */
    Void insert(HybridTime ht, const HybridApproximatePoint& hpt);
    /** @brief Return the requested evolution property. */
    SizeType size() const;
    /** @brief Perform the specified evolution operation. */
    const InterpolatedCurve& curve(Nat m) const;
    const List<HybridInterpolatedCurve>& curves() const { return *this->_curves_ptr; }
    /** @brief Perform the specified evolution operation. */
    Void draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d& v) const;

  private:
    SharedPointer<List<HybridInterpolatedCurve>> _curves_ptr;
};

template<>
/** @brief Orbit represents discrete and continuous orbit collections for hybrid systems. */
class Orbit<Vector<HybridApproximatePoint>>
    : public HybridDrawableInterface {
  public:
    /** @brief Construct the object. */
    Orbit(const Vector<HybridApproximatePoint>& vhpt);
    /** @brief Perform the specified evolution operation. */
    Void insert(HybridTime ht, const HybridApproximatePoint& hpt, SizeType setNumber);
    //! The number of independent trajectories represented by this orbit.
    SizeType number_of_trajectories() const { return _curves_ptr->size(); }
    /** @brief Return the requested evolution property. */
    SizeType size(SizeType setNumber);
    /** @brief Perform the specified evolution operation. */
    const InterpolatedCurve& curve(Nat setNumber, Nat step);
    const List<HybridInterpolatedCurve>& curves(Nat setNumber) { return this->_curves_ptr->at(setNumber); }
    /** @brief Perform the specified evolution operation. */
    Void draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d& v) const;

  private:
    SharedPointer<Vector<List<HybridInterpolatedCurve>>> _curves_ptr;
};

/** @brief Perform the specified evolution operation. */
Orbit<HybridApproximatePoint> extend_auxiliary(const Orbit<HybridApproximatePoint>& horb, const HybridAutomatonInterface& ha);

template<>
OutputStream&
    /** @brief Apply the requested operation. */
    operator<<(OutputStream& os, const Orbit<HybridApproximatePoint>& orb);

template<>
/** @brief Orbit represents discrete and continuous orbit collections for hybrid systems. */
class Orbit<HybridEnclosure>
    : public HybridDrawableInterface {
  public:
    typedef HybridEnclosure EnclosureType;
    typedef ListSet<HybridEnclosure> EnclosureListType;

    Orbit() {}
    Orbit(const EnclosureType& set) : _initial(set) {}
    Void adjoin_reach(const EnclosureType& set) { this->_reach.adjoin(set); }
    Void adjoin_intermediate(const EnclosureType& set) { this->_intermediate.adjoin(set); }
    Void adjoin_final(const EnclosureType& set) { this->_final.adjoin(set); }

    Void adjoin_reach(const EnclosureListType& set) { this->_reach.adjoin(set); }
    Void adjoin_intermediate(const EnclosureListType& set) { this->_intermediate.adjoin(set); }
    Void adjoin_final(const EnclosureListType& set) { this->_final.adjoin(set); }

    const EnclosureType& initial() const { return this->_initial; }
    const EnclosureListType& reach() const { return this->_reach; }
    const EnclosureListType& intermediate() const { return this->_intermediate; }
    const EnclosureListType& final() const { return this->_final; }

    /** @brief Perform the specified evolution operation. */
    Void draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d& v) const;

  private:
    EnclosureType _initial;
    EnclosureListType _reach;
    EnclosureListType _intermediate;
    EnclosureListType _final;
};

template<>
OutputStream&
    /** @brief Apply the requested operation. */
    operator<<(OutputStream& os, const Orbit<HybridEnclosure>& orb);

} // namespace Ariadne

#endif // ARIADNE_HYBRID_ORBIT_HPP
