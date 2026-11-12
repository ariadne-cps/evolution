/***************************************************************************
 *            orbit.hpp
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

/** @file evolution/orbit.hpp
 *  @brief Orbit containers for continuous and discrete systems.
 */

#ifndef ARIADNE_ORBIT_HPP
#define ARIADNE_ORBIT_HPP

#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

#include "evolution/enclosure.hpp"
#include "labelled/curve.hpp"
#include "labelled/graphics_interface.hpp"
#include "labelled/point.hpp"
#include "numeric/numeric.hpp"
#include "set/curve.hpp"
#include "set/function_set.hpp"
#include "set/list_set.hpp"
#include "threading/using.hpp"

namespace Ariadne {

#ifdef DOXYGEN
//! \brief Class for storing evolution data.
template<class E>
class Orbit {
  public:
    //! The type used to store a single enclosure set.
    typedef E EnclosureType;
    //! The type used to store a list of enclosure sets.
    typedef EL EnclosureListType;
    //! The initial set of the orbit.
    const EnclosureType& initial() const;
    //! The set of points reached by evolution from the given initial set over the evolution time.
    const EnclosureListType& reach() const;
    //! The set of points reached by evolution from the given initial set at the final evolution time.
    const EnclosureListType& final() const;
};
#endif

template<class ES>
class Orbit;
/** @brief Apply the requested operation. */
template<class ES>
OutputStream& operator<<(OutputStream&, const Orbit<ES>&);

template<class BS>
class ListSet;

template<class X>
class Point;
/** @brief InterpolatedCurve represents orbit containers for continuous and discrete systems. */
class InterpolatedCurve;

template<class F>
/** @brief Orbit represents orbit containers for continuous and discrete systems. */
class Orbit<Point<F>> {
  public:
    /** @brief Construct the object. */
    Orbit(const Point<F>& pt);
    /** @brief Perform the specified evolution operation. */
    Void insert(F t, const Point<F>& pt);
    const InterpolatedCurve& curve() const { return *this->_curve; }

  private:
    std::shared_ptr<InterpolatedCurve> _curve;
};

template<class F>
/** @brief Orbit represents orbit containers for continuous and discrete systems. */
class Orbit<Point<Approximation<F>>> {
  public:
    Orbit(const Point<Approximation<F>>& pt) : _curve(new InterpolatedCurve(0, pt)) {}
    Void insert(F t, const Point<Approximation<F>>& pt) { this->_curve->insert(t, pt); }
    const InterpolatedCurve& curve() const { return *this->_curve; }

  private:
    std::shared_ptr<InterpolatedCurve> _curve;
};

template<class F>
/** @brief Orbit represents orbit containers for continuous and discrete systems. */
class Orbit<List<Point<Approximation<F>>>> {
  public:
    /** @brief Construct the object. */
    Orbit(const List<Point<Approximation<F>>>& ptLst);
    /** @brief Perform the specified evolution operation. */
    Void insert(F t, const Point<Approximation<F>>& pt, SizeType curveNumber);
    const List<InterpolatedCurve>& curves() { return *_curves; }
    const InterpolatedCurve& curve(SizeType i) { return _curves->at(i); }

  private:
    std::shared_ptr<List<InterpolatedCurve>> _curves;
};

template<class F>
/** @brief Orbit represents orbit containers for continuous and discrete systems. */
class Orbit<LabelledPoint<Approximation<F>>>
    : public LabelledDrawable2dInterface {
  public:
    Orbit(const LabelledPoint<Approximation<F>>& pt) : _curve(new LabelledInterpolatedCurve(pt)) {}
    Void insert(F t, const LabelledPoint<Approximation<F>>& pt) { this->_curve->insert(t, pt); }
    const LabelledInterpolatedCurve& curves() const { return *this->_curve; }
    virtual Void draw(CanvasInterface& canvas, const Variables2d& axes) const override { _curve->draw(canvas, axes); }
    virtual LabelledDrawable2dInterface* labelled_clone() const override { return new Orbit<LabelledPoint<Approximation<F>>>(*this); }

  private:
    std::shared_ptr<LabelledInterpolatedCurve> _curve;
};

template<class F>
/** @brief Orbit represents orbit containers for continuous and discrete systems. */
class Orbit<List<LabelledPoint<Approximation<F>>>>
    : public LabelledDrawable2dInterface {
  public:
    /** @brief Construct the object. */
    Orbit(const List<LabelledPoint<Approximation<F>>>& ptLst);
    /** @brief Perform the specified evolution operation. */
    Void insert(F t, const LabelledPoint<Approximation<F>>& pt, SizeType curveNumber);
    const List<LabelledInterpolatedCurve>& curves() { return *_curves; }
    const LabelledInterpolatedCurve& curve(SizeType i) { return _curves->at(i); }
    /** @brief Perform the specified evolution operation. */
    virtual Void draw(CanvasInterface& canvas, const Variables2d& axes) const override {
        for (SizeType i = 0; i < _curves->size(); i++) {
            _curves->at(i).draw(canvas, axes);
        }
    }
    virtual LabelledDrawable2dInterface* labelled_clone() const override { return new Orbit<List<LabelledPoint<Approximation<F>>>>(*this); }

  private:
    std::shared_ptr<List<LabelledInterpolatedCurve>> _curves;
};

template<>
/** @brief Orbit represents orbit containers for continuous and discrete systems. */
class Orbit<Enclosure> {
    typedef Enclosure ES;
    typedef ListSet<ES> ESL;

  public:
    typedef ES EnclosureType;
    typedef ListSet<ES> EnclosureListType;

    Orbit(const ES& set) : _initial(set) {}
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

  private:
    ES _initial;
    ESL _reach;
    ESL _intermediate;
    ESL _final;
};

template<>
/** @brief Orbit represents orbit containers for continuous and discrete systems. */
class Orbit<LabelledEnclosure>
    : public LabelledDrawable2dInterface {
    typedef LabelledEnclosure ES;
    typedef ListSet<LabelledEnclosure> ESL;

  public:
    typedef ES EnclosureType;
    typedef ESL EnclosureListType;

    Orbit(const ES& set) : _initial(set) {}
    virtual Void adjoin_reach(const EnclosureType& set) { this->_reach.adjoin(set); }
    virtual Void adjoin_intermediate(const EnclosureType& set) { this->_intermediate.adjoin(set); }
    virtual Void adjoin_final(const EnclosureType& set) { this->_final.adjoin(set); }

    virtual Void adjoin_reach(const EnclosureListType& set) { this->_reach.adjoin(set); }
    virtual Void adjoin_intermediate(const EnclosureListType& set) { this->_intermediate.adjoin(set); }
    virtual Void adjoin_final(const EnclosureListType& set) { this->_final.adjoin(set); }

    const EnclosureType& initial() const { return this->_initial; }
    const EnclosureListType& reach() const { return this->_reach; }
    const EnclosureListType& intermediate() const { return this->_intermediate; }
    const EnclosureListType& final() const { return this->_final; }

    virtual Void draw(CanvasInterface& canvas, const Variables2d& axes) const override { _reach.draw(canvas, axes); }
    virtual LabelledDrawable2dInterface* labelled_clone() const override { return new Orbit<LabelledEnclosure>(*this); }

  private:
    ES _initial;
    ESL _reach;
    ESL _intermediate;
    ESL _final;
};

/** @brief Apply the requested operation. */
template<class ES>
OutputStream& operator<<(OutputStream& os, const Orbit<ES>& orb);

template<class ES>
OutputStream&
    operator<<(OutputStream& os, const Orbit<ES>& orb) {
    os << "Orbit(\n  initial=" << orb.initial()
       << "\n  intermediate=" << orb.intermediate()
       << "\n  reach=" << orb.reach()
       << "\n  final=" << orb.final()
       << ")\n";
    return os;
}

/** @brief Perform the specified evolution operation. */
template<class ES>
Void draw(FigureInterface& figure, const Orbit<ES>& orbit) {
    for (const auto& enclosure: orbit.reach()) {
        draw(figure, enclosure);
    }
    /** @brief Perform the specified evolution operation. */
    draw(figure, orbit.initial());
    for (const auto& enclosure: orbit.final()) {
        draw(figure, enclosure);
    }
}

/** @brief Apply the requested operation. */
template<class ES>
FigureInterface& operator<<(FigureInterface& figure, const Orbit<ES>& orbit) {
    draw(figure, orbit);
    return figure;
}

} // namespace Ariadne

#endif // ARIADNE_ORBIT_HPP
