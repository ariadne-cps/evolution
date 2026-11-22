/***************************************************************************
 *            evolution/located_projection.hpp
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

/** @file evolution/located_projection.hpp
 *  @brief Project continuous values selected from a hybrid location.
 */

#ifndef ARIADNE_EVOLUTION_LOCATED_PROJECTION_HPP
#define ARIADNE_EVOLUTION_LOCATED_PROJECTION_HPP

#include "hybrid/value.hpp"
#include "rendering/rendering_interface.hpp"

namespace Ariadne {

// The caller provides the continuous variable ordering for a Located<T>.
// No hybrid-aware figure or paving-specific drawing interface is involved.
template<class T, class Visitor>
Void with_projected_location(const Located<T>& item, const RealSpace& space,
                             const RealVariable& x, const RealVariable& y,
                             Visitor&& visitor) {
    with_location(item, space, [&](const RealSpace& selected, const T& value) {
        const Projection2d projection(selected.dimension(),
                                      selected.index(x), selected.index(y));
        visitor(projection, value);
    });
}

// Select an existing value, without lazily materialising another location.
template<class T, class Visitor>
Void with_projected_location(const Hybrid<T>& family, const DiscreteLocation& location,
                             const RealSpace& space,
                             const RealVariable& x, const RealVariable& y,
                             Visitor&& visitor) {
    with_location(family, location, space,
        [&](const RealSpace& selected, const T& value) {
            const Projection2d projection(selected.dimension(),
                                          selected.index(x), selected.index(y));
            visitor(projection, value);
        });
}

// Use the location-dependent RealSpace stored by the hybrid family.
template<class T, class Visitor>
Void with_projected_location(const Hybrid<T>& family, const DiscreteLocation& location,
                             const RealVariable& x, const RealVariable& y,
                             Visitor&& visitor) {
    with_location(family, location,
        [&](const RealSpace& selected, const T& value) {
            const Projection2d projection(selected.dimension(),
                                          selected.index(x), selected.index(y));
            visitor(projection, value);
        });
}

} // namespace Ariadne

#endif // ARIADNE_EVOLUTION_LOCATED_PROJECTION_HPP
