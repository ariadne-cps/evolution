/***************************************************************************
 *            evolution/enclosure_recombiner.hpp
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

/** @file evolution/enclosure_recombiner.hpp
 *  @brief Polymorphic handle and reference continuous enclosure recombiner.
 */

#ifndef ARIADNE_ENCLOSURE_RECOMBINER_HPP
#define ARIADNE_ENCLOSURE_RECOMBINER_HPP

#include "evolution/enclosure_recombiner_interface.hpp"
#include "utility/handle.hpp"

namespace Ariadne {

/** @brief Shared owning handle to an immutable recombination strategy. */
class EnclosureRecombiner : public Handle<const EnclosureRecombinerInterface> {
  public:
    /** @brief Own a polymorphic clone of the supplied implementation. */
    explicit EnclosureRecombiner(const EnclosureRecombinerInterface& recombiner);
    /** @brief Dispatch recombination and return its outputs and provenance. */
    EnclosureRecombinationResult recombine(const List<Enclosure>& enclosures,
                                           const EnclosureRecombinationContext& context = {}) const;
    /** @brief Write the represented strategy. */
    friend OutputStream& operator<<(OutputStream& os, const EnclosureRecombiner& recombiner) {
        return recombiner._ptr->_write(os);
    }
};

/** @brief Recombine all surviving sources around the first source's nonlinear map.
 *  @details Contracts each source independently and encloses state/clock residuals
 *  in extra parameters. Reference predicates are relaxed over the source union.
 *  At most one output is produced; this strategy does not subdivide to meet a
 *  requested radius. Such a request is checked and reported in the result.
 */
class ReferenceEnclosureRecombiner : public EnclosureRecombinerInterface {
  public:
    /** @brief Create an independently owned reference recombiner. */
    ReferenceEnclosureRecombiner* clone() const override;
    /** @brief Construct the reference upper enclosure and record source provenance.
     *  @pre Sources have initialised models and equal state dimensions; surviving
     *       sources have bounded reduced domains and identical symbolic auxiliary maps.
     */
    EnclosureRecombinationResult recombine(const List<Enclosure>& enclosures,
                                           const EnclosureRecombinationContext& context = {}) const override;
    /** @brief Write the strategy name. */
    OutputStream& _write(OutputStream& os) const override;
};

} // namespace Ariadne

#endif // ARIADNE_ENCLOSURE_RECOMBINER_HPP
