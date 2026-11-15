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

#include <cstdint>

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

/** @brief Opt-in, per-call reference recombination phase durations in nanoseconds.
 *  @details Phases are exclusive of each other but not exhaustive: control
 *  flow, provenance, hull merging, and timer overhead remain unaccounted for.
 *  A profiled call clears these counters before processing any sources.
 */
struct EnclosureRecombinationProfile {
    std::uint64_t filter_ns = 0;          //!< Validate and filter input sources, including emptiness checks.
    std::uint64_t normalise_ns = 0;       //!< Normalise source parameterisations and model constraints.
    std::uint64_t image_setup_ns = 0;     //!< Construct joint residual images for each source.
    std::uint64_t image_bounds_ns = 0;    //!< Calculate image ranges and tighten by source predicates.
    std::uint64_t graph_setup_ns = 0;     //!< Assemble and prepare constraint graphs.
    std::uint64_t graph_reduce_ns = 0;    //!< Run the prepared constraint solver on each graph.
    std::uint64_t output_assembly_ns = 0; //!< Build the final residual parameterisation and predicates.
    std::uint64_t output_reduce_ns = 0;   //!< Contract the final recombined enclosure.
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
    /** @brief Run the same reference algorithm while collecting optional per-call phase times.
     *  @details A null profile disables all clock reads; the strategy keeps no profiling state.
     *  This overload is reference-specific and does not change the polymorphic interface.
     */
    EnclosureRecombinationResult recombine(const List<Enclosure>& enclosures,
                                           const EnclosureRecombinationContext& context,
                                           EnclosureRecombinationProfile* profile) const;
    /** @brief Write the strategy name. */
    OutputStream& _write(OutputStream& os) const override;
};

} // namespace Ariadne

#endif // ARIADNE_ENCLOSURE_RECOMBINER_HPP
