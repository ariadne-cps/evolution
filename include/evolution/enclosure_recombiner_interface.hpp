/***************************************************************************
 *            evolution/enclosure_recombiner_interface.hpp
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

/** @file evolution/enclosure_recombiner_interface.hpp
 *  @brief Contract and provenance for continuous enclosure recombination.
 */

#ifndef ARIADNE_ENCLOSURE_RECOMBINER_INTERFACE_HPP
#define ARIADNE_ENCLOSURE_RECOMBINER_INTERFACE_HPP

#include <optional>

#include "evolution/enclosure.hpp"
#include "utility/writable.hpp"

namespace Ariadne {

/** @brief Numerical limits supplied for one recombination call.
 *  @details Values belong to the current evolution configuration; recombiners
 *  must not retain a reference to this context or modify global configuration.
 */
class EnclosureRecombinationContext {
    std::optional<ExactDouble> _maximum_enclosure_radius;

  public:
    /** @brief Request recombination without a radius limit. */
    EnclosureRecombinationContext() = default;
    /** @brief Request a finite, nonnegative maximum radius. */
    explicit EnclosureRecombinationContext(ExactDouble maximum_enclosure_radius);
    /** @brief Return this call's radius limit, or no value when unrestricted. */
    const std::optional<ExactDouble>& maximum_enclosure_radius() const { return _maximum_enclosure_radius; }
};

/** @brief One output enclosure and its contributing original input indices.
 *  @details A source may contribute to multiple outputs when a strategy
 *  subdivides it. Indices refer to the input list before empty-source removal.
 */
struct EnclosureRecombinationOutput {
    Enclosure enclosure;           //!< An upper enclosure for the represented source portions.
    List<SizeType> source_indices; //!< Original sources represented by this output.
};

/** @brief Recombined collection, provenance and the outcome of a radius request. */
struct EnclosureRecombinationResult {
    List<EnclosureRecombinationOutput> outputs; //!< Outputs jointly covering the input union.
    List<SizeType> discarded_sources;           //!< Original indices discarded only after proving emptiness.
    Bool radius_limit_satisfied = true;         //!< All outputs are certified within the requested radius, or no limit was requested.
};

/** @brief Strategy constructing an upper covering collection of continuous enclosures.
 *  @details Implementations preserve the union of the input sets and do not
 *  modify the sources. Each output records its contributing original indices;
 *  only proven-empty inputs may be listed as discarded. An empty input produces
 *  an empty result. Failure to meet a requested radius is reported explicitly,
 *  without sacrificing coverage. The radius is measured by Enclosure::radius().
 *  Work and diagnostics belong to each call, not to mutable shared strategy state.
 */
class EnclosureRecombinerInterface : public virtual WritableInterface {
  public:
    /** @brief Destroy a recombiner through its interface. */
    ~EnclosureRecombinerInterface() override = default;
    /** @brief Create an independently owned copy of this strategy. */
    virtual EnclosureRecombinerInterface* clone() const = 0;
    /** @brief Recombine sources under this call's numerical limits. */
    virtual EnclosureRecombinationResult recombine(const List<Enclosure>& enclosures,
                                                   const EnclosureRecombinationContext& context = {}) const = 0;
};

} // namespace Ariadne

#endif // ARIADNE_ENCLOSURE_RECOMBINER_INTERFACE_HPP
