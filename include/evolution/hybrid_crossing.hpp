/***************************************************************************
 *            evolution/hybrid_crossing.hpp
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

/** @file evolution/hybrid_crossing.hpp
 *  @brief Typed geometric and numerical crossing information.
 */

#ifndef ARIADNE_EVOLUTION_HYBRID_CROSSING_HPP
#define ARIADNE_EVOLUTION_HYBRID_CROSSING_HPP

#include <cstdint>
#include <variant>

#include "function/function_patch.hpp"
#include "geometry/box.hpp"

namespace Ariadne {

/** @brief Certified guard geometry; UNRESOLVED makes no degeneracy claim. */
enum class CrossingGeometry : std::uint8_t {
    UNRESOLVED,
    NEGATIVE,
    POSITIVE,
    INCREASING,
    DECREASING,
    CONVEX,
    CONCAVE
};

enum class CrossingEvidenceKind : std::uint8_t { NONE,
                                                 CROSSING_TIME,
                                                 GUARD_MAXIMUM };

/** @brief Lightweight provenance without retaining numerical patches. */
struct CrossingSummary {
    CrossingGeometry geometry = CrossingGeometry::UNRESOLVED;
    CrossingEvidenceKind evidence = CrossingEvidenceKind::NONE;
};

/** @brief A validated first-entry time map on the original state chart.
 *  The numerical producer certifies the root and first-entry properties.
 *  Construction rejects missing, empty-domain or unbounded maps.
 */
class CrossingTimeCertificate {
    ValidatedScalarMultivariateFunctionPatch _time;

  public:
    explicit CrossingTimeCertificate(const ValidatedScalarMultivariateFunctionPatch&);
    const ValidatedScalarMultivariateFunctionPatch& time() const;
};

/** @brief A certified concave guard maximum; it need not be a tangency.
 *  All three required data are owned on the original state chart. The producer
 *  certifies their mathematical meaning; construction checks finite bounds and
 *  matching state-chart dimensions. Critical times may precede the step.
 */
class GuardMaximumCertificate {
    ValidatedScalarMultivariateFunctionPatch _time;
    UpperIntervalType _guard_range;
    UpperBoxType _state_bounds;

  public:
    GuardMaximumCertificate(const ValidatedScalarMultivariateFunctionPatch&,
                            const UpperIntervalType& guard_range,
                            const UpperBoxType& state_bounds);
    const ValidatedScalarMultivariateFunctionPatch& time() const;
    const UpperIntervalType& guard_range() const;
    const UpperBoxType& state_bounds() const;
};

/** @brief Independent geometry and optional typed numerical evidence.
 *  Evidence is immutable and may be absent when an equation solve failed.
 *  No mutable tag can claim a certificate whose map was not supplied.
 */
class CrossingData {
    CrossingGeometry _geometry;
    std::variant<std::monostate, CrossingTimeCertificate, GuardMaximumCertificate> _evidence;

  public:
    explicit CrossingData(CrossingGeometry = CrossingGeometry::UNRESOLVED);
    CrossingData(CrossingGeometry, const CrossingTimeCertificate&);
    explicit CrossingData(const GuardMaximumCertificate&);
    CrossingGeometry geometry() const;
    const CrossingTimeCertificate* crossing_time() const;
    const GuardMaximumCertificate* guard_maximum() const;
    CrossingSummary summary() const;
};

OutputStream& operator<<(OutputStream&, CrossingGeometry);
OutputStream& operator<<(OutputStream&, CrossingEvidenceKind);
OutputStream& operator<<(OutputStream&, const CrossingSummary&);
OutputStream& operator<<(OutputStream&, const CrossingData&);

} // namespace Ariadne
#endif
