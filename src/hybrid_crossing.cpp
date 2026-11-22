/***************************************************************************
 *            hybrid_crossing.cpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#include "evolution/hybrid_crossing.hpp"

#include <stdexcept>

namespace Ariadne {
namespace {
void validate_geometry(CrossingGeometry geometry) {
    switch (geometry) {
        case CrossingGeometry::UNRESOLVED:
        case CrossingGeometry::NEGATIVE:
        case CrossingGeometry::POSITIVE:
        case CrossingGeometry::INCREASING:
        case CrossingGeometry::DECREASING:
        case CrossingGeometry::CONVEX:
        case CrossingGeometry::CONCAVE:
            return;
        default:
            throw std::invalid_argument("Unknown crossing geometry");
    }
}

void validate_time_map(const ValidatedScalarMultivariateFunctionPatch& time) {
    if (time.raw_pointer() == nullptr)
        throw std::invalid_argument("A crossing certificate requires a time map");
    if (time.domain().is_empty() || !time.domain().is_bounded() || cast_exact_interval(time.range()).is_empty() || !definitely(time.range().is_bounded()))
        throw std::invalid_argument("A crossing certificate requires a nonempty finite state chart and time range");
}
} // namespace

CrossingTimeCertificate::CrossingTimeCertificate(const ValidatedScalarMultivariateFunctionPatch& time)
    : _time(time) { validate_time_map(_time); }
const ValidatedScalarMultivariateFunctionPatch& CrossingTimeCertificate::time() const { return _time; }

GuardMaximumCertificate::GuardMaximumCertificate(const ValidatedScalarMultivariateFunctionPatch& time,
                                                 const UpperIntervalType& guard_range, const UpperBoxType& state_bounds)
    : _time(time), _guard_range(guard_range), _state_bounds(state_bounds) {
    validate_time_map(_time);
    if (_state_bounds.dimension() != _time.argument_size() || cast_exact_box(_state_bounds).is_empty() ||
        !definitely(_state_bounds.is_bounded()) || cast_exact_interval(_guard_range).is_empty() || !definitely(_guard_range.is_bounded()))
        throw std::invalid_argument("A guard maximum requires matching finite state bounds and a finite guard range");
}
const ValidatedScalarMultivariateFunctionPatch& GuardMaximumCertificate::time() const { return _time; }
const UpperIntervalType& GuardMaximumCertificate::guard_range() const { return _guard_range; }
const UpperBoxType& GuardMaximumCertificate::state_bounds() const { return _state_bounds; }

CrossingData::CrossingData(CrossingGeometry geometry) : _geometry(geometry) { validate_geometry(_geometry); }
CrossingData::CrossingData(CrossingGeometry geometry, const CrossingTimeCertificate& time)
    : _geometry(geometry), _evidence(time) {
    validate_geometry(_geometry);
    if (_geometry != CrossingGeometry::INCREASING && _geometry != CrossingGeometry::CONVEX &&
        _geometry != CrossingGeometry::CONCAVE)
        throw std::invalid_argument("A first-entry certificate requires compatible crossing geometry");
}
CrossingData::CrossingData(const GuardMaximumCertificate& maximum)
    : _geometry(CrossingGeometry::CONCAVE), _evidence(maximum) {}
CrossingGeometry CrossingData::geometry() const { return _geometry; }
const CrossingTimeCertificate* CrossingData::crossing_time() const { return std::get_if<CrossingTimeCertificate>(&_evidence); }
const GuardMaximumCertificate* CrossingData::guard_maximum() const { return std::get_if<GuardMaximumCertificate>(&_evidence); }
CrossingSummary CrossingData::summary() const {
    return {_geometry, crossing_time() != nullptr   ? CrossingEvidenceKind::CROSSING_TIME
                       : guard_maximum() != nullptr ? CrossingEvidenceKind::GUARD_MAXIMUM
                                                    : CrossingEvidenceKind::NONE};
}

OutputStream& operator<<(OutputStream& os, CrossingGeometry geometry) {
    switch (geometry) {
        case CrossingGeometry::UNRESOLVED:
            return os << "UNRESOLVED";
        case CrossingGeometry::NEGATIVE:
            return os << "NEGATIVE";
        case CrossingGeometry::POSITIVE:
            return os << "POSITIVE";
        case CrossingGeometry::INCREASING:
            return os << "INCREASING";
        case CrossingGeometry::DECREASING:
            return os << "DECREASING";
        case CrossingGeometry::CONVEX:
            return os << "CONVEX";
        case CrossingGeometry::CONCAVE:
            return os << "CONCAVE";
        default:
            return os << "unknown";
    }
}
OutputStream& operator<<(OutputStream& os, CrossingEvidenceKind kind) {
    switch (kind) {
        case CrossingEvidenceKind::NONE:
            return os << "NONE";
        case CrossingEvidenceKind::CROSSING_TIME:
            return os << "CROSSING_TIME";
        case CrossingEvidenceKind::GUARD_MAXIMUM:
            return os << "GUARD_MAXIMUM";
        default:
            return os << "unknown";
    }
}
OutputStream& operator<<(OutputStream& os, const CrossingSummary& crossing) {
    return os << "{geometry=" << crossing.geometry << ",evidence=" << crossing.evidence << "}";
}
OutputStream& operator<<(OutputStream& os, const CrossingData& crossing) {
    os << "{geometry=" << crossing.geometry() << ",evidence=" << crossing.summary().evidence;
    if (const auto* time = crossing.crossing_time())
        os << ",crossing_time_range=" << time->time().range() << ",crossing_time_error=" << time->time().error();
    if (const auto* maximum = crossing.guard_maximum())
        os << ",critical_time_range=" << maximum->time().range() << ",critical_time_error=" << maximum->time().error()
           << ",guard_range_at_critical_time=" << maximum->guard_range() << ",state_bounds_at_critical_time=" << maximum->state_bounds();
    return os << "}";
}
} // namespace Ariadne
