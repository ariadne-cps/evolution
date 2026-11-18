/***************************************************************************
 *            enclosure_recorded_draw.hpp
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 */

#ifndef ARIADNE_EVOLUTION_ENCLOSURE_RECORDED_DRAW_HPP
#define ARIADNE_EVOLUTION_ENCLOSURE_RECORDED_DRAW_HPP

#include <algorithm>
#include <atomic>
#include <vector>

#include "rendering/recording_canvas.hpp"
#include "threading/workload.hpp"

namespace Ariadne {

// Keep drawing progress over the entire collection, independently of the
// batches used to bound recording memory.
class RecordedDrawProgress {
  public:
    explicit RecordedDrawProgress(SizeType total) : _total(total) {}

    Void complete_one() { _completed.fetch_add(1u, std::memory_order_relaxed); }

    Void update(ProgressIndicator& indicator) const {
        indicator.update_final(static_cast<double>(_total));
        indicator.update_current(static_cast<double>(_completed.load(std::memory_order_relaxed)));
    }

  private:
    SizeType _total;
    std::atomic<SizeType> _completed{0u};
};

// Compute individual enclosure drawings independently, then replay them in
// their original order onto the sole mutable destination canvas. Keep only one
// bounded batch of recordings in memory, irrespective of orbit size.
template<class EnclosureList, class Draw>
Void draw_enclosures_recorded(CanvasInterface& canvas, const EnclosureList& enclosures, Draw draw) {
    constexpr SizeType batch_size = 256u;
    const SizeType total = enclosures.size();
    RecordedDrawProgress progress(total);
    std::vector<RecordingCanvas> recordings;
    SizeType offset = 0u;

    // One workload retains the indicator and completed count across batches.
    using Workload = DynamicWorkload<SizeType>;
    Workload workload(
        [&](const SizeType&, SharedPointer<ProgressIndicator> indicator) {
            progress.update(*indicator);
        },
        [&](Workload::Access&, const SizeType& index) {
            // The sentinel only refreshes the final 100% after all replay.
            if (index == total)
                return;
            draw(enclosures[offset + index], recordings[index]);
            progress.complete_one();
        });

    for (; offset < total; offset += batch_size) {
        const SizeType count = std::min(batch_size, total - offset);
        recordings = std::vector<RecordingCanvas>(count);
        for (SizeType index = 0; index < count; ++index)
            workload.append(index);
        workload.process();
        for (const auto& recording: recordings)
            recording.replay(canvas);
    }
    if (total != 0u) {
        // Progress callbacks run before tasks, so explicitly report completion.
        workload.append(total);
        workload.process();
    }
}

} // namespace Ariadne

#endif // ARIADNE_EVOLUTION_ENCLOSURE_RECORDED_DRAW_HPP
