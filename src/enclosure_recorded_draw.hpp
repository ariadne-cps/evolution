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
#include <vector>

#include "rendering/recording_canvas.hpp"
#include "threading/workload.hpp"

namespace Ariadne {

// Compute individual enclosure drawings independently, then replay them in
// their original order onto the sole mutable destination canvas. Keep only one
// bounded batch of recordings in memory, irrespective of orbit size.
template<class EnclosureList, class Draw>
Void draw_enclosures_recorded(CanvasInterface& canvas, const EnclosureList& enclosures, Draw draw) {
    constexpr SizeType batch_size = 256u;
    for (SizeType offset = 0; offset < enclosures.size(); offset += batch_size) {
        const SizeType count = std::min(batch_size, enclosures.size() - offset);
        std::vector<RecordingCanvas> recordings(count);
        StaticWorkload<SizeType> workload([&](const SizeType& index) {
            draw(enclosures[offset + index], recordings[index]);
        });
        for (SizeType index = 0; index < count; ++index)
            workload.append(index);
        workload.process();
        for (const auto& recording: recordings)
            recording.replay(canvas);
    }
}

} // namespace Ariadne

#endif // ARIADNE_EVOLUTION_ENCLOSURE_RECORDED_DRAW_HPP
