/***************************************************************************
 *            carried_expansion_diagnostics.hpp
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

/** @file src/carried_expansion_diagnostics.hpp
 *  @brief Diagnostics of carried Taylor expansions.
 */

#pragma once

#include "function/taylor_function.hpp"

namespace Ariadne::EvolutionDetail {

/** @brief Print a snapshot of the carried Taylor expansion to the given stream.
 *  An explicit output stream avoids redirecting process-wide diagnostics.
 */
Void print_carried_expansion_snapshot(
    OutputStream& output, const String& stage, SizeType step,
    const ValidatedVectorMultivariateTaylorFunctionModelDP& mapping,
    const Sweeper<FloatDP>& sweeper);

} // namespace Ariadne::EvolutionDetail
