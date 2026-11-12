/***************************************************************************
 *            log_verbosity_scope.hpp
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

/** @file tests/log_verbosity_scope.hpp
 *  @brief Scoped logging verbosity for tests.
 */

#pragma once
#include "logging/logging.hpp"

// Exercise diagnostics without changing the caller's logging configuration.
class LogVerbosityScope {
    unsigned _previous;

  public:
    explicit LogVerbosityScope(unsigned verbosity)
        /** @brief Perform the specified evolution operation. */
        : _previous(Ariadne::Logger::instance().configuration().verbosity()) {
        Ariadne::Logger::instance().configuration().set_verbosity(verbosity);
    }
    ~LogVerbosityScope() { Ariadne::Logger::instance().configuration().set_verbosity(_previous); }
};
