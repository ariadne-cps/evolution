/***************************************************************************
 *            test_logging.hpp
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

/** @file tests/test_logging.hpp
 *  @brief Silencing verbose diagnostics in evolution tests.
 */

#pragma once
#include "logging/logging.hpp"

/** @brief Redirect verbose diagnostic logging away from CTest output.
 *  Call before worker creation to avoid redirecting assertion output.
 */
inline void redirect_test_logging_to_null() {
#ifdef _WIN32
    Ariadne::Logger::instance().redirect_to_file("NUL");
#else
    Ariadne::Logger::instance().redirect_to_file("/dev/null");
#endif
}
