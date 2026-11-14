/***************************************************************************
 *            concise_logging.hpp
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

/** @file tests/concise_logging.hpp
 *  @brief Discard Ariadne Logger output without modifying verbosity or test reporting.
 */

#pragma once

#include "ariadne-testing.hpp"
#include "logging/logging.hpp"

/** @brief Redirect Ariadne Logger diagnostics independently of the test reporter.
 *  @details Logger output goes to the platform's null device on every test
 *  execution, just as in constraint-solving. Do not change logger verbosity:
 *  the same logging expressions are evaluated, including dedicated verbose
 *  fixture variants, so coverage is preserved. --concise and --full still
 *  control only Ariadne Testing's console reporter.
 */
template<class Argv>
int run_evolution_test_suite(const Ariadne::Testing::Suite& suite, int argc, Argv argv) {
#ifdef _WIN32
    Ariadne::Logger::instance().redirect_to_file("NUL");
#else
    Ariadne::Logger::instance().redirect_to_file("/dev/null");
#endif
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
