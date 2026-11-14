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
 *  @brief Preserve logging execution and coverage while discarding console output in CI.
 */

#pragma once

#include <string_view>

#include "ariadne-testing.hpp"
#include "logging/logging.hpp"

/** @brief Run an Evolution test suite with ordinary Testing CLI semantics.
 *  @details The --concise flag controls the Testing reporter and redirects
 *  Ariadne Logger output to the platform's null device. Unlike muting or
 *  setting verbosity to zero, this keeps high-level log expressions active
 *  for code coverage. Without --concise, no Logger setting is changed.
 *  CTest enables --concise only when EVOLUTION_CTEST_CONCISE is ON.
 */
template<class Argv>
int run_evolution_test_suite(const Ariadne::Testing::Suite& suite, int argc, Argv argv) {
    bool concise = false;
    for (int i = 1; i < argc; ++i) {
        if (argv[i] != nullptr && std::string_view(argv[i]) == "--concise")
            concise = true;
    }
    if (concise) {
#ifdef _WIN32
        Ariadne::Logger::instance().redirect_to_file("NUL");
#else
        Ariadne::Logger::instance().redirect_to_file("/dev/null");
#endif
        // Continue evaluating all ordinary logging expressions, including
        // levels otherwise suppressed by the default logger verbosity.
        Ariadne::Logger::instance().configuration().set_verbosity(100u);
    }
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
