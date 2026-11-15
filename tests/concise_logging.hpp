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

/** @brief Discard logger output while preserving log expression evaluation.
 *  @details The nonblocking scheduler can leave a large diagnostic queue to
 *  drain when the test runner has already printed its final summary. The
 *  blocking scheduler handles messages when they are emitted, avoiding
 *  lengthy asynchronous shutdown. No verbosity or mute level is changed.
 *  The console reporter's --concise option is independent of this policy.
 */
class ScopedEvolutionTestLogging {
  public:
    ScopedEvolutionTestLogging() {
        auto& logger = Ariadne::Logger::instance();
        // The test has not yet started; switching here cannot interfere with
        // a running diagnostic fixture. It drains any startup queue.
        logger.use_blocking_scheduler();
#ifdef _WIN32
        logger.redirect_to_file("NUL");
#else
        logger.redirect_to_file("/dev/null");
#endif
        // The sink discards text: terminal-dependent line wrapping adds
        // work without observable output. LOGGING_* expressions still run.
        logger.configuration().set_handles_multiline_output(false);
    }

    ~ScopedEvolutionTestLogging() {
        // std::clog must no longer reference Logger's owned file buffer
        // when the test process reaches static destruction.
        Ariadne::Logger::instance().redirect_to_console();
    }

    ScopedEvolutionTestLogging(const ScopedEvolutionTestLogging&) = delete;
    ScopedEvolutionTestLogging& operator=(const ScopedEvolutionTestLogging&) = delete;
};

/** @brief Execute the complete suite with synchronous discarded log output. */
template<class Argv>
int run_evolution_test_suite(const Ariadne::Testing::Suite& suite, int argc, Argv argv) {
    ScopedEvolutionTestLogging logging;
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
