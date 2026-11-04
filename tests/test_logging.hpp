#pragma once
#include "logging/logging.hpp"

// Execute diagnostic formatting without filling CTest output. Redirect once,
// before the test creates workers; assertion output on stdout/stderr stays visible.
inline void redirect_test_logging_to_null() {
#ifdef _WIN32
    Ariadne::Logger::instance().redirect_to_file("NUL");
#else
    Ariadne::Logger::instance().redirect_to_file("/dev/null");
#endif
}
