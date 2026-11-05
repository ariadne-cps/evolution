#pragma once
#include "logging/logging.hpp"

// Exercise diagnostics without changing the caller's logging configuration.
class LogVerbosityScope {
    unsigned _previous;
public:
    explicit LogVerbosityScope(unsigned verbosity)
        : _previous(Ariadne::Logger::instance().configuration().verbosity()) {
        Ariadne::Logger::instance().configuration().set_verbosity(verbosity);
    }
    ~LogVerbosityScope() { Ariadne::Logger::instance().configuration().set_verbosity(_previous); }
};
