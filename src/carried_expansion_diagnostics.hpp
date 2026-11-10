#pragma once

#include "function/taylor_function.hpp"

namespace Ariadne::EvolutionDetail {

// Diagnostics consume the concrete Taylor representation produced by the
// preconditioned integrator. An explicit stream also permits other clients
// to retain the profile without redirecting process-wide stderr.
Void print_carried_expansion_snapshot(
    OutputStream& output, String const& stage, SizeType step,
    ValidatedVectorMultivariateTaylorFunctionModelDP const& mapping,
    Sweeper<FloatDP> const& sweeper);

}
