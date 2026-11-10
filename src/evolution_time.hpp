#pragma once

#include "numeric/dyadic.hpp"

namespace Ariadne::EvolutionDetail {

// Limit a step to a positive remaining horizon while preserving the
// integrator's requirement of exact binary64 representability.
Dyadic limit_step_to_remaining_time(Dyadic const& configured,Dyadic const& remaining);

}
