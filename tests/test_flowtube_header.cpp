#include "ariadne-evolution.hpp"

#include <type_traits>
#include <utility>

// The interval returned by time_range() may be temporary: endpoints must not
// escape as references to that temporary.
using Ariadne::FlowTube;
static_assert(!std::is_reference_v<decltype(std::declval<FlowTube const&>().initial_time())>);
static_assert(!std::is_reference_v<decltype(std::declval<FlowTube const&>().final_time())>);

int main() { return 0; }
