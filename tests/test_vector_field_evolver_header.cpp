#include "evolution/vector_field_evolver.hpp"

#include <type_traits>

static_assert(std::is_base_of_v<Ariadne::VectorFieldEvolver::Interface,
                               Ariadne::VectorFieldEvolver>);

int main() { }
