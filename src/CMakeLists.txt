add_library(EVOLUTION_SRC OBJECT
    enclosure.cpp
    orbit.cpp
    differential_inclusion_evolver.cpp
    iterated_map_evolver.cpp
    vector_field_evolver.cpp
    vector_field_simulator.cpp
    first_order_pde.cpp
    hybrid_orbit.cpp
    hybrid_enclosure.cpp
    hybrid_simulator.cpp
    hybrid_evolver.cpp
)
target_include_directories(EVOLUTION_SRC PRIVATE "${PROJECT_SOURCE_DIR}/include/evolution")
