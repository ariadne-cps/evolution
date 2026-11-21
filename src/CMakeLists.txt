add_library(EVOLUTION_SRC OBJECT
    enclosure.cpp
    continuous_integration.cpp
    enclosure_recombiner.cpp
    enclosure_augmented_paving.cpp
    orbit.cpp
    discrete_time_evolver.cpp
    continuous_time_evolver.cpp
    vector_field_simulator.cpp
    first_order_pde.cpp
    hybrid_orbit.cpp
    hybrid_enclosure.cpp
    hybrid_simulator.cpp
    hybrid_evolution_backend.cpp
    hybrid_evolution_tree.cpp
    hybrid_continuous_time_evolver.cpp
)
target_include_directories(EVOLUTION_SRC PRIVATE "${PROJECT_SOURCE_DIR}/include/evolution")
