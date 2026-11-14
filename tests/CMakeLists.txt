include(CTest)
option(EVOLUTION_CTEST_CONCISE "Use --concise for Evolution CTest executables (CI only)" OFF)
set(UNIT_TESTS
    test_orbit_collections
    test_scalings
    test_hybrid_time
    test_hybrid_enclosure
    test_hybrid_evolver
    test_hybrid_evolution_backend
    test_hybrid_evolution_tree
    test_hybrid_step_contracts
    test_hybrid_simulator
    test_first_order_pde
    test_first_order_pde_evolution_header
    test_flowtube_header
    test_orbit_evolution_header
    test_differential_inclusion_evolver
    test_enclosure
    test_enclosure_augmented_paving
    test_enclosure_operations
    test_enclosure_recombination
    test_enclosure_precision
    test_iterated_map_evolver
    test_vector_field_evolver
    test_vector_field_evolver_header
    test_vector_field_simulator
    test_evolution_dependency_boundary
)
foreach(TEST IN LISTS UNIT_TESTS)
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} PRIVATE ariadne-evolution evolution ariadne-testing)
    if(EVOLUTION_CTEST_CONCISE)
        add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}> --concise)
    else()
        add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    endif()
    set_property(TEST ${TEST} PROPERTY LABELS "evolution")
endforeach()
add_custom_target(tests)
add_dependencies(tests ${UNIT_TESTS})
