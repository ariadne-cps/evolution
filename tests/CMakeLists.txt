include(CTest)
set(UNIT_TESTS
    test_hybrid_enclosure
    test_hybrid_evolution
    test_hybrid_evolver
    test_hybrid_simulator
    test_first_order_pde
    test_first_order_pde_evolution_header
    test_flowtube_header
    test_orbit_evolution_header
    test_differential_inclusion_evolver
    test_enclosure
    test_enclosure_precision
    test_iterated_map_evolver
    test_vector_field_evolver
    test_vector_field_simulator
)
foreach(TEST IN LISTS UNIT_TESTS)
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-evolution evolution)
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    set_property(TEST ${TEST} PROPERTY LABELS "evolution")
endforeach()
add_custom_target(tests)
add_dependencies(tests ${UNIT_TESTS})
add_test(NAME test_evolution_dependency_boundary
    COMMAND ${CMAKE_COMMAND}
        "-DEVOLUTION_SOURCE_DIR=${PROJECT_SOURCE_DIR}/include/evolution"
        "-DEVOLUTION_CMAKE_FILE=${PROJECT_SOURCE_DIR}/src/CMakeLists.txt"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/check_dependency_boundary.cmake")
