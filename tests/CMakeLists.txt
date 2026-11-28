include(CTest)
set(UNIT_TESTS
    test_approximate_inclusions
    test_domain_bisection_sampling
    test_initial_point_sampling_contract
    test_approximate_constant_parameters
    test_hybrid_constant_parameters
    test_constant_parameter_maps
    test_constant_parameter_auxiliaries
    test_evolver_paradigms
    test_orbit_collections
    test_scalings
    test_hybrid_time
    test_global_time
    test_hybrid_enclosure
    test_enclosure_recorded_draw
    test_hybrid_continuous_time_evolver
    test_hybrid_evolution_backend
    test_hybrid_evolution_tree
    test_hybrid_step_contracts
    test_hybrid_discrete_time_step_planner
    test_hybrid_discrete_time_step_image
    test_hybrid_discrete_time_step_executor
    test_hybrid_discrete_time_predicate_evaluator
    test_hybrid_discrete_time_upper_scheduler
    test_hybrid_discrete_time_evolver
    test_hybrid_discrete_time_approximate_evolver
    test_hybrid_crossing
    test_hybrid_inclusions
    test_hybrid_inclusion_impacts
    test_hybrid_continuous_time_approximate_evolver
    test_first_order_pde
    test_first_order_pde_evolution_header
    test_flowtube_header
    test_orbit_evolution_header
    test_continuous_inclusions
    test_enclosure
    test_enclosure_augmented_paving
    test_enclosure_operations
    test_enclosure_recombination
    test_enclosure_precision
    test_discrete_time_evolver
    test_discrete_time_inclusion_evolver
    test_discrete_time_approximate_evolver
    test_continuous_time_evolver
    test_continuous_time_evolver_header
    test_continuous_time_approximate_evolver
    test_evolution_dependency_boundary
    test_ordinary_dynamics_contract
    test_continuous_integration
)
foreach(TEST IN LISTS UNIT_TESTS)
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} PRIVATE ariadne-evolution evolution ariadne-testing)
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}> --concise --logtonull)
    set_property(TEST ${TEST} PROPERTY LABELS "evolution")
endforeach()

add_custom_target(tests)
add_dependencies(tests ${UNIT_TESTS})
