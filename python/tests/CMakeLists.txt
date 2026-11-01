add_test(NAME test_python_evolution_import
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_evolution.py")
set_property(TEST test_python_evolution_import PROPERTY
    ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:pyariadne>")
set_property(TEST test_python_evolution_import PROPERTY LABELS "python;evolution")
add_dependencies(tests pyariadne)
