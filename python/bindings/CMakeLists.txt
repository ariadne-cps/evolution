include_directories(${Python3_INCLUDE_DIRS})
get_target_property(PYTHON_COMMON_SOURCE_DIR pyariadne-python-common SOURCE_DIR)

add_library(pyariadne-evolution-submodule-obj OBJECT evolution_submodule.cpp)
target_include_directories(pyariadne-evolution-submodule-obj PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${PYTHON_COMMON_SOURCE_DIR}/include"
    $<TARGET_PROPERTY:pyariadne-integration,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(pyariadne-evolution-submodule-obj PRIVATE
    evolution pyariadne-python-common pybind11::module pyariadne-system pyariadne-integration)
set_target_properties(pyariadne-evolution-submodule-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(pyariadne-evolution INTERFACE)
target_sources(pyariadne-evolution INTERFACE
    $<TARGET_OBJECTS:pyariadne-evolution-submodule-obj>)
target_include_directories(pyariadne-evolution INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<TARGET_PROPERTY:pyariadne-integration,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(pyariadne-evolution INTERFACE
    pyariadne-equation_solving
    pyariadne-integration
    pyariadne-system)
ariadne_add_python_module(NAME evolution MODULE_SOURCE evolution_module.cpp)
