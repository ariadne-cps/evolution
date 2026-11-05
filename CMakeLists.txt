cmake_minimum_required(VERSION 3.15)
project(Evolution VERSION 1.0)
option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

if(NOT TARGET ariadne-build)
    find_package(Git REQUIRED)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_CURRENT_SOURCE_DIR}"
                    submodule update --init -- ss/bui RESULT_VARIABLE _result)
    if(_result)
        message(FATAL_ERROR "Cannot initialize build submodule")
    endif()
    include("${CMAKE_CURRENT_SOURCE_DIR}/ss/bui/cmake/ProjectOptions.cmake")
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    require_same_dependency_commit(Build CHAINS
        "build" "integration|build" "equation-solving|build"
        "system|build" "system|hybrid|build")
    require_same_dependency_commit(Function CHAINS
        "integration|function" "equation-solving|function"
        "system|hybrid|labelled|expression-set|symbolic|function")
endif()

setup_project()
include_directories("${PROJECT_SOURCE_DIR}/include")
add_subdirectory(ss)

add_library(evolution INTERFACE)
target_include_directories(evolution INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/ariadne-evolution>)
target_link_libraries(evolution INTERFACE system integration equation_solving)

add_subdirectory(src)
target_link_libraries(EVOLUTION_SRC PRIVATE evolution)

ariadne_register_public_headers(
    TARGET evolution
    DIRECTORIES "${PROJECT_SOURCE_DIR}/include/evolution"
    FILES "${PROJECT_SOURCE_DIR}/include/ariadne-evolution.hpp")

setup_project_library(
    TARGET ariadne-evolution
    OBJECTS EVOLUTION_SRC
    LINK_LIBRARIES evolution)

setup_standalone_project_tests(
    TARGET ariadne-evolution
    TEST_DIRECTORY tests
    TEST_TARGET tests
    SOURCES "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/src")

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    ariadne_install_dependency_bundle(
        TARGET evolution
        DESTINATION include/ariadne-evolution)
endif()

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()
