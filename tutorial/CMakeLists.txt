cmake_minimum_required(VERSION 3.15)
project(AriadneEvolutionTutorial LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build independently against an installed Evolution library, not build-tree
# headers or objects. CI checks this path in every Release configuration.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(AriadneEvolution REQUIRED)

add_executable(tutorial_evolution tutorial_evolution.cpp)
target_link_libraries(tutorial_evolution PRIVATE AriadneEvolution::ariadne-evolution)
