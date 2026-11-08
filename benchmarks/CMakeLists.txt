# Benchmarking registers the aggregate benchmarks target for standalone builds.
# Production Evolution does not link the Benchmarking framework.
add_executable(benchmark_watertank benchmark_watertank.cpp)
target_link_libraries(benchmark_watertank PRIVATE
    ariadne-evolution evolution ariadne-benchmarking threading)
add_dependencies(benchmarks benchmark_watertank)
