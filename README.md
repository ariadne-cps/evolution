# Ariadne Evolution

[![macOS Status](https://github.com/ariadne-cps/evolution/actions/workflows/macos.yml/badge.svg)](https://github.com/ariadne-cps/evolution/actions/workflows/macos.yml)
[![Ubuntu Status](https://github.com/ariadne-cps/evolution/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/ariadne-cps/evolution/actions/workflows/ubuntu.yml)
[![Windows Status](https://github.com/ariadne-cps/evolution/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/evolution/actions/workflows/win.yml)
[![Sanitizers Status](https://github.com/ariadne-cps/evolution/actions/workflows/sanitizers.yml/badge.svg)](https://github.com/ariadne-cps/evolution/actions/workflows/sanitizers.yml)
[![Coverage Status](https://github.com/ariadne-cps/evolution/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/evolution/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/ariadne-cps/evolution/branch/main/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/evolution)

Ariadne Evolution is the standalone C++20 numerical evolution layer for continuous, discrete, and hybrid dynamical systems: enclosures, orbits, simulators, evolvers, numerical propagation methods, and first-order PDE evolution.
Although independently maintained and usable, this repository is also part of the [Ariadne repository structure](https://github.com/ariadne-cps/ariadne/blob/main/doc/structure.md).

## Layout

Public headers: `include/evolution`; C++ implementations: `src`; C++ tests: `tests`; Python bindings and smoke tests: `python`; pinned dependencies: `ss`. The umbrella header is `include/ariadne-evolution.hpp`.

## Dependencies

Evolution uses the Git submodules [System](https://github.com/ariadne-cps/system), [Integration](https://github.com/ariadne-cps/integration), [Equation Solving](https://github.com/ariadne-cps/equation-solving), and [Build](https://github.com/ariadne-cps/build).

## Build and test

```bash
git clone --recurse-submodules https://github.com/ariadne-cps/evolution.git
cd evolution
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DPYTHON_ENFORCE_TESTS=ON
cmake --build build --parallel --target everything
ctest --test-dir build --output-on-failure
```

The standalone build provides `ariadne-evolution`; the exported module target is `evolution`.

## Python bindings

The `pyariadne-evolution` target composes the upstream Python bindings with the Evolution registrations. A standalone build creates a combined `pyariadne` extension and runs an Evolution smoke test with `ctest --test-dir build -L python --output-on-failure`. Python bindings are compiled across CI platforms; sanitizer jobs exclude Python import tests while retaining their compilation.

## License

Ariadne Evolution is distributed under the GNU General Public License v3.0. See [LICENSE](LICENSE).

## Watertank benchmark

The benchmark measures validated upper evolution of a tank with opening, open,
closing and closed valve modes. The dynamics are `height' = -height/50 +
3*aperture/10`; the valve takes four time units to open or close and switches
at heights 5.5 and 8. The initial valve is opening, with height in [0,1/16]
and aperture in [0,1/64]. The default stopping limits are time 80 or 8 events.

Build and run in Release mode for timing comparisons:

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=OFF
cmake --build build-release --parallel --target benchmark_watertank
./build-release/benchmarks/benchmark_watertank --samples 5 --warmups 1
./build-release/benchmarks/benchmark_watertank --help
```

Options include `--duration`, `--events`, `--step`, `--radius`, `--accuracy`,
`--samples`, `--warmups`, plus the standard logging and threading options
`--verbosity`, `--scheduler`, `--theme` and `--concurrency`. Each measured sample
reports elapsed milliseconds and reach/intermediate/final enclosure counts;
the summary reports the median, minimum and maximum. At verbosity 2 the final
locations, event histories, time ranges and state bounds are also printed.
Setup and final-result formatting are outside the timed region; integration,
orbit construction and any enabled internal logging are inside. Keep build,
configuration, concurrency and verbosity fixed when comparing results.

The benchmark is separate from CTest: correctness tests assert explicit
expectations and do not require generated plots or manual inspection.
