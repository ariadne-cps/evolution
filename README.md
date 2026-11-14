# Ariadne Evolution

[![Format Status](https://github.com/ariadne-cps/evolution/actions/workflows/format.yml/badge.svg)](https://github.com/ariadne-cps/evolution/actions/workflows/format.yml)
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

Hybrid evolution separates numerical step results from execution and Orbit
accumulation. The default `FlatHybridEvolutionBackendFactory` retains the current
workload scheduling; `HybridEvolverBase::set_evolution_backend_factory()` selects
another factory without changing the `orbit()` overloads or return type.
`TreeHybridEvolutionBackendFactory` selects serial parent-first processing with
trajectory histories and genealogy backed by `DerivationTree`. Each call creates
independent run state. See the [backend and trajectory-tree
design](doc/evolution_tree_design.md) for the current API and deferred operations.

## Dependencies

Evolution uses the Git submodules [System](https://github.com/ariadne-cps/system), [Integration](https://github.com/ariadne-cps/integration), [Equation Solving](https://github.com/ariadne-cps/equation-solving), [Tree](https://github.com/ariadne-cps/tree), and [Build](https://github.com/ariadne-cps/build). The C++ test suites additionally use [Ariadne Testing](https://github.com/ariadne-cps/testing), pinned as the test-only submodule `ss/tst`. Each registered test runs with a fresh fixture, retaining the original mathematical checks and diagnostic variants.

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

The `pyariadne-evolution` target composes the upstream Python bindings with the Evolution registrations. A standalone build creates a combined `pyariadne_evolution` extension and runs an Evolution smoke test with `ctest --test-dir build -L python --output-on-failure`. Python bindings are compiled across CI platforms; sanitizer jobs exclude Python import tests while retaining their compilation.

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

## Recombination benchmark

The `benchmark_recombination` executable compares the reference continuous
`Enclosure` recombiner on **real post-reset enclosures** generated by the hybrid
evolver. Two nonlinear two-location systems have source dynamics
`x'=1, y'=x^2`, target dynamics `x'=1+y^2, y'=1`, identity reset and guard
`x+y>=0`. The urgent-crossing case starts in `[-2,-1] x [-1/2,1/2]`;
the permissive case starts in `[1,2] x [-1/2,1/2]`. The generator
terminates when the transition fires, then selects adjacent checkpoints in
the middle of the jump-time order. This is **not** a benchmark of hybrid
trajectory-tree recombination, which is not yet implemented.

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=OFF
cmake --build build-release --parallel --target benchmark_recombination
./build-release/benchmarks/benchmark_recombination --samples 3 --warmups 1 --step 0.0625 --group-size 4
./build-release/benchmarks/benchmark_recombination --help
```

Run again with `--step 0.015625` and `--step 0.00390625` for the other two
baseline integration-step settings. `--accuracy` varies the integrator's
step error limit. The benchmark generates each dataset once, **outside timing**,
then measures fresh-copy constraint reduction, reference recombination, their
sum, independent continuation of each selected post-reset enclosure, and
continuation of the single recombined enclosure. Both continuations use the
same target-only dynamics, stopping at absolute time 4 (urgent) or 2
(permissive). Run with fixed concurrency and verbosity for comparisons.

The report includes source and output counts, provenance, final-continuation
counts, parameter/constraint counts, and outward-rounded interval widths of
state, absolute time and dwell time, before and after recombination. The widths
are **diagnostics of wrapping and numerical quality**, not a proof that the
recombined set covers the union. Certified coverage and clock/constraint
semantics belong to the existing C++ recombination tests. The reference
algorithm produces at most one enclosure and may inflate distant inputs;
this baseline deliberately selects nearby sources rather than merging a full
heterogeneous trajectory set. Benchmark runs are not registered in CTest.
