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
then measures reduction of fresh source copies **once per sample**. Both
continuations start from the same contracted continuous sets, wrapped as new
initial enclosures with empty event histories in the target-only automaton.
Both stop at absolute time 4 (urgent) or 2 (permissive). The end-to-end totals
are contraction + independent continuation, versus contraction + recombination
+ recombined continuation. Copying, projection evaluation and reporting are
outside timing. Keep concurrency and verbosity fixed.

Both the initial and the continued enclosure collections are compared by
outward-rounded coordinate projection widths for state, absolute time and dwell
time. The benchmark also uses **independent closed-form source and target flow
equations** to reconstruct initial coordinates from every generated post-reset
checkpoint and propagate its entire reduced parameter box to the final horizon.
This yields analytic interval projections for comparison with both evolutions.
An analytically impossible source/guard intersection or disjoint analytic and
computed outer projections fails the run. These are **necessary consistency
checks, not certified coverage of the source union by the recombined set**.
All checks and analytic projections run outside timing. The report includes
independent/recombined final counts, source provenance, constraint/parameter
counts, phase timings and end-to-end totals.
Widths close to roundoff at a fixed absolute time or dwell are reported
with an undefined (n/a) ratio. Each interval also reports both endpoints and
their outward-hull shortfalls so that near-boundary differences are visible;
these are **not** evidence that a source set escapes the output, since the
source hull may include infeasible points. In addition to coordinate projections,
the benchmark bounds the nonlinear target-flow invariant
`I = x - y - y^3/3`, both at recombination and after continuation. This
provides a correlation-sensitive numerical diagnostic, **not** a proof of
constrained-set inclusion. Certified numerical coverage and clock semantics
remain the responsibility of C++ correctness tests.

For scaling, run with --step 0.0625, 0.015625, 0.00390625 and --group-size
2, 4, 8 (if enough checkpoints are produced); vary --accuracy independently.
The generation orbit's reach/intermediate/final counts are printed, but they
are not an integrator accepted-step count. The benchmark is not in CTest.

### Reproducible C++ benchmark sweep

The benchmark executable runs the 3 x 3 step/group matrix **in process**,
without Python, subprocesses or parsing its own text output. Each scenario
produces one typed result, which is written directly to CSV. Failures carry
their original C++ exception message, and the sweep continues with other
configurations while returning a nonzero exit code if any failed.

```bash
./build-release/benchmarks/benchmark_recombination \
  --sweep 1 --samples 5 --warmups 1 \
  --csv build-release/recombination-sweep.csv
```

Default sweep axes are steps `0.0625,0.015625,0.00390625` and group
sizes `2,4,8`. Restrict or change them directly using
`--steps 0.0625,0.015625 --groups 2,4` (comma-separated).
The single-configuration mode remains unchanged: use `--step` and
`--group-size` without `--sweep 1` to obtain verbose numerical
diagnostics. Run a single configuration first when a finer step may be expensive.

The CSV includes configuration, selected checkpoint count, end-to-end means,
speedup, final-state width ratios, invariant inflation and analytic-reference
ratios. All ratios are diagnostics of outward bounds, not inclusion proofs.
Failed rows contain an error message; incomplete runs are never reported as
successful timing measurements. No benchmark run is registered in CTest.

### Fixed birth-time window (independent resolution study)

For the separate resolution experiment, keep the physical jump-time window
fixed while varying the integration step, e.g.:

```bash
./build-release/benchmarks/benchmark_recombination \
  --sweep 1 --steps 0.0625,0.015625,0.00390625 \
  --time-window 0.625,0.75 --samples 5 --warmups 1 \
  --csv build-release/recombination-fixed-window.csv
```

This selects all actual post-reset checkpoints whose validated time ranges
intersect the window, then imposes the window **on each source's joint
state/time function** and contracts its constraints. It does not synthesize
a boxed input. Source counts may differ between steps. In this mode
`--groups` and `--group-size` do not control selection, and the
CSV records the explicit selection mode and window endpoints (blank group
size). The older group-size sweep is unchanged and answers a different
question. A reported `ok` still means only the implemented
numerical consistency checks passed; it is not a general set-coverage proof.

### Phase breakdown and fixed-centre window-width sweep

The C++ sweep now writes mean contraction, recombination,
contraction-plus-recombination, independent continuation and recombined
continuation times separately, alongside both end-to-end means. Dataset
generation remains outside the timers. The CSV includes `window_width`
and reports phase times in nanoseconds.

To vary only the temporal window width at a fixed step and midpoint, run:

```bash
./build-release/benchmarks/benchmark_recombination \
  --sweep 1 --steps 0.015625 \
  --time-window 0.625,0.75 \
  --window-widths 0.0625,0.125,0.25 \
  --samples 5 --warmups 1 \
  --csv build-release/recombination-window-widths.csv
```

The centre is 0.6875 for every window. The selected checkpoint count is
an experimental result, not fixed across window widths. A separate sweep
with a fixed width and multiple steps isolates resolution. An `ok` result
means the numerical consistency checks passed, **not** a general
constrained-set inclusion certificate. The benchmark tests necessary
analytic overlap for every individual generated source as well as for
the union; both remain weaker than a mathematical inclusion proof.

In fixed-window mode, contraction for the selected input sets is timed.
The selection phase uses a separate feasibility copy to discard only
provably empty restrictions. The C++ correctness suite also tests a
whole-domain algebraic witness for two guard-constrained curves,
preserving state, absolute time and dwell time jointly.
