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
`TreeHybridEvolutionBackendFactory` records trajectory histories and genealogy
backed by `DerivationTree`. Its single checkpoint scheduler supports parent,
immediate and location admission with independent recombination policies.
Both backends support sequential and concurrent execution through ThreadManager.
Each call creates independent run state. See the [backend and trajectory-tree
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

## Location-by-location hybrid evolution

To finish all currently admitted work in one location before recombining pending
discrete successors, select phased admission:

```cpp
evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory(
    LocationCheckpointAdmissionPolicy{}, AllCompatibleRecombinationPolicy{}));
auto orbit = evolver.orbit(initial, termination, Semantics::UPPER);
```

Splits and continuous continuations execute in the active location. Discrete
successors wait until both its ready queue and its in-flight worker count are
empty. Compatible waiting sets are then recombined across trajectory parents;
the next waiting location is admitted as one batch. Different destinations or
event histories remain separate, and self-loops enter a later batch. The policy
does not synchronise physical time or modify the genealogy to determine admission.
Omit the second argument (or pass `NoRecombinationPolicy{}`) to retain this
admission without numerical recombination. All Tree configurations use the same
constructor: admission, grouping, trigger, recombiner, context, dispatch and
acceptance. Defaults select parent admission without recombination. Sibling
restrictions are expressed with `HybridRecombinationScope::SIBLINGS` in the common
grouping policies.
The default trigger evaluates before admission, so location admission merges at
its barrier, parent admission before releasing eligible descendants, and immediate
admission after each publication. See the design document for configuration examples.

For large frontiers, spatial grouping limits the diameter of every group and
previews its state hull before numerical recombination. Coordinate scales divide
state coordinates when comparing centre distances; the radius remains unscaled.
For example, this configuration limits position differences to 12 units and
other coordinate differences to one unit, with at most eight sets per group:

```cpp
const Scalings scales(1.0_x, {{RealVariable("x"), 12.0_x}, {RealVariable("y"), 12.0_x}});
evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory(
    LocationCheckpointAdmissionPolicy{}, ProximityRecombinationPolicy(1.0, 8u, 12.0, scales),
    EveryNStepsRecombinationTriggerPolicy(128u), ReferenceEnclosureRecombiner{},
    EnclosureRecombinationContext(12.0_x)));
```

Waiting sets are considered every 128 completed steps and before admission.
Replacement seeds wait for evolution before being eligible again, avoiding
repeated accumulation of numerical parameters. The context also checks the
actual output radius; unsuccessful merges preserve the original work.

## Hybrid evolution progress

Flat and Tree measure temporal progress toward the requested final time across
ready, waiting and executing checkpoints. Each new branch inherits its elapsed
time prefix. With a 200-second horizon, splitting a set at time 120 or creating
permissive alternatives at that time preserves 60% progress. Recombination
replaces source credits and remaining-time estimates by those of its outputs.
Different branch times and clock uncertainty can still change the percentage.

Backend `progress()` exposes the credited and remaining time, cumulative numerical
work and successful completion separately. Recombination emits begin and outcome
messages at the normal structural logging level; admission logs show actual hybrid
checkpoint counts. See the design document for the accounting and diagnostics.

## Recombination benchmark

The `benchmark_recombination` executable compares the reference continuous
`Enclosure` recombiner on **real post-reset enclosures** generated by the hybrid
evolver. Two nonlinear two-location systems have source dynamics
`x'=1, y'=x^2`, target dynamics `x'=1+y^2, y'=1`, identity reset and guard
`x+y>=0`. The urgent-crossing case starts in `[-2,-1] x [-1/2,1/2]`;
the permissive case starts in `[1,2] x [-1/2,1/2]`. The generator
terminates when the transition fires, then selects adjacent checkpoints in
the middle of the jump-time order. This is **not** a benchmark of hybrid
trajectory-tree recombination; it measures the numerical enclosure operation.

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

### Distributed fixed-cardinality subsets and verified quality diagnostics

The C++ recombination benchmark can select a fixed number of existing
post-reset enclosures spread across the **whole** fixed time window,
including the first and last eligible enclosures. The input enclosures
retain their original joint state/time constraints; no interpolation
or synthetic enclosures are introduced.

```bash
./build-release/benchmarks/benchmark_recombination \
  --sweep 1 --steps 0.015625 \
  --time-window 0.625,0.75 \
  --distributed-counts 2,4,8 \
  --samples 5 --warmups 1 \
  --csv build-release/recombination-distributed.csv
```

The CSV distinguishes `window-distributed` (a subset) from
`time-window` (the entire collection of eligible enclosures) and
records `generated`, `eligible` and `selected` separately. If a
requested cardinality exceeds the eligible count, that configuration
fails explicitly instead of duplicating sources. With a fixed step and
window, this experiment separates selected cardinality from the
physical time span, **but does not measure coverage of the complete
window** when it omits eligible sources.

In addition to interval-width ratios, the CSV now records four
sufficient **outer-projection containment** checks (joint coordinate
and clock intervals, plus the nonlinear target-flow invariant, before
and after continuation). A value `1` certifies containment of the
corresponding computed validated outer intervals; `0` is
**inconclusive**, not a counterexample to actual constrained-set
inclusion. The separate `guard_universally_verified` flag certifies
that the *entire output outer approximation* satisfies the urgent
guard equality `x+y=0` or permissive guard inequality `x+y>=0`.
`guard_violation_upper` bounds the potential excursion of that
outer approximation from the guard (absolute magnitude for the urgent
case, negative-side excess for permissive). Both guard metrics may be
conservative because they do not solve all constraints.

None of these checks, individually or collectively, constitutes a
general proof that the recombined **constrained set** contains the
union of the source constrained sets. Such a proof requires verified
joint feasibility/inclusion reasoning, not interval-hull comparisons.

### Recombination-only cardinality scaling

The existing `benchmark_recombination` executable also supports
`--scaling-recombine 1`, which measures **only the reference recombination**
on increasing numbers of real post-reset enclosures. One source pool is
generated per scenario outside timing, and reused for all counts; no sources
are duplicated to manufacture a larger dataset.

```bash
./build-release/benchmarks/benchmark_recombination \
  --scaling-recombine 1 --step 0.015625 \
  --scaling-counts 2,4,8 --samples 5 --warmups 1 \
  --profile-recombine 1 --csv build-release/recombination-scaling.csv
```

For 1024 inputs, supply `--scaling-counts 2,4,8,16,32,64,128,256,512,1024`
and a sufficiently fine `--step` to produce at least 1024 **distinct**
checkpoints. The benchmark fails explicitly if the source pool is too small.
Generating that pool can itself be expensive, but is not included in timing.

With `--time-window START,END`, each cardinality selects time-distributed
sources spanning the same fixed physical window. Otherwise it uses adjacent
prefixes of the pool, so the reference stays fixed but time spread grows
with cardinality. The CSV contains source counts, mean/median/min/max
recombination nanoseconds, output parameter/constraint counts and optionally
the same eight internal phase means. Output/provenance checks and source
reduction happen outside timing; there is no continuation, so this mode
measures recombination latency, **not** the end-to-end benefit of merging.
Profiling adds overhead: compare profiled and unprofiled runs separately.
The original benchmark modes are unchanged, and the default CSV name for
scaling is `recombination-scaling.csv`. This mode is not in CTest.
### Default parallel source recombination

`ReferenceEnclosureRecombiner{}` now defaults to
`EnclosureRecombinationExecution::PARALLEL`; explicitly construct
`ReferenceEnclosureRecombiner(EnclosureRecombinationExecution::SEQUENTIAL)`
to retain the original serial strategy. The parallel path is still subject
to regression and performance validation: both source emptiness filtering
and independent source-image graph propagation use ThreadManager/TaskGroup,
in batches of at most 64 sources. For each batch, at most one task per
configured worker is enqueued, with dynamic assignment of source indices to
balance irregular contractor costs and reduce scheduling overhead.
The filtering workers contract private copies of the inputs; the caller
moves the survivors into its ordered input list and reconstructs discarded
indices and provenance in original input order. During image propagation the
reference model is read-only; each task owns its normalised model, graph,
contractor, partial bounds and profiling counters. The final hull is
accumulated in original input order after joining tasks. Final output
reduction remains sequential. At concurrency zero the full serial path
is used automatically. With fewer than four inputs, filtering is serial;
with fewer than four surviving sources, image propagation is serial.

The scaling benchmark defaults to `--parallel-recombine 1`; use `--parallel-recombine 0` to compare the explicit serial strategy.
For example, run separate processes with identical settings (including
dataset-generation step and instrumentation):

```bash
./build-release/benchmarks/benchmark_recombination \
  --scaling-recombine 1 --parallel-recombine 0 --concurrency 4 \
  --step 0.015625 --scaling-counts 8,16,32 --samples 5 --warmups 1 \
  --csv build-release/recombination-sequential.csv

./build-release/benchmarks/benchmark_recombination \
  --scaling-recombine 1 --parallel-recombine 1 --concurrency 4 \
  --step 0.015625 --scaling-counts 8,16,32 --samples 5 --warmups 1 \
  --csv build-release/recombination-parallel.csv
```

The aggregate CSV includes execution mode and configured concurrency.
Both the phase profiler and per-source profiler remain opt-in. In
parallel mode the aggregate `filter_ns` counter measures elapsed wall-clock
time for the full parallel filter, including scheduling and joins; the
normalisation and graph phase counters sum per-task durations and can
overlap across workers. Per-source timers also overlap, so they must not
be summed into an exact wall-clock duration. An equivalence regression
compares serial and parallel output models, provenance and discarded
inputs for concurrency 0, 1, 2 and 4; a second regression exercises
130 inputs spanning three filter batches with discarded and surviving
sources on both sides of the batch boundaries. This is not a substitute
for release-mode performance measurements or thread-sanitizer testing.

### Per-source recombination diagnostics

The cardinality benchmark can additionally collect diagnostics for **every
original source index** (including filtered inputs), while leaving the ordinary
reference algorithm and aggregate timing CSV unchanged:

```bash
./build-release/benchmarks/benchmark_recombination \
  --scaling-recombine 1 --step 0.015625 \
  --scaling-counts 16,20,24,28,32 --samples 5 --warmups 1 \
  --profile-recombine 1 --profile-sources 1 \
  --source-csv build-release/recombination-sources.csv \
  --csv build-release/recombination-scaling.csv
```

The source CSV has one row per source **per cardinality** with the original
input index and selected pool index, input parameter/constraint counts, graph
dimension/constraint/active-image counts, and per-source means of filtering,
normalisation, image construction/ranging, graph preparation and reduction.
The prepared constraint solver's hull/shaving/monotonicity counters are averaged
over measured rounds. These diagnostics are only collected when explicitly
requested. The per-source timings omit shared result assembly and final
contraction, so they must not be summed into an exact total. As with the
aggregate profiler, profiling adds overhead and must not be mistaken for an
unprofiled latency measurement. Use the source indices to identify difficult
enclosures before changing contractor algorithms.

### Opt-in profiling inside reference recombination

Use `--profile-recombine 1` to measure internal stages in addition to the
existing end-to-end and phase timings. For example:

```bash
./build-release/benchmarks/benchmark_recombination \
  --sweep 1 --steps 0.015625 \
  --time-window 0.625,0.75 --distributed-counts 2,4,8 \
  --samples 5 --warmups 1 --profile-recombine 1 \
  --csv build-release/recombination-profile.csv
```

The CSV includes eight `recombine_*_ns` mean durations: input filtering
(including emptiness checks), parameter normalisation, residual-image
construction, direct image bounding, constraint-graph setup, constraint-solver
reduction, output assembly, and final output reduction. When profiling is
disabled (default), these columns are blank and the reference phase timers do
not read a clock. The constraint solver independently collects call-local
propagation statistics. Each profiled call resets the phase counters.

The internal durations are exclusive but **not exhaustive**: hull accumulation,
control flow, provenance and profiling overhead are not assigned to a stage.
Consequently they need not sum to `recombination_mean_ns`. Compare profiled
and unprofiled end-to-end runs on the same input to estimate instrumentation
overhead before claiming a speedup. Time-distributed subsets omit the remaining
eligible sources, so they measure cardinality scaling, not full-window coverage.

The reference implementation reuses models whose domain already matches the
graph parameter domain. It also substitutes certified singleton image bounds
for fixed graph coordinates, retaining their predicates for propagation into
the other coordinates. This avoids identity compositions and reduces graph
dimension, including after resets that set dwell time to zero. Source constraints
and their bounds are preserved. Graph contraction uses the constraint solver's
prepared propagation API, compiling each scalar function/procedure once and
restricting shaving to its genuine variable dependencies. Hull propagation and
adaptive shaving reach a fixed point using the existing validated contractors.
Preparation and solver statistics remain local to each recombination call.

### Why reference recombination is conservative

The reference upper-enclosure recombiner has a **constructive inclusion
argument**, not just empirical overlap checks. Every feasible source
parameter is lifted into a common unit parameter box. Its joint
state/absolute-time/dwell-time map is the reference map plus a residual
that belongs to a validated outer box; widened reference predicates
admit the same lifted witness. Final constraint contraction must retain
all feasible parameters. See the
[constructive upper-coverage proof](doc/evolution_tree_design.md#constructive-upper-coverage-proof-and-numerical-obligations)
for the exact assumptions, treatment of degenerate and non-dyadic
parameter domains, and the distinction between upper and lower
semantics. The unit suite includes a specific regression for affine
normalisation with non-dyadic endpoints and padded parameter dimensions.

Projection containment and guard-width metrics in the benchmark are
**quality diagnostics**. Their inconclusive values are not evidence of
unsound recombination, nor do positive projection checks independently
establish whole constrained-set inclusion. The proof instead relies on
validated model arithmetic and feasible-point-preserving contractors.

### Test reporter versus Ariadne Logger

The two output channels are independent. CTest passes `--concise` to all
Evolution C++ test executables, as in `threading`, `integration`, `system`
and `constraint-solving`. The flag controls **only** the Ariadne Testing
reporter; running a test executable directly without it uses the full
reporter. No CI-specific CMake switch is needed.

The Ariadne Logger output is separately redirected to `/dev/null` on POSIX
or `NUL` on Windows for every test run, including `--full`. This matches
the logging redirect in the `constraint-solving` test suite. **Verbosity is
never changed by the test runner.** Logging expressions at enabled levels
still execute (including explicit high-verbosity fixture variants); only the
terminal output is discarded, preserving logging-code coverage.

```bash
cmake --build build-release --target test_hybrid_evolver
ctest --test-dir build-release -R '^test_hybrid_evolver$' --output-on-failure
./build-release/tests/test_hybrid_evolver --concise
./build-release/tests/test_hybrid_evolver --full
```

The three header-only compilation tests do not use the testing framework and
continue to ignore their arguments.
