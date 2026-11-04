# Coverage verification notes

The coverage target measures the public Evolution headers and implementation
files using the normal test suite. Tests assert numeric values, bounds,
parameter dimensions, event histories, stopping conditions and public contract
failures. The watertank workload is a separate benchmark.

## Enclosure representation invariants

The former private `Enclosure::_check` repeatedly tested map dimensions after
operations that already establish them. Public constructors validate incoming
argument dimensions; map application also checks result compatibility with the
auxiliary mapping. Parameterised steps check both map dimensions and domains
before mutation. Constraint entry points check arity. Function factories create
maps on the supplied domain; composition and restriction preserve the parameter
arity. Parameter extension embeds the time, dwell and constraint maps while
appending the state coordinate and its variable kind together. Products and
reconditioning rebuild these corresponding members together.

These facts make the old post-operation dimension-failure branches unreachable
for supported function factories and successful public operations. They have
been removed. The public legacy `HybridEnclosure::_check` symbol remains as a
compatibility hook; the contained labelled enclosure enforces its contracts at
creation and update boundaries.

The reduced domain starts equal to the parameter domain. Restriction intersects
it with a subdomain, extension takes the same Cartesian product for both domains,
and reconditioning constructs corresponding new domains before installing the
new representation. The pinned hull contractor modifies a coordinate only via
`max(old_lower,candidate_lower)` and `min(old_upper,candidate_upper)`
(`function/procedure.cpp` and `procedure.tpl.hpp`). `FloatDP` implements these
operations using `std::max`/`std::min` with the old endpoint first. Thus a valid
endpoint remains valid even if an intermediate candidate is NaN, and contraction
preserves inclusion. The redundant NaN restoration and subset-failure checks
have been removed on this basis.

The local Taylor reconditioning helper receives sorted, distinct parameter
indices selected from `[0,argument_size)`. It is called only after establishing
that at least one parameter must be discarded. Its error index is either a state
component index below the number of error variables or that number itself (the
uniform-error sentinel). Assertions for contrary internal states were redundant.

Generic `Function` domains are `EuclideanDomain(argument_size)`; their dimension
checks therefore follow from the preceding arity preconditions. Function-patch
domains are bounded boxes and retain explicit domain compatibility checks.

## Template and macro interpretation

LLVM can emit constructor entry variants and multiple template instantiations
whose individual counters differ. The official function/line/branch summary is
the completion criterion; structured per-symbol data and macro expansions are
used to explain residuals. Region coverage is recorded but is not a target.

## Hybrid event and timing regressions

`test_hybrid_evolver` checks a validated flow against analytic endpoints and
checks event kind, location, history and event-time bounds. It includes convex,
concave and degenerate guards, upper/lower semantics, invariant and progress
blocking, and uncertain initial times with the exact relation x(t0)=t0.

A correlated near-miss test uses x'=y'=1 and y0=x0-0.251+epsilon, with
|epsilon|<=0.0001. The maxima of y-x^2 and y-x^2+0.0001 at x=1/2 are
strictly negative. The crossing operation receives a genuinely const candidate
set containing both events and returns no crossings. It previously erased the
current iterator through const_cast; it now leaves its input unchanged, and the
caller updates active candidates from the returned map after iteration.

Progress predicates have no reset function, just like invariants. Transition
formatting now tests whether an event activates a transition before formatting
its reset. The progress scenario is executed with full logging as a regression
for the former null-function dereference.

The unused local arithmetic overloads had no call sites, as confirmed by clean
compilation after removal. The old time-creep strategy and alternate final-time
path were guarded by local constants `TIME_CREEP=false` and `if(false)`;
neither had a public configuration path. Their removal preserves the selected
guard-creep strategy and the existing final-time behavior.

## Portable constants and quiet diagnostics

The `_x` suffix requires an exactly representable double. On platforms with a
wider `long double`, non-dyadic literals such as `0.01_x` violate that contract.
Expected boxes and error intervals deliberately use `ExactDouble(0.01)` to name
the rounded binary64 value, matching the former values on Apple ARM. Large PDE
inputs instead use exact decimal literals (`1e60_dec`, `1e308_dec`). The critical
time logging threshold uses `_pr`, consistently with the crossing-time threshold.

Each runtime test redirects the logger once, at entry, to `NUL` on Windows and
`/dev/null` elsewhere. This keeps verbose formatting and its coverage active
while leaving assertion diagnostics on stdout and stderr visible. Public-header
compile checks remain independent of the logging helper. Benchmarks continue to
use the caller's logging options.

## Guard splitting and vector iterator lifetime

`List<HybridEnclosure>` is backed by `std::vector`. The grazing-guard path can
append hitting and past components while iterating over the original sets.
Appending invalidates the cached end iterator even without reallocation; when
reallocation occurs it also invalidates element iterators and references. The
former loop therefore had undefined behavior, detected by GCC's checked iterators
in `test_grazing_and_terminating_event`. It now snapshots the original count and
reacquires each original set by index. All accesses through the current set
reference finish before either append; local hitting/past copies remain valid.
Newly appended sets are not processed again, as required by the original algorithm.
The same grazing test runs on every supported platform and checks both optional
outcomes, numeric endpoint bounds, event histories and final times.

## Exported factory ownership in enclosure tests

The enclosure operations test only needs a validated Taylor function factory;
it does not exercise implementation-specific factory construction. It now uses
`make_taylor_function_patch_factory(sweeper)` exported by Function, retaining the same
DP precision and threshold sweeper. This keeps concrete model/factory vtable
instantiation in the provider library rather than forcing an incidental copy in
the test translation unit. GCC coverage builds previously diagnosed hidden
virtual `_unchecked_compose` overloads while emitting that unnecessary copy.
Warning policies and the enclosure expectations are unchanged.

`FunctionPatch::operator()` considers Taylor model argument types during overload
resolution, even when the call evaluates an approximate point. GCC requires
those candidate types to be complete when checking convertibility. The test
therefore includes `taylor_model.hpp` as well as the factory declarations; the
factory itself is still constructed through the exported interface.

## Hybrid solver failures and step-size contracts

UnknownSolutionException derives from EquationSolverException; their former
crossing handlers performed the same conservative INCREASING fallback. A single
EquationSolverException handler covers both. With a zero proof budget, the
standard interval Newton solver exercises that fallback and the corresponding
CONCAVE fallback for critical times. Tests verify the target, reset values,
event history and event-time bounds under both upper and lower semantics.
An unexpected backend runtime error is rethrown with its diagnostic preserved;
the formerly following assignment after a guaranteed throw was unreachable.

IntegratorInterface explicitly allows the selected step to differ from its
suggestion, including a larger step. HybridEvolver now restricts the returned
validated flow to min(available_step,configured_maximum). Restricting the time
domain preserves validation and the configured limit. A conforming integrator
that chooses twice the suggestion checks both the cap and analytic orbit bounds.

Initial invariant filtering is checked for partially valid and wholly invalid
sets, including termination before any events. The public location/space
contracts are tested with invalid locations and variable lists. The diagnostic
summary now includes the synchronized intermediate-set count.

Asynchronous initial times exercise straddling and unwinding policies. A final
slice can preserve its time as an equality constraint without a point-valued
interval hull. The test verifies the numerical final-time equality, conservation
of x-t, and the state and clock at an admissible parameter point. Configuration
and enum constants use portable binary64 or decimal constructors.

## Dependent timing and progress from asynchronous clocks

The three dependent-evolution policies are exercised by a derived evolver with
actual state-, time- and parameter-dependent step lengths. Two intermediate
clocks distinguish the policies; the exact linear solution checks the common
endpoint. StepKind streaming now names SPACE_DEPENDENT_EVOLUTION_TIME and
TIME_DEPENDENT_EVOLUTION_TIME instead of reporting UNKNOWN for valid values.
A nonlinear concave guard tests a nonzero critical-time error and both jumping
and non-jumping trajectories. A partial urgent crossing without creep checks
reset values, event identity and the analytic event-time interval.

UBSan on Linux exposed a progress-reporting defect for valid negative initial
clocks: ProgressIndicator converted a negative percentage to unsigned. Hybrid
progress now measures elapsed time from the initial lower clock bound, clamps
it to the requested duration, and uses 1/1 for an already reached horizon. The
existing straddling-time regression is retained, including its negative initial
times and verbose diagnostics. A zero-horizon case also checks the otherwise
undefined 0/0 percentage. No evolution time or state is altered by this change.

The zero-horizon expectation exposed a separate boundary defect:
_process_jumped_set discarded clocks equal to the final time. It now rejects
strictly later clocks and retains the invariant-valid final-time slice at the
closed horizon. No additional flow or reset is taken from that slice.

## Grazing partitions must remember earlier crossings

For x'=1, y'=0 with urgent guard y-x*x>=0 and positive y, every trajectory
crosses at t=-x0-sqrt(y0). Close to tangency, that root can lie beyond the
15/16-reduced solver interval, so the standard solver legitimately chooses the
GRAZING fallback. The old early exit tested elapsed_time<=critical_time. At an
endpoint after the maximum this is false, even for a trajectory which crossed
and returned; the endpoint was consequently retained without a jump.

The early exit now tests 0<=critical_time instead: a future crossing is measured
from the starting state, not the current endpoint. The existing three-way
partition then excludes trajectories that already crossed. Before the change,
the new analytic regression failed on the event count and reset state. It checks
both the pre-crossing horizon and the post-crossing endpoint. The mixed grazing
test additionally checks that non-jumping trajectories have y<=0 (within numeric
tolerance), rather than accepting their full original y interval. This removes
an avoidable spurious continuous path under upper semantics; it does not imply
that the previous upper enclosure omitted a true reachable state.
