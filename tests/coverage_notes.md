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
