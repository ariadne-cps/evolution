# System families: simulation and validated evolution

This is the implementation roadmap adopted in design v6.20 on 22 November
2026. Planned names in this document do not declare implemented APIs.

| Native family | Validated evolution | Numerical simulation | Current status |
| --- | --- | --- | --- |
| Discrete time (`DifferenceAlgebraicEquation` / `DifferenceAlgebraicInclusion`) | `DiscreteTimeEvolver<ValidatedTag>` for equations | `DiscreteTimeEvolver<ApproximateTag>` for equations and sampled inclusions | Validated equation UPPER/LOWER; approximate equation/inclusion with explicit tick-input policy |
| Continuous time (ODAE / ODAI) | `ContinuousTimeEvolver<ValidatedTag>` | `ContinuousTimeEvolver<ApproximateTag>` | Validated accepts both with explicit time; approximate accepts ODAE and structured ODAI with explicit time, fixed parameters and an explicit input policy |
| Hybrid discrete time | `HybridDiscreteTimeEvolver<ValidatedTag>` | `HybridDiscreteTimeEvolver<ApproximateTag>` | Validated deterministic UPPER public; bounded-input inclusion progression and state-only transitions to deterministic modes supported; inclusion targets/auxiliaries, LOWER and approximate pending |
| Hybrid continuous time | `HybridContinuousTimeEvolver<ValidatedTag>` | `HybridContinuousTimeEvolver<ApproximateTag>` | ODAE/ODAI timed UPPER evolution; approximate ODAE/structured ODAI supports external time, fixed parameters and explicit input policies; further contact handling planned |

`DifferenceAlgebraicEquation` keeps its existing name and contract. Hybrid discrete time now has its native model in `ariadne-cps/system`.
Its deterministic UPPER scheduler is now consumed by the public validated
hybrid discrete evolver. LOWER and approximate evolution are unsupported.
The tick semantics must define guard evaluation before/after the update,
priorities and simultaneous events, resets, multiple jumps at the same tick,
termination and the domain of integer indices. Continuous event semantics must
not be silently reused on a time grid.

## Shared responsibilities and separate numerical kernels

Share system ownership, initial-coordinate handling, labels and auxiliary
outputs, configuration, termination, trajectory collection, progress and
diagnostics, branch provenance and result production/publication where the
contracts agree. Reuse flat/tree scheduling and concurrency infrastructure for
both hybrid families when appropriate; do not duplicate entire hybrid engines.
Subdivision and reconditioning require a valid numerical representation.

Keep the step kernels explicit: discrete iteration applies a map/relation at
an integer index, whereas continuous integration produces a real-time reach
tube and endpoint. First-entry maps, Lie derivatives, within-step maxima and
creep belong to continuous localization, not discrete ticks. Generalization
must expose these capabilities rather than implying equivalent certificates.

Validated evolvers produce sets with their supported approximation semantics.
Evolvers produce approximate sampled trajectories. Sampling a trajectory
does not provide a validated enclosure, and API similarity cannot grant LOWER,
implicit DAE/DAI, multipatch or preconditioning support.

## Explicit continuous-to-discrete conversion

A discrete algorithm accepts a continuous system only through a converter
explicitly selected by its caller. No implicit ODAE/ODAI constructor or inverse
discrete-to-continuous conversion is planned.

A converter declares the step duration, time origin/grid, integration or
discretization scheme, parameter treatment and input policy. Approximate
simulation and validated step propagation have separate contracts. An ODAI
conversion must retain nondeterminism through a suitable typed discrete
relation/multimap, or reject the operation if `DifferenceAlgebraicEquation` cannot represent it.

Hybrid conversion additionally declares how within-step events, precedence,
resets, zero-time behavior and tangencies are handled. Endpoint guard samples
alone may miss crossings or contacts and cannot justify validated conservative
propagation. An explicitly approximate simulation conversion may use a weaker
contract, with its limitations visible to the caller.

## Implementation order and acceptance

1. Consolidate current hybrid continuous evolution, including ODAI UPPER
   IMPACT contacts with zero derivative, competing blockers, the closed horizon
   and parameter-dependent arrival. Keep unsupported LOWER cases explicit.
2. Complete continuous/hybrid time, fixed-parameter, input and auxiliary
   capabilities only for representations the numerical methods support.
3. Define and test the native hybrid discrete-time system in System.
4. Extract compatible hybrid infrastructure and implement the discrete hybrid
   evolver with declared UPPER/LOWER capabilities.
5. Add `DiscreteTimeEvolver<ApproximateTag>`; extend `ContinuousTimeEvolver<ApproximateTag>` with ODAI,
   explicit time, trajectory-constant parameters and reproducible input policies.
6. Implement `HybridDiscreteTimeEvolver<ApproximateTag>` using the same native tick model;
   extend `HybridContinuousTimeEvolver<ApproximateTag>` with the corresponding continuous
   capabilities and contact handling.
7. Implement the optional converters with separate simulation/validated
   contracts, testing within-step events and rejection of lost guarantees.
8. Verify the four-family/two-operation matrix: compatible and incompatible
   systems, integer/real horizons, reproducible seeds, inputs and parameters,
   simultaneous events, branch/reset behavior and numerical guarantees.
9. As the final Evolution activity, add missing Python bindings and give the
   Python tutorial functional parity with C++ for every available capability.
   Update installation, exports and documentation without compatibility aliases.

IMPACT, external global time and trajectory-constant physical parameters advance
the first two items. Further approximate evolver capabilities remain open. The discrete hybrid System model and a public deterministic UPPER validated
hybrid discrete evolver are implemented. Inclusion images, LOWER,
approximate hybrid discrete evolution and converters are still planned. The detailed
[continuous hybrid handler design](hybrid_transition_handlers_design.md)
records the implemented event capabilities and their mathematical boundaries.

Ordinary validated equations also support [auxiliary outputs depending on fixed
parameters](constant_parameter_auxiliaries_design.md), with explicit argument
layouts and rejection of conversions that cannot preserve those bindings.
Named [physical parameter maps](constant_parameter_maps_design.md) now support
recombination and compatible ordinary continuation without resampling.
[Hybrid constant parameters](hybrid_constant_parameters_design.md) extend the
same binding to ODAE/ODAI UPPER modes, guards, invariants, resets and equation
auxiliaries, including parameter-free intermediate modes. Partial parameter
tuples, LOWER parameterized evolution and preconditioning remain unsupported. Paving
decoupling is developed on the separate `decoupling` branch.


## Explicit parameter values for approximate ODAE evolution

[Approximate constant parameters](approximate_constant_parameters_design.md)
are selected by a RealValuation supplied to either continuous approximate
constructor. Named values remain fixed across trajectories, repeated calls,
clones and hybrid mode changes. Domains and required names are checked;
RK4 retains external stage time and the physical state dimension. Hybrid
modes can use different subsets of the selected names. Structured ODAI now also support [explicit input sampling](approximate_inclusions_design.md).
Random/grid parameter policies and new event localization remain separate work.

## Implemented paradigm axis

The explicit paradigm migration is implemented in the existing continuous,
hybrid continuous and validated discrete algorithms. See
[evolver paradigm contracts](evolver_paradigms_design.md) for actual input/result
interfaces, supported template specializations, Python registrations and tests.
This changes the naming of the planned algorithms above; it does not implement
missing discrete models, sampling policies or converters. Python tutorial
parity remains the final Evolution activity.

## Explicit continuous hybrid automata (v6.29)

The System model names are now `ContinuousTimeHybridAutomaton` and
`CompositeContinuousTimeHybridAutomaton`. Their shared continuous-time API is
`ContinuousTimeHybridAutomatonInterface`, derived from the time-independent
`HybridAutomatonInterface` structural contract. The common base provides only
name, mode existence, events, transitions, targets and state dimensions.

Continuous hybrid evolvers, factories and traits explicitly consume the
continuous-time interface, never the structural base. Existing C++/Python
consumer code and tests are migrated without aliases. This does not implement
any discrete automaton, nor does it modify continuous flow, guard, reset or
numerical semantics. The matching System gitlink must be updated jointly.

Next, extract shared graph/transition validation and composition code beyond
virtual interfaces before introducing `DiscreteTimeHybridAutomaton` and
its time-specialized interface. Its tick order and event semantics require a
separate, testable definition. Keep `decoupling` paving independent.

## Integer-time hybrid tick scheduler — first Evolution slice (v6.39)

The native System contract fixes a discrete step `n -> n+1`. Guards and
progression invariants are classified on the same pre-step source state.
Every possible event yields an independent reset successor in its target
mode; progression is permitted only where **all** event invariants hold.
Both reset and progression consume one tick. There are no zero-time chains,
implicit event priorities or compositions of simultaneous resets.

`HybridDiscreteTimeStepPlanner` provides a numerical-kernel-independent
structural scheduler: its evidence callback supplies certified
`FALSE`, `TRUE` or `UNKNOWN` for each source-set event predicate.
For urgent events, the invariant evidence is derived as the logical
complement of the guard evidence (and **not** independently tested).
For permissive events, an explicit invariant is evaluated. Each possible
branch is returned separately; uncertainty never silently discards a
branch. The `certain` field certifies only universal satisfaction of
that branch's **precondition** over the source set, and is explicitly
**not** a validated LOWER reachability or image certificate.

The planner returns source and successor integer indices and separate
transition/progression branch metadata, but **does not** compute set-valued
images, evaluate symbolic predicates, choose inclusion inputs, integrate
flows, or construct a numerical Orbit. The current continuous hybrid backend
uses continuous HybridEnclosure and HybridTime and is intentionally not
reused for discrete timing. Future validated and approximate hybrid
discrete evolvers must supply independently justified evaluation and
image kernels and reuse only semantics-compatible work scheduling and
publication components. Direct free-input references in event expressions
remain unsupported by the current model.

Tests verify empty event sets, all urgent guard truth values, overlapping
permissive advancement and reset, multiple simultaneous events, uncertain
progression, blocked progression and missing source locations. All existing
continuous evolution kernels are unchanged. Format remains deferred.

## Reused deterministic discrete hybrid image kernel (v6.40)

`HybridDiscreteTimeStepImage::apply` is the first numerical step adapter for
one branch returned by the existing integer-time planner. It uses
`Enclosure::apply_discrete_time_map_step` for deterministic progression,
`HybridEnclosure::apply_reset` for event-labelled transitions, and
System's `make_reset_function` for symbolic reset composition. Transition
resets are ordered by the *target state space* rather than the order of
assignments in user input, preserving the target's coordinate layout.

The adapter handles immutable, independently declared
`DifferenceAlgebraicEquation` source and destination modes. An event reset
consumes one tick and preserves its event provenance; the map is not applied
before a reset. Progression applies the equation map and advances the tick
through the existing enclosure primitive. The incoming source enclosure
remains unmodified; results are independent branches. The planner continues
to determine branch possibilities from pre-step evidence and the adapter
must not be mistaken for a numerical proof of those predicates. Applying
the branch image to the full source enclosure is suitable only as an
UPPER over-approximation until source restrictions are computed.

Inclusion source and inclusion destination modes are **explicitly rejected**
by this first adapter: their set-valued images and input-dependent auxiliary
correlations require dedicated rigor. No nominal input is sampled and no
claim of LOWER certification is made. No scheduling, work queues, numerical
integration engine, Orbit type or existing continuous evolver is changed.

Tests exercise progression/reset at exactly one tick, independent event
provenance and preservation of the source enclosure, coordinate order for
reversed reset assignments, and explicit rejection of unsupported inclusion
branches. The Evolution CI must validate this adapter before it is used by
a high-level discrete hybrid evolver. Format remains deferred.

## Deterministic hybrid discrete one-tick execution (v6.41)

`HybridDiscreteTimeStepExecutor::execute` composes the already existing
`HybridDiscreteTimeStepPlanner` with
`HybridDiscreteTimeStepImage::apply`. Its input is one source
`HybridEnclosure`, the source integer index and a certified predicate-evidence
callback. It returns a plan index pair and a vector of independent successor
enclosures paired with the original branch metadata. No single successor
is preferred, and every possible transition or progression branch is
computed from an unmodified copy of the same source set. An empty vector
indicates that the supplied evidence excludes every branch.

This is an *UPPER deterministic one-tick adapter*, not a full
`HybridDiscreteTimeEvolver<ValidatedTag>` implementation or a new
numerical engine. It preserves event provenance and relies on existing
enclosure map and reset operations to advance global time by one.
The caller must keep the supplied integer source index consistent with
the source enclosure time, and the callback must classify all predicates
over that enclosure, preserving correlation with declared auxiliary
outputs. A `certain` branch precondition is never a LOWER image certificate.
The executor intentionally does not intersect branch-specific predicates
before applying their images; unconditional source-set images are conservative
UPPER supersets. No lower semantics or reachability proof is inferred.

The executor rejects inclusion source modes and any unsupported target
image through existing kernel checks; it never chooses nominal inputs.
Failure of any branch aborts the entire execution rather than returning
an incomplete set of possible successors. It does not schedule repeated
ticks, subdivide sets, assemble an Orbit or change continuous-time evolution.

The planner now also validates every returned three-valued predicate
evidence, including permissive guard and invariant callbacks, instead of
silently interpreting an invalid enum as UNKNOWN or TRUE. Focused tests
cover overlapping permissive successors, unknown evidence, urgent
transition/progression alternatives, fully blocked steps, inclusion
rejection, event provenance and source-set immutability. The existing
regression matrix must validate the new code. Format remains deferred
until completion of the full redesign.

## Certified deterministic discrete hybrid predicate evidence (v6.42)

`HybridDiscreteTimePredicateEvaluator::classify` evaluates a System
`ContinuousPredicate` over the *entire* source `HybridEnclosure`.
For a deterministic difference equation mode, it constructs a scalar
sign-indicator function using System's existing
`make_constraint_function` on the physical state space and
`algebraic_sort` of the mode's explicit auxiliary assignments.
This substitutes dependent auxiliary expressions *before* interval
evaluation, avoiding the unsound independence of an auxiliary
output and its state inputs. The indicator is then evaluated using
`HybridEnclosure::range_of`, retaining the source parameterization
and constrained domain.

Evidence is TRUE only if the validated indicator lower bound is strictly
positive; FALSE only if its upper bound is strictly negative. All remaining
ranges, including those touching zero, yield UNKNOWN. This deliberately
avoids falsely certifying strict predicate truth at a boundary; it can
produce additional conservative UPPER branches. It does not construct
branch-specific restricted source enclosures, refine uncertainty by
subdivision, certify LOWER reachability or handle bounded-input
inclusion input/auxiliary dependencies.

`HybridDiscreteTimeStepExecutor::execute_certified` couples that evidence
to the existing planner/executor, so clients no longer need to supply
independent ad-hoc predicate sampling to obtain a deterministic UPPER
single-tick fan-out. The original callback interface remains available
for separately certified evidence. The caller must still keep the
explicit integer index consistent with the `HybridEnclosure` time map;
the adapter does not yet provide multi-step scheduling or an Orbit.
Inclusion modes and incompatible physical state spaces are rejected.
No numerical kernels in continuous evolution are changed.

Tests cover uniformly true and false predicates, uncertain boundaries,
symbolically correlated auxiliary assignments, urgent branching,
time advancement, immutability, inclusion rejection and state-space
mismatches. Functional CI is required to validate this new component.
Format remains deferred as requested.

## Bounded deterministic hybrid UPPER exploration (v6.43)

`HybridDiscreteTimeUpperScheduler::evolve` composes the certified
single-tick executor into a finite-horizon breadth-first exploration.
Every possible successor is retained as its own `HybridEnclosure` at
the next integer index. The result
`HybridDiscreteTimeUpperExploration` is a tree of nodes containing
the sample index, enclosure, stable parent index, and incoming branch
metadata (event, target and predicate certainty). Leaves are explicitly
classified as `HORIZON` or `BLOCKED`: an eventless progression is not
a blocked state, while a state admitting neither progression nor jump
is retained as a terminal at its original sample index.

The initial enclosure is never modified. Scheduling is breadth-first,
with independent source copies for all branches and no implicit
priority, reset chaining or zero-time event. The requested terminal
index is inclusive. A caller-specified positive `maximum_nodes` is a
hard limit: exceeding it throws, rather than discarding potentially
reachable branches and silently returning an unsound UPPER result.
No automatic merging or deduplication is performed, so combinatorial
growth remains explicit.

The caller supplies an integer `initial_index` that must correspond to
the initial enclosure's global-time map. The scheduler does not yet
enforce that consistency numerically, nor does it implement
subdivision, constrained branch images or robust inclusion modes.
This is *not* yet the public
`HybridDiscreteTimeEvolver<ValidatedTag>`, and no LOWER
semantics are claimed. The existing certified guard evaluator
and numerical image kernel preserve conservative UPPER branches.
Unsupported inclusion source modes throw instead of being sampled.
Test coverage includes linear progression, inclusive horizon, stable
parent links, permissive branching, event provenance, blocked leaves,
hard budget rejection, invalid horizons and inclusion rejection.

The crossing analysis files are intentionally kept separate:
`hybrid_crossing` exposes the public typed numerical certificate
contract, whereas `hybrid_crossing_analysis` is an internal
continuous-time numerical analysis algorithm based on Lie derivatives
and flow-localized roots. Neither is reused for integer-time ticks.

All changes are confined to Evolution on `multifunction`.
Format is deliberately deferred until the full redesign.
Functional CI must validate the new scheduler before advancing.

## v6.53 — Ordinary discrete LOWER witness audit and semantic terminology

Ariadne's terms are **LOWER** (every published enclosure contains a
reachable witness) and **UPPER** (each represented image encloses the
appropriate possibilities). They are deliberately distinct from
**INNER** (all represented points are reachable) and **OUTER**
(exhaustive union of UPPER images covers every reachable point).
For bounded LOWER enclosures, any represented point is within the
enclosure's diameter of a reachable witness; no separate distance
evaluation is needed.

For deterministic `DifferenceAlgebraicEquation` evolution from a known nonempty
initial enclosure, the sound map image propagates a witness.
UPPER-only subdivision and optional LOWER early termination preserve
the witness contract. A focused regression verifies affine and
nonlinear examples, prevents LOWER subdivision even when upper
subdivisions are enabled, and preserves existing stopping-policy
tests. See [the audit](discrete_time_lower_witness_audit.md).
Arbitrary constrained sources with unknown feasibility are **not**
yet proved to contain a witness; avoid calling the entire general
LOWER API fully audited without this qualification.

No hybrid or inclusion algorithm is modified by this audit;
the continuous and discrete hybrid regression suites remain required.

## v6.55 — Remove the obsolete ordinary discrete IteratedMap API

System's public deterministic discrete model is now
`DifferenceAlgebraicEquation`, which directly owns symbolic updates,
auxiliary declarations, compiled update and auxiliary functions.
The removed `IteratedMap` type is not exported in the umbrella header,
the hybrid automaton interface, Python, or the C++ implementation.

`DiscreteTimeEvolver<ValidatedTag>` accepts
`DifferenceAlgebraicEquation` as its public `SystemType`,
including `EvolutionTraits`, the Python constructor, tutorials,
and all ordinary discrete tests. The Evolution `ss/sys` gitlink is
pinned to the matching System migration commit. The numerical map
step and the LOWER witness proof are unchanged: this is an API
migration, **not** an additional approximation capability.
System and Evolution CI on the corresponding commits, including
both hybrid test families, must pass before accepting the migration.

The four-notion matrix (LOWER/INNER/UPPER/OUTER) is unchanged.
Ordinary discrete equation LOWER remains audited for known nonempty
initial sources, with arbitrary constrained-source witness evidence
still open. Format remains intentionally outside acceptance.

## v6.56 — Paradigm test migration after removing IteratedMap

The v6.55 API migration removed the legacy public `IteratedMap` class
from System and from Evolution, but `test_evolver_paradigms.cpp`
still contained two negative `static_assert` checks naming that
obsolete type. This is a **test compilation regression**, not a
failure of validated numerical evolution or the LOWER witness contract.

The assertions now explicitly check that
`EvolutionTraits<ValidatedTag, DifferenceAlgebraicEquation>` exists
and selects `DiscreteTimeEvolver<ValidatedTag>`, while
`EvolutionTraits<ApproximateTag, DifferenceAlgebraicEquation>`
is unavailable and a continuous approximate evolver cannot be
constructed from a discrete difference equation. The previous
negative capability checks are not relaxed.

An audit of the current C++ source, public headers, test suites and
Python/tutorial entry points in Evolution identified no other
operational references to `IteratedMap` in the inspected files.
The change keeps System at gitlink `7df54dc`, avoiding unnecessary
submodule churn while the coordinated CI is still running.
System's functional checks and all Evolution suites (continuous,
discrete and hybrid) remain required; Format remains non-blocking.
No matrix capability changes are claimed by this compile-only fix.

## v6.57 — Ordinary discrete validated LOWER: reject uncertified initial witnesses

The previously documented open initial-feasibility condition is now
enforced at the `DiscreteTimeEvolver<ValidatedTag>::orbit` boundary.
An unsplit bounded parameter box without constraints is a valid
nonempty source. A constrained source is accepted for LOWER only
if validated evaluation at a candidate parameter point certifies
**all** constraints simultaneously. Uncertified sources throw even
at horizon zero. Failure of finite candidate search does **not**
imply actual emptiness; this is a sufficient, deliberately
incomplete certificate. UPPER is not restricted by this check.

The mathematical invariant remains that each published LOWER
enclosure contains at least one actually reachable trajectory point;
the enclosure diameter bounds the distance to that witness. LOWER
is not INNER, and UPPER is not automatically OUTER. For a
deterministic difference equation, the unsplit validated image
propagates the certificate.

**Capacity matrix delta:** ordinary discrete equation / Validated
LOWER stays **partial** because some nonempty constrained initial
sets lack a certificate under this finite-point strategy. All other
cells are unchanged. Regression tests for known feasible and
infeasible constraints are added; the continuous and hybrid
discrete suites must also remain green. Format stays excluded.

## v6.58 — LOWER admission parity with ordinary continuous time

After comparing the ordinary continuous and discrete `orbit`
entry points, the additional v6.57 finite-point witness
admission gate was removed from the discrete path. Both
`ValidatedTag` ordinary equation evolvers use the same
initial-source nonemptiness assumptions; their numerical
UPPER enclosures preserve an actual witness under LOWER as
long as a witness exists initially and each published
LOWER branch retains it. The discrete path explicitly
never subdivides under LOWER; the analytic witness
regressions stay in place. A new regression protects
constrained but feasible sources that the finite-point
search would have rejected.

**Capacity matrix:** discrete ordinary equation / Validated
LOWER is now marked *implemented*, on equal footing with
continuous ordinary equation / Validated LOWER, not a
stronger separately certified paradigm. This does not
turn LOWER into INNER or UPPER into OUTER.

The continuous `_append_initial_set` preprocessing has
an unconditional size-based split, including for LOWER.
This existing divergence from witness-preserving branching
has been flagged for a focused separate review and must
not be hidden by an extra admission gate on discrete
evolution. All existing continuous and hybrid suites
remain mandatory in CI.

## v6.59 — Correct lifetime in constrained discrete LOWER regression

The initial-state parity test introduced in v6.58 indexed an ephemeral
identity vector function when making a scalar coordinate. A CI failure
in the function library's `get(i)` precondition revealed this test
construction error. The identity vector function is now retained as a
local owner while its component is used to build two constraints.
The feasible source, analytical witness, orbit assertions, and LOWER
admission policy remain unchanged. The complete capability matrix is
unchanged; ordinary discrete deterministic LOWER has the same source
contract as ordinary continuous deterministic LOWER. This is a test
repair, not a new capability. Functional CI is required before validation.

## v6.60 — LOWER for event-free timed/parameterised hybrid continuous equations

The deterministic validated hybrid continuous-time evolver now accepts
`Semantics::LOWER` in equation modes with explicit global time,
fixed physical parameters, or both **when the source mode has no
events**. Previously the mode gate rejected all such requests even
at zero horizon. An event-free deterministic flow has no ambiguous
branch: for an initially reachable witness, the validated
`compute_clocked_equation_step` and `apply_integration_step`
enclose its trajectory for every accepted step. The existing outer
reach/endpoint representation therefore also satisfies Ariadne's
LOWER witness contract; no INNER approximation or distance
calculation is needed. The parameter tuple remains fixed and
correlated across steps.

**The limitation is not cosmetic.** Clocked/parameterised hybrid
equations with event guards still use `publish_outer_jump_candidates`,
which may emit a reset from an intersected UPPER tube even if no
reachable trajectory satisfies the guard. Such branches cannot
be published under LOWER. The new gate explicitly rejects
`LOWER` in these modes with any events until a certified
crossing/guard mechanism is available. Ordinary autonomous
guarded LOWER keeps its existing event-aware implementation;
inclusion LOWER remains unsupported. Preconditioned integration
restrictions for clocked/parameterised modes are unchanged.

Tests now cover zero horizon and finite LOWER flow with fixed
parameters (timed and autonomous, across flat/tree backends),
global-time-only LOWER without physical parameters, analytical
state/time witnesses, physical-parameter/auxiliary correlation,
no history of spurious events, and explicit rejection of
event-bearing timed and parameterised LOWER modes.
The complete continuous/hybrid-discrete/impact suites remain
mandatory. The capacity matrix's hybrid continuous equation
LOWER cell **remains partial** until guarded timed and
parameterised branches can be handled without spurious witnesses.
Format remains outside the acceptance gate. The v6.60 CI result
must be checked on its own commit.

## v6.61 — Certified no-crossing and immediate-reset LOWER paths with hybrid events

Previously, timed and constant-parameter hybrid equation modes with any
event were rejected under LOWER, despite the autonomous equation path
already supporting whole-source immediate guards. This iteration adds
a **witness-preserving guarded subset**.

Initial event processing uses the existing certified
`classify_initial_guard` and `permits_initial_reset`. A fully
activated guard may publish an immediate LOWER reset; an unresolved
blocking guard prevents flow continuation. Invariant conditions must
hold throughout the original initial source. During each clocked
equation integration step, every validated reach tube is checked
against all blocking conditions via `HybridEnclosure::range_of`,
preserving dependencies on state, time and fixed parameters:
INVARIANT requires upper bound <= 0; URGENT, PROGRESS and IMPACT
require a strictly negative guard upper bound. If a proof is
unavailable, LOWER terminates that branch without creating a
potentially spurious reset or restricted flow enclosure.
PERMISSIVE events can be ignored as a valid nondeterministic
choice; a possible UPPER candidate is never advertised as LOWER.

The new regression covers both flat and tree backends, with
time-explicit or parameterized dynamics: excluded URGENT guards
permit LOWER progress, full-source URGENT guards cause a genuine
immediate reset and event history, and unresolved initial
activation produces no fabricated LOWER final state. Existing
event-free clocked and parameterized tests remain unchanged, as
do the autonomous geometry-aware crossing tests. System's
DifferenceAlgebraicEquation API is unaffected.

This does **not** establish completeness of LOWER for all
continuous-time hybrid equations: a blocking guard that becomes
active during a step currently causes conservative termination
rather than a certified timed reset. Providing a first-entry
certificate in the enlarged state/time/parameter chart is the
remaining mathematical work. Thus the capacity matrix cell
"Hybrid continuous / equation / Validated LOWER" stays ◐, now
covering more event-bearing cases, pending full functional CI
on the exact v6.61 checkpoint. Format is non-blocking.

## v6.62 — Global-time regression aligned with witnessed LOWER for nonautonomous equations

The `test_global_time::nonautonomous_equations_and_auxiliaries_in_both_evolvers`
case retained a pre-v6.60 expectation that
`hybrid.orbit(initial, HybridTime(0,1), LOWER)` throws. This is obsolete:
event-free deterministic hybrid equations with explicit time now admit
LOWER, as do the corresponding ordinary continuous equations.

Replace the negative assertion with positive analytical witness
checks for both flat and tree evolution backends, at zero horizon
and at time `1/8`. For `x'=t, x(0)=0`, prove the computed
validated state interval contains `x(1/8)=1/128`, the
auxiliary `a=x+t` contains `17/128`, and the final time is
`1/8`. Also verify no spurious discrete events and nonempty
reach output. Keep all existing UPPER accuracy assertions.

This is a **stale test contract**, not a reason to undo the
LOWER extension, suppress the failing case, or weaken mathematical
coverage. Production numerical code and the capability matrix are
unchanged. Hybrid continuous equation LOWER stays partial for
unresolved timed/parameterized first-entry crossings, even though
the no-event nonautonomous case is now supported. The v6.61 CI
detected this test failure; v6.62 must pass the functional CI
on its own commit before acceptance. Format remains nonblocking.

## v6.63 — Certified first boundary witness for a unique timed/parameterised URGENT crossing

The clocked/parameterized equation LOWER path now tests a narrow
constructive first-crossing certificate before its fallback to
no-blocker propagation or conservative termination. It requires:

1. Exactly one validated reach tube for the step, preserving a
   correlated parameter/time chart (no unstructured box image).
2. A source URGENT guard strictly negative over the entire initial
   enclosure and strictly positive over the entire *hypothetical
   uninterrupted* endpoint enclosure.
3. Every competing blocking event excluded from the entire tube
   (INVARIANT <= 0; other URGENT/IMPACT/PROGRESS strictly negative).

Continuity of the deterministic flow then ensures that every real
source trajectory has a first boundary crossing within the step,
and no competing blocker invalidates it. The guard-equality slice
of the validated tube contains at least one genuinely reached
boundary point. This is a valid LOWER witness, even though the
slice remains an outer enclosure; there is no INNER calculation or
numerical distance evaluation. The reset is applied to that
witness-bearing slice, and the successor retains its event history
and correlated time/parameter functions. The code does not infer
a precise crossing time or enumerate merely possible UPPER jumps.

A failure to prove the strict sign change or exclude competing
blocking events falls back to the existing conservative stopping
policy. An optional PERMISSIVE guard remains optional and may be
ignored. The old autonomous geometry-aware crossing path and
all UPPER candidate generation remain unchanged.

New dedicated tests prove an actual URGENT crossing for a
time-dependent equation `x'=1+t` and a fixed-parameter equation
`x'=p`, for both flat/tree backends. The initial condition is
`x=0`, the certified guard is `x>=1/64`, and the reset is
`x'=x+1`: each final LOWER enclosure must contain the exact
state `1+1/64`, with the expected event history and binding.
A second negative regression uses competing URGENT boundaries
`x>=1/128` and `x>=1/64` and requires that LOWER never
invent a target final state from conflicting tube intersections.

**Matrix status:** hybrid continuous / equation / Validated LOWER
remains ◐ because cases without strict endpoint sign change,
tangent IMPACT, concurrent crossings and other complex guard
geometries still terminate conservatively. The additional
crossing case has been implemented but is not accepted until
all functional CI passes on this checkpoint. No unrelated matrix
cell or tested numerical semantics have been weakened.


## v6.64 — Joint explicit-clock/physical-parameter LOWER first-crossing regression

This checkpoint adds test coverage, not a new numerical algorithm or a
promotion of hybrid continuous equation LOWER from partial (◐) to complete.
The v6.63 first-boundary witness implementation is exercised when **both**
explicit global time and fixed physical parameters occur in the *same*
equation mode. The dynamics are `x'=p+t`, `p in [2,3]`, with the
time-dependent URGENT guard `x>=t+1/64`. A reset `x'=x-t+1`
must preserve a genuinely reached first-boundary witness at
`x=1+1/64` for both flat and tree evolution backends. The final
enclosure must retain the physical parameter binding, target location
and exact event history.

A second test case for each backend adds an earlier competing URGENT
boundary `x>=t+1/128` and demands that no LOWER final target be
published. These tests specifically exercise the state/time/parameter
argument layout and guard/reset correlation. The prior timed-only and
parameter-only positive/negative regressions remain unchanged, as do
all UPPER, inclusion, autonomous, and hybrid discrete tests.

The earlier v6.62 baseline passed Ubuntu, macOS, Windows, Coverage and
Sanitizers. The v6.63 checkpoint was still being validated when this
test-only increment was prepared; do not infer its acceptance from the
v6.62 CI or infer acceptance of this v6.64 commit before its own five
functional workflows complete. Format remains nonblocking.


## v6.65 — Physical-parameter and clock-dependent URGENT guard regression

This checkpoint adds one focused C++ regression and does **not** change
the validated numerical kernel, the capability matrix, or the UPPER path.
For the joint equation `x'=p+t`, with `p in [2,3]` and initially
`x=t=0`, the URGENT condition `x>=t+p/256` must certify a LOWER
first crossing. The reset `x'=x-t-p/256+1` maps that real boundary
to `x=1`; the final result must contain this witness, retain the
physical parameter binding, target mode, single event, and exact
horizon time. Both flat/tree evolution backends are exercised.

The same two backends also exercise a second, earlier URGENT condition
`x>=t+p/512`, with a distinct reset, to ensure that the first
boundary cannot be incorrectly selected as a certified unique
crossing when another blocking event is present. The expected LOWER
final collection is empty in that deliberately unresolved case.
All previous hybrid continuous/discrete, autonomous, UPPER and
inclusion regressions remain intact.

The v6.64 functional CI was still running when this test-only commit
was prepared. Acceptance of v6.65 requires separate successful
Ubuntu, macOS, Windows, Coverage and Sanitizers workflows on its own
SHA; no prior SHA is substituted for that gate. Format remains
nonblocking. Hybrid continuous equation LOWER stays partial (◐).

## v6.66 — Fix the v6.65 regression's event-limit versus horizon-time expectation

The new `test_lower_guard_depending_on_clock_and_parameter`
reported a failed check for `t=1/8` despite reaching the certified
URGENT reset state. This was an **incorrect test expectation**, not
evidence of a LOWER numerical-kernel defect. In
`_process_jumped_set`, `HybridTime(1/8, 1)` stops upon the first
completed event with `DISCRETE_LIMIT` and records the reset
enclosure at the **crossing time**, not at the global time horizon.

The corrected test preserves that event-limited orbit and its
witness/event/parameter assertions, and additionally proves that the
reset occurs strictly before the horizon. A separate orbit with
`HybridTime(1/8, 2)` allows the single jump followed by continuous
evolution in the stationary target mode, and checks the exact horizon
time, target, one-event history, parameter binding and reached state
`x=1`. The competing-URGENT negative cases remain unchanged. No
test has been removed; no algorithm, semantics or model parameters
have been weakened or modified. Hybrid continuous/equation LOWER
remains partial (◐).

v6.64 passed Ubuntu, macOS, Windows, Coverage and Sanitizers on
`c68deeec6cfe96d19fa862f1e38d1e02e5734a1c`. v6.65
(`115d7413c000480276511dc05d4bde25599c319b`) failed
Ubuntu/macOS Debug test execution on this incorrect time assertion.
Acceptance of this fix requires separate successful Ubuntu, macOS,
Windows, Coverage and Sanitizers workflows on its new commit SHA;
Format is nonblocking.

## v6.67 — URGENT LOWER crossing certified at a closed step endpoint

The clocked/parameterized equation LOWER first-crossing certificate
no longer requires the hypothetical uninterrupted endpoint guard to
be **strictly positive**. If the guard is certified strictly negative
on the complete initial source and **nonnegative** on the complete
validated endpoint, continuity of the deterministic flow guarantees
that at least one (indeed each source trajectory's) first zero occurs
within the closed step. Other blockers must still be certified inactive
over the entire reach tube, and the single structured validated tube
and guard-equality reset mechanism are unchanged. An interval endpoint
that merely *possibly* contains zero does **not** pass this gate.

The new regression `test_lower_urgent_at_closed_step_endpoint`
uses a global-time URGENT guard `t>=1/32` and zero state flow, with
and without trajectory-constant physical parameters and on both flat
and tree backends. At the exact integration-step endpoint, LOWER must
record the real reset witness `x=1` at `t=1/32` and its event and
parameter history. The same model with the shorter horizon `1/64`
must remain in the source mode without inventing an event.

Only the single endpoint-sign condition in the validated hybrid
equation kernel and its dedicated regression are changed. All
previous positive and negative tests remain active; UPPER, autonomous,
hybrid discrete, inclusion, and Python behavior are not altered.
Continuous hybrid equation LOWER remains partial (◐): tangential
IMPACT, competing crossing geometry, undecidable endpoint signs and
other complex first-entry cases remain conservative.

At publication, v6.66 Ubuntu and macOS had passed on
`106657f2f097afee82c3a2a18390ffc9df8a9cf1`, while Windows,
Coverage and Sanitizers were still being checked. This v6.67
checkpoint requires its **own** successful functional CI across all
five workflows before acceptance. Format remains nonblocking.

## v6.68 (stage 1) — Shared paving-independent sampling infrastructure

Copied `sampling_counts.hpp`, `domain_bisection_sampling.hpp/.cpp`, and
`test_domain_bisection_sampling.cpp` from evolution/decoupling byte-for-byte.
Adapted `initial_point_sampling.hpp/.cpp` solely to the renamed current
`ContinuousTimeEvolver<ApproximateTag>` and
`HybridContinuousTimeEvolver<ApproximateTag>` API types, retaining the
same midpoint-of-domain-leaves algorithm. Copied its separate sampling
contract regression and registered both suites in CMake. Unspecified
coordinate counts default to one; flat dimensions are never multiplied.
The ongoing migration will make the existing Approximate evolvers call
these helpers, eliminating their grid/paving sampling paths. The validated
algorithms and per-step `InputSamplingPolicy` are unaffected. Neither
this stage nor a subsequent commit is accepted without its own five-way
functional CI. Format remains nonblocking.

## v6.68 (stage 2) — Approximate initializers use decoupling's sampling core

Both `ContinuousTimeEvolver<ApproximateTag>` and
`HybridContinuousTimeEvolver<ApproximateTag>` now call the same
paving-independent `initial_point_sampling` adapters brought from
`decoupling`. The prior paving/grid-based sample generation, artificial
widening and associated configuration knobs are removed. Their named
coordinate counts use `SamplingCounts` identically: the default is
a single midpoint, explicitly requested per-axis counts multiply, and
degenerate dimensions contribute one unique sample. Explicit point
inputs are not resampled. Time integration, RK4 sampling for inclusion
inputs, guards, resets and other numerical semantics are unaffected.

Both existing Approximate regression suites now check cardinality and
domain membership, including a flat-coordinate case. The two
`decoupling` unit suites are retained unchanged. The C++ tutorial and
Python configuration binding use the new sampling-count API, removing
obsolete grid-oriented methods. This commit does not introduce
Approximate discrete evolution or an initial distribution over general
nonbox constrained sets; in particular, sampling centres of a bounding
box must not be represented as satisfying arbitrary nonlinear
constraints. This stage requires its own five functional workflow
passes before being considered accepted. Format remains nonblocking.

## v6.69 — Restore test-local includes for private box samplers

The paving-independent sampler migration imported `domain_bisection_sampling.hpp`
and `initial_point_sampling.hpp` into `src/`, but left out the original
`decoupling` CMake target-local include directories for the two imported
regressions. The resulting `test_domain_bisection_sampling.cpp` compile
failed before tests could run because its private header was not on the
test target's include path. Restore both target-local
`target_include_directories(... PRIVATE "${PROJECT_SOURCE_DIR}/src")`
directives exactly as in `decoupling`. No production source,
sampling policy, test expectation or public include path changes.
The v6.67 five-workflow baseline was fully passing; v6.68 was
not validated. Accept this change only after successful functional CI
on the new commit. Format remains nonblocking.

## v6.70 — Uniform test compilation and no Windows timeout override

The two temporary sampler-specific include directories introduced in v6.69
have been removed from `tests/CMakeLists.txt`. Both sampling headers were
moved from `src/` to the normal `include/evolution/` header hierarchy;
all source and test includes now use `evolution/...`. All tests compile
through the same CMake loop and the same default public include path, with
no target-specific exceptions. The sampling implementation and tests remain
unchanged, including the no-paving behavior and one-midpoint default.

The special `if(WIN32 AND CMAKE_BUILD_TYPE STREQUAL "Debug")` block that
overrode the CTest timeout to 600 seconds for two evolver suites has also
been removed. CTest's ordinary timeout rules now apply uniformly.

The previously green v6.67 is still the validated baseline. v6.68 and
v6.69 failed compilation; this correction is only accepted after the
five functional workflows pass on the exact new SHA. Format remains
nonblocking.

## v6.71 — Approximate deterministic discrete-time equations

Implements `DiscreteTimeEvolver<ApproximateTag>` for the native
`DifferenceAlgebraicEquation`. It exposes only point-trajectory
`orbit(points, integer_horizon)`, plus an overload taking an
initial box and sampling it once, via the existing
paving-independent `sample_vector_field_initial_points` helper.
The default is exactly one initial midpoint; optional
`SamplingCounts` chooses the number of deterministic bisection
subintervals on each named, nonflat axis. Explicit labelled points
are reordered by state variables and never resampled. The result
contains the initial sample and a point for every integer tick,
with output auxiliary variables recomputed at each point.
Empty initial lists, zero/negative horizons, missing labels,
incompatible box dimensions, and flat-coordinate sampling have
dedicated regression checks. The new specialization is exposed
through the public umbrella header and `EvolutionTraits` and
does not accept `Semantics::UPPER` or `Semantics::LOWER`.

This introduces no paving dependency, new numerical approximation
semantics, extra test-specific compiler options, or Windows
timeout override. Hybrid-discrete and inclusion Approximate
remain separate future increments. The implementation is
published as a candidate and requires all five functional CI
workflows on the exact commit before its capability is marked
validated. Format remains nonblocking.

## v6.72 — Unified evolver terminology and nonconst orbit test fixes

All occurrences of the former alternative terminology have been
removed from repository-controlled filenames, source text, test names,
documentation, tutorials and Python bindings. The three renamed
sampling artifacts are now `initial_point_sampling.hpp`,
`initial_point_sampling.cpp`, and
`test_initial_point_sampling_contract.cpp`. The build registrations
and every corresponding include directive follow the new paths.
The Python binding export helper is called
`export_approximate_evolver`, without altering exported modern
evolver type names or numerical behavior. The obsolete-name negative
assertions were removed from the Python suite because preserving their
old string literals would violate the strict terminology rule;
the existing positive paradigm export checks remain.

The imported hybrid approximate sampling regression was declaring
sampled orbit objects const even though the existing hybrid
`Orbit::curve` and `Orbit::curves` readers are nonconst methods.
The test now declares those orbit variables nonconst, retaining all
cardinality, domain-membership, flat-axis and endpoint assertions.
The newly introduced approximate discrete-time test likewise
retains nonconst orbit variables when calling the corresponding
accessors. No new per-test compiler flags, target include directories
or Windows timeouts have been introduced.

This is an API spelling/compilation repair only, not a numerical
capability promotion. Previous v6.68 through v6.71 commits did not
pass the complete five-workflow gate. Acceptance requires successful
Ubuntu, macOS, Windows, Coverage and Sanitizers on this exact commit.
Format stays nonblocking.

## v6.73 — Approximate discrete tick conversion compile fix

The new `DiscreteTimeEvolver<ApproximateTag>` used
`tick + FloatDP(1, dp)` and attempted to assign the result to a
`FloatDP`. Ariadne's raw `FloatDP` addition intentionally yields
`Bounds<FloatDP>`, so compilation fails under Clang/GCC with this
conversion. The discrete time label is now computed directly from the
exact integer loop index via `FloatDP(step + Integer(1), near, dp)`,
avoiding interval-to-scalar narrowing and cumulative rounding drift.
The previous equality guard still prevents unrepresentable consecutive
floating-point tick labels. State propagation, auxiliary recomputation,
initial sampling, orbit tests, CMake test registration, and the
unified evolver terminology remain untouched. This correction must
pass the five functional workflows on its own commit; it is not yet
considered validated.

## v6.74 — Restore distinct evolver variables in installed C++ tutorial

The repository-wide terminology rename replaced two previously distinct
local identifiers in `tutorial/tutorial_evolution.cpp` with the same
`evolver` variable name in a single `main()` scope. It consequently
made the installed tutorial invalid C++: the approximate continuous
evolver and validated discrete evolver must remain distinct types and
objects. Their local identifiers are now `approximate_evolver` and
`discrete_evolver`, with each configuration and orbit call using the
corresponding object. Model definitions, result assertions, time
budgets, and evolution algorithms are untouched. This restores the
original distinction rather than removing tutorial checks.

On v6.73 the installed C++ tutorial step failed for observed
Ubuntu Release jobs; other CI jobs also failed and require their
own diagnosis. This fix does not claim to resolve those independent
test-compile issues. Only five functional workflow successes on
this exact new SHA can validate the checkpoint.

## v6.75 — Fully descriptive tutorial evolver identifiers

The installed C++ tutorial now uses `approximate_continuous_evolver`
for `ContinuousTimeEvolver<ApproximateTag>` and
`validated_discrete_evolver` for `DiscreteTimeEvolver<ValidatedTag>`.
All corresponding configuration and orbit calls use the renamed
objects. This is an identifier-only change: no model, simulation,
validation, sampling, assertion or CMake behavior was changed.
Acceptance still requires passing the five functional CI workflows
on the new commit; Format remains nonblocking.

## v6.76 — Disambiguate approximate and validated hybrid-continuous tutorial evolvers

The installed C++ tutorial contained a second name collision introduced
by the earlier repository-wide identifier normalization:
`HybridContinuousTimeEvolver<ApproximateTag>` and
`HybridContinuousTimeEvolver<ValidatedTag>` were both declared in
`main()` under the same identifier. This causes a C++ redefinition
error and incorrectly resolves the subsequent validated-only
`set_integrator`, step-size, enclosure-radius, and UPPER orbit calls
against the approximate type. The two local objects are now named
`approximate_hybrid_continuous_evolver` and
`validated_hybrid_continuous_evolver` respectively, with every
call matched to its intended paradigm. This follows the explicit
`approximate_continuous_evolver` and `validated_discrete_evolver`
conventions established in v6.75. The tutorial checks, mathematical
models, numerical algorithms, test registrations and CMake rules
are unchanged. No special include path or timeout is added.

This corrects the specific installed tutorial compilation errors
reported at the v6.75 checkpoint. Other Debug, sanitizer and platform
errors, if any, require independent diagnosis. Acceptance is conditional
on passing Ubuntu, macOS, Windows, Coverage and Sanitizers for this
exact new commit. Format remains nonblocking.

## v6.77 — Uniform tutorial output descriptions

All user-facing output prefixes in the installed C++ tutorial are now
concise semantic descriptions rather than a mixture of C++ type names
and generic text. The examples print:
`Approximate continuous equation`,
`Validated continuous equation`,
`Validated continuous inclusion`,
`Validated discrete equation`,
`Validated hybrid continuous equation`, and
`Validated hybrid continuous inclusion`.
The two-mode hybrid automaton section prints
`Hybrid transition model`, because it checks its transition structure
and budget but does not compute an orbit. The final PDE example uses
`First-order PDE` instead of its class name.

The discrete example is an *equation*, not an inclusion; the labels
therefore reflect the concrete instantiated models, not an illustrative
placeholder. The orbit sizes, enclosure counts, trajectory counts,
model construction, assertions, integrators, semantics, and CMake
registration remain untouched. This is a display-only change; no new
algorithm or capability is claimed. CI validation is required on the
new SHA; Format remains nonblocking.

## v6.78 — Keep the installed C++ tutorial demonstrative, not test-like

Removed the three explicit correctness checks in
`tutorial/tutorial_evolution.cpp`: asserting exactly one approximate
continuous trajectory, checking the manually constructed hybrid
transition relation, and checking the hybrid jump budget. The
`stdexcept` include is no longer needed. The trajectory count is
printed directly from the produced orbit instead of being compared
against a hard-coded expectation. All semantic output labels from
v6.77 remain intact.

Also removed the duplicate explicit `pde.validate()` call:
`first_order_pde(pde, profile, double_precision)` already invokes
`pde.validate()` internally, so PDE input precondition validation is
not lost. All numerical demonstrations, configurations, orbit calls,
printed results, and actual model definitions remain otherwise
unchanged. The dedicated test targets remain responsible for
assertions; no tests, CMake policies or CI checks have been weakened.
This change is not considered validated before its five functional
GitHub workflow results are successful. Format is nonblocking.

## v6.79 — Repair distinct validated/approximate time-aware test objects

The Debug build of `test_ordinary_dynamics_contract.cpp` reported two
different `ContinuousTimeEvolver` specializations declared with the
same local identifier `timed_evolver` in one scope. The resulting
redeclaration made the approximate `step_size()` assertion resolve
against `ContinuousTimeEvolver<ValidatedTag>`, which has no such
configuration method. The local objects are renamed to
`validated_timed_continuous_evolver` and
`approximate_timed_continuous_evolver`; the original assertions
and their respective numerical capabilities remain unchanged.
This is a compilation repair, not a relaxation of coverage or
a change in runtime semantics. CI must validate the exact resulting
commit on all five functional workflows. Format remains nonblocking.

## v6.80 — Approximate ordinary discrete inclusions using explicit tick-input samples

After the fully green v6.79 checkpoint, `DiscreteTimeEvolver<ApproximateTag>`
can now consume either `DifferenceAlgebraicEquation` or
`DifferenceAlgebraicInclusion`, through an owning `DiscreteDynamics`
metadata wrapper. The concrete difference declarations are preserved;
the numerical choice is isolated in Evolution and does not alter System.
The inclusion constructor requires an `InputSamplingPolicy`, and a
single input vector is chosen independently for each trajectory and
integer tick. There is no implicit inclusion-to-equation conversion.
The sample is dimension-checked and must be finite and inside the
declared input box. The input selected for tick n is used both in
the next-state update f(x[n],u[n]) and, when present, in the
auxiliary values g(x[n],u[n]). Auxiliary output at the final
index uses the independently selected final-index input.
No unselected input is represented by the numerical trajectory.

Initial boxes use the same paving-independent bisection sampling
as the approximate continuous evolvers: one midpoint by default;
explicit `SamplingCounts` allow more samples. Explicit initial
points remain unresampled. Tick labels are exact integer indices
converted with nearest rounding to representable `FloatDP`, with
an overflow/duplicate-timestamp check.

The `EvolutionTraits` specializations now map the native discrete
inclusion and the common `DiscreteDynamics` view to the single
approximate discrete evolver. Regression tests retain the existing
equation assertions and add deterministic time-indexed input
selection, per-trajectory choices, auxiliary correlation,
cloning/polymorphism, initial subdivision, midpoint input,
invalid input dimensions/bounds and empty initial lists.
No UPPER or LOWER semantics is introduced for the sampled path.
Hybrid-discrete Approximate and native validated ordinary discrete
inclusions remain separate pending work.

The checkpoint remains a candidate until Ubuntu, macOS, Windows,
Coverage and Sanitizers are all green on the exact commit. Format
stays nonblocking; test-target CMake settings remain uniform.

## v6.81 candidate — Branch-complete sampled hybrid-discrete evolution

`HybridDiscreteTimeEvolver<ApproximateTag>` advances a concrete point
through the native `DiscreteTimeHybridAutomatonInterface`. Each
enabled guard admits a separate transition at the following integer
index, and permissive progression simultaneously admits its own
separate branch when its progression invariant allows it. **No
implicit transition priority or selection** is imposed: multiple
activations continue as multiple point trajectories, preserving
their common prefix, event choices, initial-seed identity, and
individual locations. Each transition applies its reset directly
to the source pre-state without first applying the map. Ordinary
progression uses the mode's update; both take one tick.
A blocked path is retained with explicit `blocked` provenance.

An append-only branch-aware
`HybridDiscreteTimeApproximateOrbit` returns full independent
`HybridDiscreteTimeApproximatePath` objects. The prior fixed-size
hybrid `Orbit<Vector<HybridPoint>>` cannot dynamically fork the
number of trajectories; this specialized type avoids both silent
branch loss and artificial merging. In particular, the
specialization does not claim the generic single-point
`EvolverInterface<ApproximateTag,...>` contract, whose
result shape would preclude a fanout. No UPPER/LOWER set
semantics are claimed.

`InputSamplingPolicy` remains separate from initial sampling.
Every inclusion source selects an input sample per state/tick,
and input-dependent progression computes a single chosen
numerical image; choices are shared by transitions from that
same source path. A descendant can select a fresh sample at its
next index. Initial boxes use `sample_hybrid_initial_points`
and thus the same paving-independent `SamplingCounts`
as the continuous approximate methods (one midpoint default).
The per-call maximum path budget defaults to 10000 paths;
exceeding it aborts rather than truncating a branch.

The implementation uses the existing single-tick scheduler
for urgent/permissive guard and progression logic. Source
predicates evaluate the numerically sampled pre-state with
their declared algebraic assignments; an exactly zero
indicator is treated as unknown and all admissible alternatives
remain in the approximate exploration, rather than guessing
boundary priority. Tests cover simultaneous permissive
events, urgent competition, multi-tick continuation, blocked
paths, exact event histories, initial cardinalities, and
explicit budget failure. This is finite point-path enumeration,
not exhaustive exploration over an initial set or input set.

As with every implementation checkpoint, the five
functional workflows must pass on the exact commit for
acceptance. Format remains nonblocking. No test-specific
include directories or Windows timeouts are added.

## v6.82 — Disambiguate approximate hybrid path alias

The v6.81 build failed at compilation in
`src/hybrid_discrete_time_approximate_evolver.cpp` because the
anonymous-namespace alias `Path` shadowed the existing
`Ariadne::Path` class. Its uses in `Frontier`, the initial
path construction, a blocked leaf and a spawned branch were
therefore ambiguous; undeclared-variable diagnostics were
secondary failures. Rename exactly those five references to the
unambiguous `ApproximateHybridPath`. No algorithm, trajectory
fanout, reset, input sampling, test, CMake rule or runtime
expectation is changed. Acceptance requires Ubuntu, macOS,
Windows, Coverage and Sanitizers to pass on the exact v6.82 SHA.
Format stays nonblocking.
