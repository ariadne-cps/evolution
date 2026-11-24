# System families: simulation and validated evolution

This is the implementation roadmap adopted in design v6.20 on 22 November
2026. Planned names in this document do not declare implemented APIs.

| Native family | Validated evolution | Numerical simulation | Current status |
| --- | --- | --- | --- |
| Discrete time (`DifferenceAlgebraicEquation`) | `DiscreteTimeEvolver<ValidatedTag>` | `DiscreteTimeEvolver<ApproximateTag>` | Validated exists; approximate planned |
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
Simulators produce approximate sampled trajectories. Sampling a trajectory
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
