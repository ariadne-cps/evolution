# System families: simulation and validated evolution

## NON-NEGOTIABLE CONTRACT — Validated UPPER / LOWER semantics

**STATUS: CANONICAL AND BINDING. Read this section BEFORE any historical
design entry below, before changing code or tests, and before updating the
capability matrix. It overrides every inconsistent proposal, interpretation,
or incomplete experiment recorded in this version history. No evolver may
silently redefine these semantics.**

### Mathematical definitions (both are rigorous OUTER enclosures)

Let `R` be the actual reachable set, `S_U` an admissible superset of
`R`, `S_L` an admissible **nonempty** subset of `R`, and `E_U,E_L`
the numerically validated enclosures returned by the evolvers. The
required relationships are

```text
UPPER:  R ⊆ S_U ⊆ E_U
LOWER:  ∅ ≠ S_L ⊆ R  and  S_L ⊆ E_L
```

The **UPPER** result is a rigorous overapproximation of a superset of
the reached set; the **LOWER** result is a rigorous overapproximation
of a subset of the reached set. Numerical enclosures are **outer in
both cases**. In particular, `E_L ⊆ R` is **NOT** required and would
be an INNER approximation, a different contract. The existence of a
reachable witness inside `E_L` is the weakest nonempty-`S_L`
instance of LOWER, not permission to invent witness selection policies.
The same definitions apply in every ordinary/hybrid, discrete/continuous,
equation/inclusion family. UPPER completeness must cover the intended
reachable set and all admissible branches; a single per-branch outer
enclosure must not be mislabeled a global UPPER result.

### Dynamics and time advance — never invent a second LOWER dynamics

- **Same physical dynamics for both semantics.** Discrete time advances
  by application/composition of the original map. Continuous time
  advances the original flow using its configured validated integrator;
  UPPER and LOWER may legitimately use different integrators or
  numerical policies, but never different right-hand sides or event
  rules fabricated to make LOWER easier.
- For a discrete bounded-input inclusion `x[n+1]=f(x[n],u[n])`,
  apply the same set-parameterised map to the full current
  `X[n] × U[n]`. Each time index uses its own fresh, independent,
  **full admissible input box** when required by the model.
  Preserve input/state/history correlations and auxiliary bindings.
  **Never** substitute `u[n]=mid(U)`, a sampled input, a convenient
  input subbox, or input-equality constraints under LOWER. Sampling
  belongs only to a separately and explicitly specified policy/API.
- An **event-free** step has no need for event/guard certification:
  apply the map (discrete) or integrate the flow (continuous). Event
  analysis is **separate** from numerical time propagation.

### Set subdivision, aggregation, and provenance

- UPPER may subdivide if the union of propagated pieces still
  rigorously covers the intended source/reachable superset.
- LOWER must preserve the actual reachable subset `S_L` represented
  by the enclosure. Splitting an outer LOWER enclosure into children
  does **not** prove that each child contains any reachable point:
  `S_L ⊆ E_1 ∪ E_2` does not imply `S_L ∩ E_i ≠ ∅` for either
  arbitrarily chosen `i`.
- **Never publish an independently selected child merely because its
  guard is decidable.** Any subdivision/partition used by LOWER must
  retain the necessary aggregation or explicit, rigorous correspondence
  to the reachable subset and its provenance. Without that proof and
  bookkeeping, do not split. A numerical refinement or a proof that
  an arbitrary child of an *outer* enclosure satisfies a predicate is
  not by itself a reachable-witness certificate.
- Preexisting parameter constraints and source correlations must survive;
  no artificial restrictions or dropped branches may be justified
  solely by the desire to force a test or complete a matrix cell.

### Hybrid events and transition semantics

- Keep time propagation and event decisions distinct. The UPPER
  treatment must cover all possibly admissible alternatives as required
  by the model. LOWER may publish only continuations known to contain
  a nonempty genuinely reachable subset, consistent with guards,
  invariants, urgency/permissiveness, resets, time and provenance.
- Do **not** invent priorities between alternatives, assume an UNKNOWN
  guard is TRUE/FALSE, or transform an uncertain UPPER candidate into
  a LOWER claim. Conversely, **a unique globally certain branch is a
  conservative sufficient implementation condition, not the definition
  of LOWER**; richer branching is possible only with a rigorous
  subset/provenance argument.
- Compare every proposed hybrid-discrete LOWER mechanism with the
  already implemented hybrid-continuous LOWER behavior and the
  canonical mathematical definition here. Do not introduce separate
  witness-selection or parameter-splitting algorithms merely because
  discrete time makes guard inspection simpler.

### Mandatory change-control and verification

1. Before changing any evolver, write down the original dynamics,
   represented `S_U` or `S_L`, the mathematical inclusion invariant,
   and **why** the smallest proposed implementation change preserves it.
   If a proof is missing, do not claim support.
2. Do not weaken, remove, or rewrite an architectural test just to
   make CI green. Classify failures as compilation/test errors,
   defective expectations, production bugs, or semantic uncertainty;
   update existing tests when capabilities change, while retaining
   their intended coverage.
3. Update this canonical section **only when the user explicitly
   changes the agreed semantics**; historical entries must never
   override it. Keep the cumulative working design and repository
   journal updated after each intervention.
4. Declare a matrix capability implemented only for its precisely
   stated contract and after all five functional CI workflows
   (Ubuntu, macOS, Windows, Coverage, Sanitizers) succeed on the
   **exact** commit. Do not confuse published code with validated code.
5. No speculative LOWER restrictions, sampled input witnesses,
   independently published subboxes, new guard priorities, or
   extra feasibility gates may be introduced as substitutes for
   the already established semantics.

---

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

## v6.83 candidate — Validated hybrid discrete inclusion-target reset

The remaining partial UPPER support for hybrid discrete inclusions
is extended to allow transitions into inclusion targets with no
auxiliary outputs. For both equation and supported inclusion sources,
a reset consumes a tick but chooses no new input: the destination's
bounded inputs are introduced universally only when it next
progresses. Every recognized guard branch retains an independent
validated enclosure and its pre-state restriction. The shared
capability gate still rejects inclusion modes with auxiliary
outputs, and direct image application also rejects auxiliary-bearing
inclusion targets. This avoids silently violating correlations
between input, auxiliary values, guard and reset.

Tests cover equation-to-inclusion and inclusion-to-inclusion
resets followed by fresh-input progression, event provenance,
integer clock, no input at arrival, independent input at the
next tick, and rejection of an auxiliary-bearing target.
Validated LOWER is not added, and the hybrid discrete inclusion
UPPER cell remains partial until auxiliary-bearing inclusion modes
are supported or their exact capability boundary is formally
settled. This checkpoint requires its own five green CI workflows.
Format remains nonblocking.

## v6.87 candidate — Ordinary discrete validated UPPER for bounded-input inclusions

A dedicated `DiscreteTimeInclusionEvolver` implements
`EvolverInterface<ValidatedTag, DifferenceAlgebraicInclusion, LabelledEnclosure, Integer>`
without changing the already public strongly typed
`DiscreteTimeEvolver<ValidatedTag>` for deterministic equations.
It is registered in `EvolutionTraits<ValidatedTag, DifferenceAlgebraicInclusion>`
and exported in the public umbrella header.

The validated kernel reuses the same full-box construction already present in
hybrid discrete UPPER progression. Each tick appends a fresh independent
`EnclosureVariableKind::INPUT` coordinate over the entire declared input box,
composes the symbolic update with the correlated previous state and the new
inputs, and advances the enclosure by one exact tick. No input samples are
chosen; earlier input parameters and source constraints remain represented.
UPPER can subdivide, but it must retain every resulting branch; the
configuration and clone are independent.

This first slice rejects `Semantics::LOWER` unconditionally, and rejects
inclusions with auxiliary outputs even for the zero-tick orbit. It also rejects
negative horizons and incompatible state spaces. These exclusions are
deliberate, pending a coherent representation of input-dependent auxiliary
values and unification of the validated ordinary discrete interface.
Thus ordinary discrete inclusion / validated UPPER changes from absent to
**partial**; LOWER stays absent. No System or Integration changes are needed.

The new `test_discrete_time_inclusion_evolver` verifies independent input
dimensions and kinds over two ticks, enclosure coverage of both extreme
input sequences, the zero-horizon case, source immutability, independent clone
configuration, and explicit negative cases. Existing ordinary discrete
LOWER/equation suites and both hybrid families stay in the full CI gate.
Acceptance requires passing Ubuntu, macOS, Windows, Coverage and Sanitizers
on the exact v6.87 commit. Format is nonblocking.

## v6.88 — Multi-input correlation and conservative subdivision regression

The v6.87 UPPER ordinary difference-inclusion implementation is exercised
with two independently bounded inputs per tick and a genuinely constrained
source state. The test keeps the effective coordinate-function owner alive,
adds two source constraints, and verifies their retention across two ticks.
Four new INPUT coordinates represent the independently selected `u[n]`
and `v[n]` at both indices, retaining the single initial state parameter.
Explicit reachable witnesses with `x[0]=3/8` at the two extreme input
sequences must remain inside the final UPPER enclosure. A separate test
enables subdivision above a configured radius and demands that the resulting
union still covers both extreme one-step input/state choices; no child branch
may be discarded. The generic validated `EvolverInterface` route is exercised.

This checkpoint changes only tests and documentation, not production
numerics, System declarations, or the capability matrix. The underlying
auxiliary-output limitation remains: `a[n]=g(x[n],u[n])` would require
a representation that correlates the input selected for a given tick
with both the resulting update and its auxiliary value. The current
`Enclosure` auxiliary argument layouts support state, time and fixed
physical parameters but do not yet encode this per-tick INPUT dependency.
Do not install a state-only auxiliary mapping in its place or claim
validated support for such output modes.

Acceptance of this v6.88 regression checkpoint requires five successful
functional workflows on its own SHA; v6.87 Ubuntu, macOS and Coverage were
green at the latest observation, while Windows and Sanitizers were still
in progress. No CI result can be borrowed from another SHA. Format remains
outside the functional gate; existing discrete equation LOWER and both
hybrid-family suites stay active.

## v6.89 candidate — A bounded tutorial cover and an umbrella-only C++ contract

The installed C++ tutorial now has a **single** include:
`#include "ariadne-evolution.hpp"`. Its supporting function, hybrid,
integration, numerical, system and standard I/O types are made public
transitively through the Evolution umbrella; a consumer never needs to
name their implementation headers. The existing standalone tutorial
CMake target is compiled and run against installed headers and libraries
by the Release workflows, so it serves as a real packaging/umbrella test,
not merely a build-tree test.

Examples use only presently implemented operations and retain their
behavior as further families become available. The compact selection
covers autonomous and nonautonomous continuous equations, Approximate
sampled paths against Validated UPPER and LOWER for ordinary equations,
structured ordinary continuous/discrete inclusions with explicit
Approximate input policies against their available validated UPPER,
discrete hybrid point-path branching and UPPER set branching, and
continuous hybrid equations and inclusions. Exactly one first-order PDE
demonstration is retained, without multiplying PDE variants. The
tutorial does **not** enumerate or describe unimplemented combinations.
The future capability matrix and remaining design work belong exclusively
in these design documents, not in the tutorial.

The ordinary discrete inclusion UPPER example uses the v6.87 public
`DiscreteTimeInclusionEvolver`. Its model has no auxiliary outputs, in
line with the currently enforced numerical capability boundary.
Discrete hybrid variants use deterministic two-mode permissive branches
and distinguish the sampled branching result from UPPER set orbits.

The umbrella's direct public includes make the current tutorial
self-contained. The README is rewritten around the supported, concise
worked cases and retains separate instructions for installed C++ and
Python tutorials, without promising not-yet-implemented Python parity.
No new algorithm, System or Integration code is changed. Acceptance
requires passing Ubuntu, macOS, Windows, Coverage and Sanitizers on the
exact checkpoint, including the installed C++ tutorial build/run.

## v6.90 — Restore installed tutorial compilation after precise const-correctness diagnosis

The v6.89 installed C++ tutorial failed in the **compilation** phase on Ubuntu,
macOS and Windows Release. The Ubuntu Clang 16 log supplied by the user gives
three explicit diagnostics: `Orbit<List<LabelledPoint<FloatDPApproximation>>>::curves()`
is not a `const` member, but three local orbit results were declared
`const auto` before calling `.curves().size()`. They are
`timed_sample`, `sampled_inclusion_orbit` and `sampled_ticks`.
This is a C++ const-qualification mismatch in the newly added tutorial
code, **not** evidence of a failed validated evolver, runtime numerical
error, or CI infrastructure outage.

Fix exactly the three local declarations to `auto` while preserving the
existing `curves()` API and all tutorial functionality. The tutorial retains
its single `#include "ariadne-evolution.hpp"`; no examples, calculations,
semantic assertions, algorithms, test coverage or platform workflows
are removed or weakened. In particular, this change does not make
`Orbit::curves()` const, because that would be an unrelated public API
change requiring independent design and tests.

The v6.89 Coverage and Sanitizers checks were successful, but Ubuntu and
macOS Release failed building the installed tutorial; Windows Release
exhibited the same failed CI step. Acceptance of v6.90 requires the
five functional workflows to pass on **this exact commit**. The
pre-existing v6.88 baseline passed all five.

## v6.91 — Current-input auxiliary binding in ordinary discrete UPPER inclusions

The ordinary `DifferenceAlgebraicInclusion` contract specifies
`a[n]=g(x[n],u[n])` and `x[n+1]=f(x[n],a[n],u[n])`, with an independent
whole-box input choice at each tick. It is not correct to bind `a[n]`
to the previous input `u[n-1]` or to treat an input-dependent auxiliary
as a state-only function.

The shared `Enclosure` explicitly represents this relation with the new
`EnclosureFunctionArguments::STATE_INPUTS` argument layout. Its
`input_parameter_function(count)` chooses the most recent `count`
`EnclosureVariableKind::INPUT` parameter axes in declared order while
leaving all earlier input coordinates and constraints intact. The active
INPUT coordinates are protected against Kuhn reconditioning, so an
auxiliary's current-input binding cannot silently disappear while
historical input coordinates are reduced conservatively. Physical
constant-parameter layouts remain disjoint from the new input layout.

`DiscreteTimeInclusionEvolver` now supports explicit algebraic auxiliary
assignments. At the initial index, it appends a full independent `u[0]`
box and binds `a[0]=g(x[0],u[0])`, including on zero-horizon orbits.
Progression uses this same `u[n]` to compute `x[n+1]`, increments
the validated time by one and appends a *fresh* `u[n+1]` for the
new auxiliary `a[n+1]`. This makes the auxiliary view consistent
at every reached/final index without input sampling. An auxiliary-bearing
resume requires an existing STATE_INPUTS binding and current INPUT axes
whose domains lie inside the model's declared input box; an input-free
initial enclosure is automatically extended. Existing zero-auxiliary
models retain their original one-new-input-per-progression behavior.

The new regression traces `a=x+u`, `x'=x+a` for two ticks:
`x[1]=u[0]`, `x[2]=2u[0]+u[1]`, `a[2]=x[2]+u[2]`,
including exact current-input correlation, fresh input kinds, full
range coverage, zero horizon, immutable source and one-tick resume.
A LOWER request remains explicitly rejected. The hybrid discrete
inclusion capability gate is unchanged: auxiliary-bearing hybrid
modes still require separate analysis of guards, event resets and
the no-new-input-on-reset rule. No tutorial example is added for an
unimplemented hybrid auxiliary case.

The ordinary discrete inclusion Validated UPPER capability is now
implemented for declared structured bounded-input systems with explicit
algebraic outputs; this does **not** claim LOWER, alternate unstructured
multifunctions or comprehensive solver guarantees. Its existing
`EvolverInterface<ValidatedTag, DifferenceAlgebraicInclusion,...>`
contract remains distinct from the equation-only validated
`DiscreteTimeEvolver<ValidatedTag>`. Acceptance requires all five
functional workflows on this exact commit, including the installed
umbrella-only tutorial. At the latest pre-commit observation, v6.90
Ubuntu and macOS had completed successfully and the other three
workflows had not yet finished.

## v6.92 — Event-free hybrid discrete validated inclusion auxiliaries

The validated UPPER hybrid discrete evolver now supports algebraic outputs
`a[n]=g(x[n],u[n])` for inclusion modes **without outgoing events**.
The initial enclosure attaches the whole independent `u[0]` box even at
zero horizon; subsequent progression consumes exactly that `u[n]` for
`x[n+1]=f(x[n],a[n],u[n])` and introduces an independent `u[n+1]`
for the next output. Input correlations and state constraints are preserved
by the shared `EnclosureFunctionArguments::STATE_INPUTS` representation
introduced in v6.91. Every supported index has consistent labelled
auxiliary values. The input values are never sampled.

The private `prepare_discrete_inclusion_source` normalizes both public
hybrid-orbit entry and direct numerical step-image entry. It validates
the source auxiliary labels, active current-input binding and declared
input domains on resumed enclosures. Without outputs, input introduction
remains the previous one-new-box-per-progression mechanism. The hybrid
integer clock, UPPER fan-out, exact event provenance and source
immutability remain unchanged.

**Capability boundary:** outgoing events from an inclusion source with
input-dependent outputs are rejected, even at zero horizon. Source guard,
progression invariant and reset expressions may mention the auxiliary;
allowing them before a dedicated correlated predicate/reset implementation
would be mathematically unsound. Transitions into inclusion modes with
outputs remain rejected too. Ordinary equation modes and the existing
state-only inclusion transitions are unaffected. Accordingly the hybrid
discrete inclusion Validated UPPER cell is still **partial**, and LOWER
remains absent.

Tests exercise `a=x+u`, `x'=x+a`, including two ticks,
`x[2]=2u[0]+u[1]`, the final identity
`a[2]=x[2]+u[2]`, current-input coordinate kinds,
zero-horizon auxiliary outputs, a resumed index, source immutability
and explicit rejection when an auxiliary-bearing source has an event.
No unsupported case is added to the user-facing tutorial, which still
has exactly the sole `#include "ariadne-evolution.hpp"`.

This checkpoint is **published but requires all five workflows** on
the exact v6.92 commit before it is accepted. At v6.91 observation,
Ubuntu/macOS succeeded; Windows, Coverage and Sanitizers were pending.

## v6.93 — State-only hybrid discrete events from input-auxiliary inclusions

UPPER validated hybrid discrete inclusion modes with algebraic
`a[n]=g(x[n],u[n])` may now have outgoing events whose **guard,
progression invariant and reset expressions all depend on state alone**.
The previous blanket rejection of every event from an
auxiliary-bearing inclusion was unnecessarily restrictive:
a pre-state condition `p(x[n])` and state-only reset
`x[n+1]=r(x[n])` do not require an input-dependent auxiliary
substitution. They can be validated with the already certified
state-only predicate and reset kernels.

`validate_discrete_inclusion_source` now checks the syntactic
arguments of each outgoing guard, invariant and reset assignment
against the source state space. For a source inclusion with
auxiliaries, any reference to one of its auxiliary variables in
any of those three expression classes causes an explicit
`std::invalid_argument`, including at zero horizon and in direct
step-image entry points. No dependency is silently eliminated by
numerical sampling. Events to inclusion target modes carrying
input-dependent auxiliaries remain unsupported.

In the accepted case, the hybrid source is first prepared with the
full `u[n]` box and the correlated mapping `a[n]=g(x[n],u[n])`.
The independently restricted branches use source-state guards and
invariants. A state-only reset consumes one tick but **does not
allocate another current input**; its event provenance and target
equation-mode auxiliary labels follow the existing reset kernel.
A progression uses the same current input as `a[n]`, advances the
state and attaches independent `u[n+1]`. The two UPPER branches
remain distinct; neither can be discarded for an UNKNOWN guard.
No changes are made to the System declaration or equation dynamics.

New tests cover direct step images and a public
`HybridDiscreteTimeEvolver<ValidatedTag>::explore/orbit` over a
permissive event with state-only guard, invariant and reset:
two distinct successors, one reset and one inclusion progression,
time-one advancement, exact input coordinate counts, final image
coverage, provenance and nonmutation of the original source.
Separate negative tests exercise auxiliary-dependent *guard*,
*invariant* and *reset*, preserving a strict semantic boundary.
The existing negative test for auxiliary-bearing inclusion targets
remains active.

The UPPER hybrid discrete inclusion cell remains **partial**:
input-auxiliary-dependent events and inclusion reset targets with
auxiliaries still need a separate correlation-preserving
implementation. Validated LOWER is unchanged and absent.
The user-facing tutorial stays at one public include,
`ariadne-evolution.hpp`, without any unimplemented examples.
The journal is the only place that records future work.

CI qualification: the v6.92 commit
`90ab0fc41ce8d4f9d06d3a8ccb97c47ecea561de` had no check
runs or commit statuses and was not reflected in the PR #32 head
as of the 25 November 2026 inspection, despite being the actual
`multifunction` branch tip. The GitHub workflow triggers are
`push: main` and `pull_request`, not direct pushes to
`multifunction`. Thus v6.92 is **not** CI-validated, and neither
will v6.93 be considered accepted without genuine CI evidence on
its own SHA. Do not borrow the five green v6.91 checks.

## v6.94 — Regression repair after actual v6.93 Debug/ASan failures

The v6.93 functional CI is **red on Ubuntu, macOS, Windows, Coverage and
Sanitizers**, although several Release matrix jobs succeeded. The reported
Debug failures must not be mistaken for CI infrastructure problems.
Classification:

1. **Incorrect legacy test expectation**, predicate evaluator
   `supported_inclusion_and_rejected_cases`: the old test expected
   `classify` to throw for an event-free `DifferenceAlgebraicInclusion`
   with `a=x+u`, even for a state-only predicate `x>=0`. Since
   the v6.92/v6.93 capability supports state-only predicates in such
   modes, the correct verified contract is three-valued evidence,
   not rejection. Test TRUE/FALSE/UNKNOWN on the same constrained
   state set and retain the mismatched-state-space exception.

2. **Incorrect legacy test expectation**, UPPER scheduler
   `inclusion_rejected_without_sampling`: auxiliary-bearing
   event-free inclusions are now supported at horizon zero and
   during progression, with a full independent input box per index.
   Test the zero-horizon output, the independent input at the
   next tick, and explicit rejection at both horizons **when a
   newly declared outgoing event references the input-dependent
   auxiliary**. No artificial weakening of coverage.

3. **Crash in a new Debug test**, step-image
   `event_free_input_dependent_auxiliary`: after the expected
   STATE_INPUTS check, the test used
   `set.input_parameter_function(1u)[0u]`, indexing a temporary
   vector-function patch while composing a scalar residual.
   This risks borrowing a component whose owner has already been
   destroyed. Retain the vector patches for inputs, auxiliaries and
   states as named local values before extracting scalar components.
   Keep the exact residual check
   `a[2]-x[2]-u[2]=0` and bounds; do not skip it.
   This is a suspected lifetime root cause, **not yet proven
   until CI on this SHA passes**.

The previously described v6.93 direct step-image state-only
branch regression was also missing from the committed test file,
despite being referenced in design notes. Register it explicitly:
an inclusion with `a=x+u` and a permissive, state-only
guard/invariant/reset must preserve **both** branches, time,
event provenance, independent INPUT dimensions, reset image,
and original source immutability.

This checkpoint changes tests and their design notes, not production
numerics or hybrid capability gates. The installed C++ tutorial
continues to use **only** `#include "ariadne-evolution.hpp"`,
and contains exclusively implemented examples. Do not declare the
v6.94 fix accepted without full Debug/Release and sanitizer checks
on the exact commit. Continue to investigate any crash or
constraint failure reported by subsequent CI runs; do not change
the underlying mathematics or suppress the test for a green metric.

## v6.95 — Validated UPPER reset arrival into inclusion targets with input-dependent outputs

The hybrid discrete UPPER reset kernel now supports a transition **into**
a bounded-input `DifferenceAlgebraicInclusion` that declares algebraic
outputs `a[n]=g(x[n],u[n])`. The target input is a *new independent
whole-box input at the arrival index*; it is neither sampled nor used
while evaluating the reset.

A one-tick jump from source mode q to target mode r proceeds in
strict order:

1. Restrict the source to the closed validated guard candidate, with
   its existing source variables, constraints, history and correlated
   auxiliaries intact.
2. Construct and evaluate the reset `x_r[n+1]=R(x_q[n],a_q[n])`
   from the **source model**. Existing capability gates reject any
   source-auxiliary dependency requiring unsupported input-aware
   symbolic reset evaluation.
3. Apply the reset with an initially empty target auxiliary space,
   record the event, switch location and increment the integer tick.
   The source's historical INPUT axes remain present when they may
   correlate with the reset state.
4. If the destination inclusion has outputs, append a new independent
   `u_r[n+1]` box and bind `a_r[n+1]=g_r(x_r[n+1],u_r[n+1])`
   using the existing `STATE_INPUTS` mechanism. No target input
   influences the preceding reset, and no next-state progression
   occurs at this arrival tick.
5. A subsequent progression from r consumes precisely
   `u_r[n+1]`, computes the validated UPPER state image, and appends
   another fresh target input for index n+2.

This follows the same time semantics for equation-to-inclusion and
inclusion-to-inclusion jumps. A source inclusion with input-dependent
outputs may jump only under already supported **state-only** guard,
invariant and reset expressions. A source equation retains its
existing symbolic auxiliary substitution.

Added step-image regressions cover equation-to-auxiliary-inclusion,
auxiliary-inclusion-to-auxiliary-inclusion, correct arrival labels,
independent INPUT axis counts, exact algebraic residual correlation,
reset-state independence of new target inputs, clock advancement,
event provenance, following progression and original-source
immutability. A public hybrid validated UPPER regression verifies
the corresponding integer-time branch history and final orbit.
Earlier negative expectations that rejected all such destinations
are deliberately replaced by equivalent positive semantic coverage;
the unsupported source-auxiliary-dependent guard/invariant/reset
rejections remain active.

The matrix still shows **partial** Validated UPPER for discrete
hybrid inclusions: source event predicates and resets that
reference input-dependent auxiliary outputs remain unsupported
(and LOWER is not implemented). This patch does not introduce a
new algorithmic semantics for those unsupported branches.

The user-facing tutorial is left unchanged, with exactly the
single `#include "ariadne-evolution.hpp"` and only implemented
worked examples. Acceptance is conditional on a full passing
Ubuntu, macOS, Windows, Coverage and Sanitizers matrix for the
**exact v6.95 SHA**, including Debug and installed Release
tutorial runs. The preceding v6.94 had Ubuntu, macOS,
Coverage and Sanitizers green on observation, Windows in progress;
its checks are not transferable to the new commit.

## v6.96 — Correct exact extremal witnesses in validated UPPER reset tests

The v6.95 macOS Debug CI reported failure in
`test_hybrid_discrete_time_step_image` case
`reset_into_auxiliary_inclusion_with_fresh_arrival_input`:
the final-state lower-bound expectation was `<=0`, which
was not a reachable witness for the fixture.

For the **exact fixture** `x[0]=0`, reset
`x[1]=x[0]+1=1`, output `a[1]=x[1]+u[1]`,
`u[1] in [-1,1]`, and progression `x[2]=x[1]+a[1]`,
the final set is the entire interval `[1,3]`.
The prior check `lower_bound<=0` was mathematically
incorrect and the prior upper check `upper_bound>=2`
did not cover the attainable maximum. The corrected
test verifies **both** exact extremal witnesses:
`lower_bound<=1` and `upper_bound>=3`.
The source immutability, arrival auxiliary relation
`a[1]-x[1]-u[1]=0`, parameter kinds, event
provenance and time advancement checks are unchanged.

The public hybrid evolver's related v6.95 scenario
starts with `x[0] in [0,1]`, resets
`x[1]=x[0]+3 in [3,4]`, then uses
`a[1]=x[1]+u[1]` and `x[2]=2*x[1]+u[1]`,
with `u[1] in [-1,1]`. Its exact reachable interval
is `[5,9]`, **not** `[5,8]`.
Strengthen the public UPPER regression upper extremum
to `>=9` while keeping the lower `<=5`.
These corrections improve semantic fidelity rather
than weakening tests to satisfy a regression.

No production algorithm, interface, tutorial, or
capability matrix is changed in this checkpoint.
In particular the tutorial retains exactly one
`#include "ariadne-evolution.hpp"` and contains
only currently supported examples. Guard and reset
expressions that refer to input-dependent source
auxiliaries remain unsupported.

The v6.94 predecessor had all five workflows green
on its own SHA. The v6.95 macOS Release completed
successfully while Debug failed in this test;
the other v6.95 workflows were still progressing
at the time of investigation. The v6.96 commit
requires all five workflows (including Debug,
Sanitizers and the installed Release tutorial)
to pass on its own SHA before acceptance.

## v6.97 — Shared correlated state/auxiliary layout for hybrid discrete UPPER events

The public `EnclosureFunctionArguments` now includes `STATE_AUXILIARY`,
mapping directly to the existing validated `state_auxiliary_function()`.
This mapping is `s -> (x(s), a(s))`; for bounded-input inclusions,
`a(s)=g(x(s),u_current(s))`. It does not introduce independent
auxiliary coordinates and retains all state/input/history constraints.

The existing System functions have already eliminated explicit
algebraic dependencies from the evolution update. The new change
completes the same composition for hybrid discrete event operations:
- For an inclusion with algebraic outputs, the predicate indicator
  is built over `join(state_space, auxiliary_space)` with no second
  substitution, and its explicit argument layout is
  `STATE_AUXILIARY`. A predicate without auxiliary references
  follows the same layout. Inclusion modes with no outputs and
  deterministic equation modes retain the existing state-only path.
- The public predicate evaluator prepares the current full input
  box even when invoked directly on an unprepared initial enclosure.
  The certified step executor uses that identical prepared source
  for three-valued predicate evaluation, branch restriction and
  subsequent image construction. The restriction to a closed
  necessary pre-state candidate is performed over the *same*
  correlated `(x,a)` function as classification.
- An inclusion-source reset that references its auxiliaries now
  constructs its function over `(x,a)` and invokes
  `HybridEnclosure::apply_reset` with `STATE_AUXILIARY`,
  rather than trying to substitute input-dependent auxiliaries
  into a function of state only. The reset still consumes exactly
  the current source input. Any target inclusion input is a new
  independent full-box coordinate installed *after* the jump.
- The prior capability gate rejecting any auxiliary reference
  in the source guard, invariant or reset is removed. The System
  declaration continues to enforce that every event expression
  references only declared source state/auxiliary variables:
  free input symbols cannot bypass the auxiliary contract.

Direct numerical tests are extended for input-dependent resets and
invariants. The certified executor exercises both URGENT and
PERMISSIVE events with correlated auxiliary guards, independent
UNKNOWN branches, one necessary closed constraint per branch,
source immutability, event provenance and source/target input
dimension counts. A public orbit and scheduler regression now
admit the same input-dependent guard and demand both UPPER
successors. Direct classifier tests prove TRUE/FALSE/UNKNOWN on
`a=x+u` over the full `u in [-1,1]` input box.

The regression cases which previously expected auxiliary-guard
rejection are **replaced by positive behavior tests**, not skipped.
The old equation-only and state-only tests remain unchanged, and
hybrid LOWER remains unsupported. The installed C++ tutorial is
unchanged and continues to contain exactly the single
`#include "ariadne-evolution.hpp"` with implemented examples only.

Acceptance requires Ubuntu, macOS, Windows, Coverage and
Sanitizers on this exact v6.97 commit. At preparation time the
v6.96 Ubuntu, macOS, Coverage and Sanitizers workflows were green;
Windows was still running. Earlier CI results do not validate the
new STATE_AUXILIARY paths.

## v6.98 — Fix test call sites for explicit discrete indicator layout

macOS Clang (and other v6.97 functional builds) exposed a compile-time
regression, not a certified-numerics failure: the test
`test_discrete_source_guard_uses_shared_certificates` still passed
the complete `DiscreteSourceIndicator` to the pre-existing
`classify_source_guard` and `restrict_source_guard` interfaces.
The indicator helper was extended in v6.97 to return both
`function` and `arguments`, so those two call sites must pass
both fields explicitly. The test additionally asserts
`arguments == EnclosureFunctionArguments::STATE` for its equation
fixture, preserving a strong check on the original state-only path.

This is a test-only compatibility correction; neither production code
nor numerical assertions are weakened. The validated correlation-aware
guard/reset implementation and the installed umbrella-only tutorial
are unchanged. The 5-job functional CI gate must pass on this v6.98
SHA; prior red v6.97 results do not indicate a numerical defect until
the test builds and executes.

## v6.99 — Correlated UPPER guard/reset witness and a compact installed tutorial example

The v6.98 commit `4d0e353b91f211742cfe481ff44597d11ce92cbb`
passed all five functional workflows on 26 November 2026:
Ubuntu, macOS, Windows, Coverage and Sanitizers. This validates
the v6.97 shared `STATE_AUXILIARY` path after fixing the obsolete
test call sites in v6.98, including the installed C++ umbrella-only
tutorial.

The v6.99 regression strengthens
`test_correlated_auxiliary_permissive_guard_and_invariant` without
modifying production code. For source `x[0] in [0,1]`,
`u[0] in [-1,1]`, `a[0]=x[0]+u[0]`, a permissive event
has guard `a>=0`, progression invariant `a<=0`,
reset `x[1]=a[0]+3`, and progression `x[1]=x[0]+u[0]`.
The two restricted UPPER successor sets must preserve precisely
the *same* source input coordinate in state and auxiliary functions:
the parameterized reset residual `x_reset[1]-x[0]-u[0]-3`
and progression residual `x_progress[1]-x[0]-u[0]`
must each enclose zero within a narrow numerical tolerance.
This explicitly rejects treating the guard's auxiliary value
and the reset's auxiliary value as independent samples.

Necessary closed constraints on each branch are retained; the
jump state enclosure must cover both attainable extrema `3`
and `5`, while the progression must cover `-1` and `0`.
The test also checks the historical `INPUT` coordinate in both
branches and the fresh next-tick `INPUT` coordinate introduced
only on the progression, together with provenance and clock checks.
All previous semantic tests remain enabled.

The installed C++ tutorial adds **one** example of this now
implemented and validated feature, reusing the existing
hybrid discrete evolver syntax rather than multiplying
paradigms and system families. It shows a bounded-input
`DifferenceAlgebraicInclusion` with `a=x+u` and one permissive
event whose guard, invariant and reset refer to `a`.
The tutorial still has exactly one include,
`#include "ariadne-evolution.hpp"`, and includes no examples
of unsupported future capabilities. The tutorial README lists
the new case as a single representative C++ syntax pattern.

The public capability matrix is unchanged. This checkpoint
adds a stronger regression and a tutorial example, not a
new numerical backend. Acceptance still requires a new
five-workflow pass on the exact v6.99 SHA, including
the installed Release tutorial on each platform.

## v7.00 — Validated LOWER, event-free ordinary equation modes in a discrete hybrid automaton

The v6.99 reference SHA `36af0b60c086673e61d35d92dac5c7845dadf3c1`
passed all five functional workflows on 26 November 2026:
Ubuntu, macOS, Windows, Coverage and Sanitizers. The newly added
auxiliary-dependent guard/reset tutorial also built and ran against the
installed `ariadne-evolution.hpp` umbrella on supported Release jobs.

This increment adds a deliberately narrow Validated LOWER path to the
public `HybridDiscreteTimeEvolver<ValidatedTag>::orbit()`:

* The initial location must have deterministic
  `DifferenceAlgebraicEquation` dynamics and **no outgoing events**.
  An input-bounded inclusion or an event-bearing source is rejected,
  even for a zero-tick request, rather than deriving a witness from a
  potentially spurious UPPER branch.
* The initial state enclosure must be **provably nonempty**.
  Undecidable nonemptiness and certified emptiness are rejected.
  No independent input selection, branch splitting, arbitrary priority,
  or UPPER-to-LOWER result conversion is used.
* The caller-provided initial state/auxiliary labels must match the
  equation. The zero-horizon result also exposes the equation's
  algebraic auxiliary mapping. The exact initial global time range
  must equal the supplied integer starting index.
* At each tick, the existing
  `HybridDiscreteTimeStepImage::apply` deterministic PROGRESSION
  branch composes the original state map and advances the validated
  absolute time exactly once. The initial witness remains feasible
  along the unique deterministic evolution. The returned
  `Orbit<HybridEnclosure>` contains its initial enclosure, one
  reach per tick, intermediate enclosures before the final index,
  and exactly one horizon-final enclosure.
* The current `maximum_nodes` limit bounds all generated nodes,
  including the initial node. Aborting at the limit produces an
  exception, not a truncated purported LOWER witness.
  `explore()` and `orbit_with_blocked()` remain UPPER-only.

This follows Ariadne's **witness-preserving** definition of LOWER;
it does not assert that every point of the returned enclosure is
reachable (an INNER approximation). In particular, nonlinear
Taylor enclosures may contain spurious numerical points while still
preserving one actual trajectory from a proved feasible source.

New public evolver tests include a deterministic equation with
auxiliary `a=2*x`, 0- and 3-tick horizons, time index,
algebraic output rebinding, initial-set immutability, absolute-index
resume, node-budget exhaustion, and a concrete `x[0]=1/2`
witness reaching `x[3]=3.5`. Rejection tests check an outgoing
permissive transition, an uncertain-input inclusion (even at zero
horizon), and a mismatched absolute starting index. The old
zero-horizon LOWER rejection for a deterministic event-free equation
is replaced by a positive regression.

The public matrix changes **validated hybrid discrete equation LOWER**
from absent to **partial (event-free equation modes only)**.
All hybrid inclusions, guarded branches, resets, and any case lacking
a proved initial witness remain outside LOWER. UPPER is unchanged.
The installed C++ tutorial is unchanged and continues to include
only `ariadne-evolution.hpp` and already implemented examples;
future capabilities are recorded here, never in the tutorial.

The v7.00 feature is published as an initial implementation, not
accepted until five functional workflows pass on its exact commit,
including Debug test runs and the installed Release tutorial.

## v7.01 — Repair validated LOWER initial witness admission

The v7.00 public hybrid discrete LOWER regression failed twice:
both lower_event_free_equation_witness and rejections raised the
"requires a provably nonempty initial enclosure" exception
for a plain unconstrained HybridExactBox with x in [0,1].

Root cause: negating Enclosure::is_empty() is not a valid positive
existential certificate. The existing emptiness semi-decision may
return INDETERMINATE even for a nonempty unconstrained exact box.
An UNKNOWN result is not evidence that a witness is unavailable.

The repaired LOWER gate uses an explicit constructive sufficient
condition on the reduced exact parameter domain:
- it must be provably nonempty;
- with no constraints, any parameter point is an initial witness;
- with constraints, the validated image of every parameter
  constraint over the *whole* reduced domain must be provably
  contained in its declared bound interval.

This rule also accepts some constrained sources when all their
constraints can be certified throughout a nonempty domain,
but it does not pretend to be a complete existential feasibility
solver. Other constrained sources are rejected conservatively
rather than silently admitted. It preserves the LOWER contract:
at least one genuine trajectory from the initial source exists.

The has_certified_lower_source helper resides in an anonymous
namespace in the discrete hybrid evolver implementation.
The deterministic progression, absolute clock, auxiliary mapping,
node budget, existing negative tests for event branches and
inclusions, and installed umbrella-only C++ tutorial remain intact.
The regression now explicitly checks the original exact-box
parameter domain and its zero-constraint status.

The v7.00 tests were red on the reported LOWER cases, so the
fix must be accepted only if all five v7.01 workflows pass
on its exact commit. The installed tutorial continues to include
only ariadne-evolution.hpp and has no unsupported examples.

## v7.02 — Constructive singleton witnesses for constrained hybrid discrete LOWER

The predecessor v7.01 `e4b53752087eb9f835cdb962ae83ec8b7c3d33d6`
has Ubuntu, macOS and Coverage green. At this checkpoint, its
Windows and Sanitizers workflows were still running: v7.01 should
not be described as fully accepted until those exact checks finish.

The preceding admission helper correctly rejects the negation
of an inconclusive `Enclosure::is_empty()` result as a witness,
but was unnecessarily incomplete for constrained initial sources:
it required **every** point in the entire reduced parameter
box to satisfy every constraint.

v7.02 first accepts that strong whole-domain certificate when
it holds, and otherwise attempts a **single constructive parameter
witness**. Provided the exact reduced domain is provably
nonempty and bounded, construct its coordinatewise midpoint
as an exact singleton parameter box. Evaluate the *full vector*
of existing constraints on this singleton with validated
interval arithmetic. Accept only when the result is
provably contained in all declared constraint bounds.

Such a proof certifies at least one feasible parameter, hence
a real initial state and a witness-preserving trajectory under
the existing deterministic equation update. Unknown comparisons
are not converted into truth. The source enclosure, its full
constraints, state function, time, auxiliary mapping and history
remain unchanged: the singleton is used only for admission,
not to replace the result with an arbitrarily chosen sample.
If midpoint verification fails, the conservative capability
gate still rejects the source. This is not a general nonlinear
constraint feasibility algorithm.

New numerical regression:
- `x[0] in [0,1]` with `x>=1/4` and `x<=3/4` has
  a certified interior midpoint `x[0]=1/2` even though
  the entire original domain does not satisfy both constraints.
  A zero-horizon LOWER result and a two-tick
  `x[n+1]=x[n]+1` evolution must exist and preserve both
  constraints; the final validated set must cover the concrete
  witness `x[2]=2.5`.
- `x[0] in [0,1]` with the impossible constraint
  `x>=2` must be rejected at zero horizon and one tick.
- Existing event/inclusion negative regressions and all
  unconstrained deterministic LOWER tests remain enabled.

No transition-level LOWER witness inference is added; the
hybrid discrete equation LOWER capability remains **partial**.
The user-facing C++ tutorial is unchanged, with exactly the
single `#include "ariadne-evolution.hpp"`, and only
implemented examples. v7.02 acceptance requires all five
functional workflow results on this exact commit.

## v7.03 — Witness-preserving LOWER for a uniquely forced urgent discrete hybrid branch

The predecessor v7.02 SHA `b366c12aac4ba536dfffe0ffd7386a0535b02973`
had Ubuntu, macOS and Coverage green when the v7.03 change
was prepared. Windows Debug and Sanitizers were still in progress,
so no full five-workflow pass is claimed for v7.02 at this stage.

The public `HybridDiscreteTimeEvolver<ValidatedTag>::orbit` LOWER
path is extended beyond an event-free deterministic equation mode
to a limited, rigorously certified discrete hybrid event sequence.

An initial parameter feasibility witness must first be certified
using the v7.02 admission helper. At each tick:

- The source mode must have a `DifferenceAlgebraicEquation`
  and **no events or exactly one URGENT event**. An inclusion,
  a permissive event, or multiple outgoing events causes a
  direct invalid_argument error; there is no hidden priority.
- No outgoing event forces the deterministic progression.
- For a unique URGENT event, the existing common validated
  `HybridDiscreteTimePredicateEvaluator` classifies its guard
  on the entire current correlated enclosure.
  Uniform `TRUE` forces the reset to a deterministic equation
  target; uniform `FALSE` forces the progression in the current
  mode. Both outcomes preserve at least one actual trajectory.
- `UNKNOWN` never allows selecting either alternative for
  LOWER. A single possible or closure-restricted UPPER branch
  is not an existence certificate. The LOWER call raises
  invalid_argument rather than sampling or treating a
  three-valued result as Boolean truth.
- The selected step is applied using the existing validated
  `HybridDiscreteTimeStepImage::apply`, preserving the original
  parameter constraints, existing state/auxiliary correlations,
  actual reset event provenance, and absolute one-tick time.
  Recompute evidence on the new location at the next tick.
- The node budget remains a hard failure, with no incomplete
  result returned; `explore` and `orbit_with_blocked` stay
  explicitly UPPER-only.

This is a **sufficient** LOWER certificate, not a general
event-existence solver. A uniformly true guard alone is not a
LOWER certificate for an empty source: the pre-existing initial
witness is required first. The only supported reset target
is another deterministic equation mode.

Public regression tests construct one urgent guard `x>=-1`
which is certainly true for `x[0] in [0,1]`: the forced reset
`x[1]=x[0]+3` followed by a target-mode progression
`x[2]=x[1]+1` must retain the witness `x[0]=1/2`,
`x[2]=4.5`. Location, one historical event, clock,
reach/intermediate/final counts, and source immutability are
checked. A separate urgent guard `x>=2` is certainly false
on the initial box: it forces progression rather than a jump.
A guard `x>=1/2` is uncertain on that box and must throw,
not return an artificial LOWER successor. Existing negative
tests for permissive events and bounded-input inclusions remain.

The module's public capability matrix now records **validated
hybrid discrete equation LOWER = partial (event-free or uniquely
forced URGENT branch sequences)**. The numerical UPPER path,
integration/system modules and C++ tutorial remain unchanged.
The tutorial still contains only
`#include "ariadne-evolution.hpp"`, with no unsupported
examples added. v7.03 is accepted only after full CI on the
exact SHA, not based on prior builds.

## v7.04 — Share certified branch planning between UPPER and LOWER

The v7.03 SHA `9b7ccd2a962ec165c83a70de42b271f335cf9cbf`
passed all five CI workflows (Ubuntu, macOS, Windows,
Coverage and Sanitizers). Thus the earlier restricted
single-URGENT-event LOWER evolution is fully verified.

The previous LOWER path duplicated the UPPER planner's
branch decisions by handling just one URGENT event manually.
This created an artificial exclusion for possible models with
several events or PERMISSIVE events whose alternatives can still
be proved unambiguously from the complete source enclosure.

The new helper `unique_certified_lower_branch`, in an
anonymous namespace local to the discrete hybrid evolver,
delegates all event predicate evidence to the existing
`HybridDiscreteTimeStepPlanner::plan`, using the existing
validated `HybridDiscreteTimePredicateEvaluator::classify`.
At each positive-index tick, LOWER requires:
- a deterministic difference equation in the current mode;
- exactly one branch in the complete three-valued planner
  result, with `branch.certain == true`;
- a reset target with equation dynamics if that unique branch
  is a transition.

All of these are sufficient LOWER conditions: a valid
feasible source witness, proven uniform preconditions, and a
deterministic map or reset preserve at least one concrete
trajectory. No priority is introduced, and a branch with merely
possible/UNKNOWN evidence cannot become a LOWER witness.

The planner already implements all event kinds and combinations:
- URGENT progression uses the logical complement of its guard.
- PERMISSIVE progression uses its explicitly declared invariant.
- Multiple events may coexist when all but one transition are
  ruled out, and there is no simultaneously possible progression.
- Competing guards, a simultaneous permissible progression,
  or an undecided single branch return zero or multiple
  admissible/certain options and are rejected, never selected
  according to enumeration order.
- At zero horizon no event is executed, so only the initial
  certified feasible set is required; a mode's future
  undecided guard does not invalidate this zero-tick witness.

New public LOWER tests verify:
(1) a PERMISSIVE transition with guard TRUE and progression
invariant FALSE forces a reset, preserving the historical
event and a two-tick witness `x[2]=4.5`;
(2) a PERMISSIVE guard FALSE with invariant TRUE forces the
ordinary progression without recording an event;
(3) a PERMISSIVE guard and invariant both TRUE leave two
alternatives, causing conservative rejection;
(4) two URGENT outgoing events with exactly one TRUE guard
admit that uniquely forced event and record its provenance;
(5) two TRUE urgent guards must be rejected as competing
alternatives; and
(6) an undecided urgent guard is rejected for a positive
horizon but admitted for a zero horizon.

The existing constrained-domain witness admission from
v7.02, hard node budget, equation-only reset targets,
and source immutability are unchanged. Bounded-input
inclusions remain excluded from LOWER. `explore` and
`orbit_with_blocked` remain UPPER-only. The tutorial C++
is unchanged and continues to include only
`#include "ariadne-evolution.hpp"`.

The public capability matrix stays **partial** for validated
hybrid discrete LOWER, now described as `deterministic
equation modes with a uniquely certified planner branch at
each tick`. Full v7.04 validation requires the five
workflow results on the exact new commit, not the v7.03 ones.

## v7.05 — Constructively certify LOWER inside an undecided guard source

Predecessor v7.04 `e97155f849cf7f200edb37000e1e61e7b4ca4a2a`
had Ubuntu and macOS fully green, while Windows, Coverage and
Sanitizers were still running during preparation of this change.
The five-workflow acceptance gate remains unchanged.

A uniformly decided guard over the *entire* enclosure was an
unnecessarily strong existential requirement for LOWER: for
x[0] in [0,1], the URGENT predicate x>=1/2 is UNKNOWN over the
full source even though x[0] in [3/4,1] provides a genuine
strictly enabled transition and a feasible initial witness.
Unlike UPPER, LOWER needs one actual trajectory, not a complete
fan-out or an inner approximation of all sources.

The v7.05 helper `certified_lower_step` searches a bounded
binary refinement of the current enclosure (maximum depth four),
using the existing `HybridEnclosure::split` operation.
Candidates retain all existing parameter constraints,
history, time and state/auxiliary map correlations.
An accepted candidate must satisfy *both*:
1. `has_certified_lower_source` proves that candidate has
   at least one feasible parameter witness, by a nonempty
   unconstrained domain, whole-domain constraint certificate
   or validated singleton-point certificate;
2. the common `HybridDiscreteTimeStepPlanner` returns exactly
   one possible branch, marked certain from validated guard
   and invariant evidence. An equation target is required for
   a transition.

The same **refined source enclosure** is used for the chosen
validated image; otherwise a branch shown certain on a subbox
could incorrectly be applied to the larger original source.
If feasibility is undecidable on a broad box, descendants
are still searched: lack of a positive certificate does not
prove infeasibility. Certifiably empty reduced domains are
pruned. The algorithm never samples to decide branch truth,
never converts UNKNOWN into TRUE, and imposes no event priority.
A bounded search failure throws rather than fabricating LOWER.

The existence proof remains explicitly **existential** and
witness-preserving. Earlier reach/intermediate orbit enclosures
may be broader than their later successful witness subbox,
without pretending they are inner approximations. Future
subboxes are derived from the actual current enclosure, so
one feasible continuous parameter history is preserved
across the selected steps and resets. Inputs/inclusions
remain excluded.

Regression updates replace the former demand that an UNKNOWN
whole-source urgent guard always throw with a positive
one-tick result on its strict enabled interior; check event
provenance, target location, time and a feasible reset
witness x[0]=7/8 mapping to x[1]=31/8.
A constrained source x[0]>=3/4 similarly remains a valid
LOWER input with constraints preserved. A singleton boundary
source x[0]=1/2 cannot certify either strict uniform
alternative and must still be rejected by the bounded search.
Earlier competing-branch, forced-transition, zero-horizon,
node-budget and source-immutability tests stay enabled.

The public capability remains *partial* Validated LOWER for
deterministic hybrid discrete equation modes, now including
bounded certified existential refinement of ambiguous
whole-enclosure event predicates. No tutorial example is
added; the installed C++ tutorial retains exclusively
`#include "ariadne-evolution.hpp"` and implemented patterns.
All five workflows must pass on the v7.05 SHA before
claiming this capability is validated.

## v7.06 — Diagnose a post-reset refined LOWER bounding-box SegFault

The v7.05 SHA `29c4efc5a4052e1470513365a652df6e1d7fad9a`
is **red** on Ubuntu and macOS Debug tests, Coverage, and
AddressSanitizer. The reported regression:
`test_hybrid_discrete_time_evolver`, case
`lower_refines_undecided_urgent_guard` (13/17), segfaulted
after check 13.4 verifying the event identity. The next
statement computes the reset successor's state bounding box.
This is a **runtime crash**, not an expected numerical mismatch.
The v7.04 SHA `e97155f849cf7f200edb37000e1e61e7b4ca4a2a`
has completed all five CI workflows successfully, so the
regression is isolated to the v7.05 addition or a newly
exercised downstream operation.

The working hypothesis is that the temporary
`ValidatedConstrainedImageSet` built by
`HybridEnclosure::state_bounding_box()` interacts badly
with a function patch and parameter domain restricted by
`HybridEnclosure::split()`. The exact failure point within
that statement is *not proven* from the last-check indicator,
which only records the most recently completed assertion.

The v7.06 candidate eliminates this extra temporary object:
`state_bounding_box()` now evaluates the already-owned
`Enclosure::state_function()` directly over
`Enclosure::reduced_domain()`, returning its validated
interval range. This is mathematically and numerically
the same outer-image computation performed by the
former `state_set().bounding_box()`: the latter constructs
`ValidatedConstrainedImageSet(domain, state_function,
constraints)` with its initial reduced domain equal to the
enclosure's domain, whereas the enclosure's reduced domain
can be more precise. In fact the new calculation is at
least as precise and never ignores reduced-domain constraints.
No guard, reset, planner, witness search, or tutorial
semantics are changed.

To prevent hiding the original regression, the failing
test *still executes the public state_bounding_box()* on
the refined reset successor and verifies that its interval
covers the concrete witness x[0]=7/8, x[1]=31/8.
Before indexing, the test now requires one state component,
a one-dimensional state function, a nonempty reduced
parameter domain, and a one-dimensional resulting box.
These checks distinguish malformed-image failures from
errors during bounding-box computation and prevent
out-of-range indexing from masking the source of failure.
All other original numerical and event-provenance
assertions remain enabled.

The v7.06 change is a **candidate fix**, not a verified
root cause or a passing regression until its five
workflows finish on the exact new SHA. If the failure
persists, inspect the new last-check position and the
sanitizer stack trace before changing production logic
or expected numerical values. The installed tutorial
continues to include only `ariadne-evolution.hpp`.

## v7.07 — Candidate SegFault repair via exact-subbox certificates without patch restriction

Both v7.05 `29c4efc5` and v7.06 `00c0a142` fail the
`lower_refines_undecided_urgent_guard` test with SegFault
after breadcrumb 13.4 (event-history identity). The testing
runtime publishes that breadcrumb at the **start** of an
assertion, so the fault occurs before the next assertion begins,
but its exact stack frame is not determined by the breadcrumb.

The v7.06 global alteration of
`HybridEnclosure::state_bounding_box()` did not fix the
crash and is reverted exactly to the implementation that
passed v7.04. We neither remove the failing assertion nor
weaken numerical correctness.

The v7.07 candidate avoids `HybridEnclosure::split()` and
`HybridEnclosure::restrict()` in LOWER witness refinement.
Those operations restrict and re-create validated state,
time, dwell-time and constraint function patches, which
are suspected of becoming invalid on a subsequent reset.

Instead search a bounded tree of **exact parameter subboxes**
(maximum four bisection levels) without restricting any
validated function patch. For each candidate, certify
existence via `has_certified_lower_source_on_box()`, using
the original constraints on that box; classify guards and
invariants with the shared planner by interval-validating
their indicators on the SAME box. Any unresolved predicate
remains UNKNOWN; only a single possible, uniformly certain
branch can be selected. The selected box is then imposed
upon a fresh copy of the original enclosure as
parameter-space interval constraints on each narrowed
coordinate; the validated state/auxiliary maps keep their
original domain, and the reset composes those maps while
preserving all original and new constraints. No numerical
sampling or event priority is introduced.

The initial zero-horizon existence test uses the same
bounded exact-subbox search, without deciding events.
An inability to prove feasibility of a parent box does
not imply its children are empty. Unresolved searches
raise an exception rather than fabricate LOWER trajectories.

The crash test continues to require correct event history,
mode, state dimension, reduced domain, public bounding
box and witness `x[0]=7/8 -> x[1]=31/8`. The constrained
variant now requires **both** original `x>=3/4` and the
selected exact parameter-subbox constraint to survive
the reset, as well as the same numerical witness.
The boundary singleton remains a conservative rejection.

This is a **candidate** repair; the actual crash cause
still requires confirmation by the precise v7.07 CI
and, if red, an AddressSanitizer stack trace.
The tutorial remains unchanged and includes only
`ariadne-evolution.hpp`.

## v7.08 — Fix compilation of validated constraints for refined parameter boxes

The v7.07 candidate `f46250b8b24fa96af2d7e49c0d7448cc3db8e21e`
failed compilation at `src/hybrid_discrete_time_evolver.cpp`,
inside `constrain_lower_source()`. Clang reports that the
three-argument `ValidatedConstraint` constructor cannot accept
`FloatDP` lower and upper bounds directly: neither implicit
`FloatDP -> ValidatedNumber` conversion nor the alternative
precision-parameter constructor is available. This is a
**C++ type mismatch**, not a failed mathematical LOWER assertion
and not a recurrence of the v7.05 runtime SegFault.

The new candidate converts each exact binary floating-point
endpoint through `Dyadic`, then `ExactNumber`, then
`ValidatedNumber`, before building the same interval constraint.
The values of the lower and upper bounds are preserved exactly
in the dyadic representation; no approximation to `double`
or deliberate weakening of the constraints is introduced.

Only this source construction and the cumulative design journal
change. The subbox witness proof, validated predicate evidence,
original state-function patch domains, constraint propagation,
existing event history, numerical tests and umbrella-only C++
tutorial are left intact. In particular the runtime SegFault
is **not considered fixed** by passing compilation: v7.08
requires all five CI workflows on its exact commit, including
the `lower_refines_undecided_urgent_guard` case and sanitizers.

## v7.09 — Diagnostic separation of a crashing LOWER refined event-history check

The v7.08 `0d08d566` Clang Debug run reaches
`lower_refines_undecided_urgent_guard` but reports a SegFault with the
last started check
`lower.final().front().previous_events()[0u] == event`.
In `ariadne-cps/testing`, the check-start notification precedes the
predicate invocation. This breadcrumb does **not** establish that the
comparison completed or that the subsequent `state_bounding_box()`
caused the crash. It also does not identify which object is invalid.

The test now names the final enclosure and its event history explicitly,
requires a one-element history, and separately compares the stored event
name with the expected event name **before retaining the original direct
event comparison**. These checks are diagnostic additions, not substitutions;
all witness, constraint, time, provenance, bounding-box and rejection
assertions remain. No production kernel, API, UPPER semantics, tutorial,
System or Integration code changes. The next crash breadcrumb or a sanitizer
stack trace must distinguish event storage access from later failures.

This is a **diagnostic candidate**, not a claimed SegFault repair. Acceptance
requires the five functional CI workflows (Ubuntu, macOS, Windows, Coverage,
Sanitizers) green on the exact new SHA; any failure is to be diagnosed rather
than bypassed or compensated for by removing the expectation.

## v7.10 — Correct LOWER witness expectations without imposing event priority

The v7.09 macOS test result is **16 passed, 1 failed**, with no
SegFault: the case `lower_refines_undecided_urgent_guard` reports
`lower.final().front().location() == target_mode` at line 550 as a
failed expectation, then terminates at the new
`recorded_events.size() == 1u` requirement. These results are
consistent with a legal, event-free progression, rather than a corrupt
event history. The previous SegFault arose when a failing history-size
check did not stop the test from indexing an empty list; it was not
proof that reset memory was damaged.

For the original x[0] in [0,1] and URGENT guard x >= 1/2, both
branches have independently certifiable strict interiors: x in [0,1/4]
forces progression and x in [3/4,1] forces reset. The bounded breadth-
first search splits the left child first, so the whole-source LOWER call
may select progression. Requiring the target mode and one recorded event
for that call was **an erroneous test expectation**; forcing a reset
in production would introduce an illegitimate priority and contradict
the existing existential LOWER contract.

The regression now preserves the unrestricted initial source and accepts
either certified alternative, checking that the mode and event history
agree and that a refined parameter constraint survives. It separately
exercises both boundary-inclusive one-sided initial sets: [1/2,1]
requires a refined urgent reset, retaining all original checks of
target mode, event identity, state dimension, public bounding box and
witness x[0]=7/8 -> x[1]=31/8; [0,1/2] requires refined progression
without an event, checking the public bounding box against
x[0]=1/8 -> x[1]=9/8. The constrained x>=3/4 variant, its two
preserved constraints, zero-horizon checks, immutable sources and the
rejection of a boundary singleton remain in place.

This changes **only the test and the design journal**. Production
LOWER/UPPER algorithms, System, Integration, the C++ umbrella-only
tutorial and all other test scenarios remain unchanged. Classification
UNKNOWN never becomes TRUE or FALSE without a validated subbox proof.
The v7.10 candidate is published but must pass Ubuntu, macOS, Windows,
Coverage and Sanitizers on its **own exact SHA**, including Debug and
the installed tutorial; no CI result is inherited from v7.09.

## v7.11 — Fair multidimensional parameter refinement for validated hybrid LOWER

The v7.10 commit `aa32618e` has passed Ubuntu and macOS at the latest
observation; Windows, Coverage, and Sanitizers were still in progress.
No CI result from v7.10 certifies this new candidate.

In `enqueue_parameter_halves()`, the bounded exact-subbox search always
bisected the first nondegenerate parameter. For multidimensional sources,
a guard depending on another parameter could remain UNKNOWN throughout
four refinement levels, even when one certifiable branch existed. This is
incompleteness of the sufficient LOWER witness search, not proof of an
unsound lower image. The splitter now starts at `depth % dimension`
and wraps around nondegenerate axes. The budget of four bisection levels,
validated source existence test, strict predicate evidence, shared tick
planner, and exact subbox passed to the parameter constraints all remain
unchanged. There is no priority between reset and progression and no
sampling of uncertain inputs.

The public regression `test_lower_refines_nonleading_guard_parameter`
uses state (x,y), with x declared first and guard URGENT y >= 1/2.
The initial y intervals [1/2,1] and [0,1/2] touch an undecided
boundary, so their interiors require a split of the **second** parameter.
The test independently proves a one-tick reset with event and
witness (1/2,7/8) -> (7/2,7/8), and a one-tick progression without
event with witness (1/2,1/8) -> (3/2,1/8). It retains exact
parameter constraints, source immutability and validated public
bounding-box checks. The boundary singleton y=1/2 remains rejected.
The previous single-parameter tests continue unchanged.

Only the LOWER axis-selection logic, one registered public regression,
and this design journal are changed. Upper semantics, System, Integration,
tutorials, and numerical enclosure representations are unchanged. The
capability matrix is unchanged: hybrid discrete deterministic Validated
LOWER remains **partial**, not a general satisfiability or reachability
algorithm.

**Acceptance:** Ubuntu, macOS, Windows, Coverage, and Sanitizers must
all pass on the precise v7.11 SHA, including Debug and the installed
C++ tutorial. This is a candidate implementation, not a CI-validated
claim.

## v7.12 — Existential validated LOWER for ordinary discrete bounded-input inclusions

**Base:** v7.11 `52461a8f`, whose Ubuntu, macOS, Windows, Coverage and
Sanitizers workflows all completed successfully on the exact SHA.

The native `DiscreteTimeInclusionEvolver::orbit(..., Semantics::LOWER)`
previously rejected every request, including zero horizon. This checkpoint
adds a **strictly existential**, validated witness path for ordinary discrete
bounded-input inclusions. It is *not* a general INNER set approximation:
a validated enclosure can be an outer numerical representation, while the
path is justified by constructive existence of an actual initial state
and a sequence of admissible exact input choices.

Admission requires a nonempty exact reduced initial parameter box and either
no constraints, validation of every constraint over the entire box, or a
validated feasible midpoint. UNKNOWN constraint feasibility is rejected,
not promoted to a witness. Every nonempty bounded input box is handled by
introducing a fresh independent INPUT parameter at each integer tick and
constraining that parameter to the exact dyadic midpoint of its declared
interval. The equality constraint is preserved alongside all original
state constraints and previous input-history choices; the validated map
and symbolic auxiliary functions remain composed over these correlated
parameters. The arrival input u[n+1] is bound independently of u[n].
For auxiliary-bearing dynamics, the current input is introduced and bound
already at zero horizon, and the arrival input is bound after progression.
Initially auxiliary-bearing resume enclosures are rejected for LOWER:
unknown preexisting input constraints need a separate feasibility proof.
UPPER keeps full-box universal propagation and its existing subdivision
setting without changes.

The C++ public tests retain all previous UPPER assertions. They establish
LOWER for two independent input channels with x[0]=1, exact midpoint
choices (u,v)=(2,1) per tick, and a two-tick witness x[2]=7; they
also retain an initial x>=1/4 constraint and witness x[2]=13/2.
A separate input-dependent auxiliary test certifies x[0]=0,
u[0]=u[1]=u[2]=2, x[1]=2, x[2]=6 and arrival auxiliary
a[2]=8. Both verify parameter/constraint counts, input independence,
source immutability and public bounding boxes. Former blanket LOWER
rejections are replaced by positive admission checks while **keeping**
negative tests for an infeasible constrained source and an unsupported
auxiliary-bearing resume. No runtime has yet been asserted to pass these
new tests.

**Capability matrix effect:** ordinary discrete inclusion / Validated LOWER
moves from *absent* to *partial*. It remains a sufficient existential
certificate only, not a complete constraint solver or an inner enclosure.
All other matrix cells are unchanged. No changes to System, Integration,
hybrid dynamics, tutorial, or approximate simulation. The new SHA must pass
all five standard workflows (Ubuntu, macOS, Windows, Coverage, Sanitizers),
including Debug and installed tutorial; no prior run counts for acceptance.

## v7.13 — Audit and correction of validated LOWER map semantics

**Binding design rule (all evolver families).** The time-propagation
operation is driven by the same original dynamics under UPPER and LOWER.
For discrete time, the propagation step is the application of the map
itself; there is no continuous integrator or LOWER-specific input-choice
policy. For a discrete inclusion `x[n+1]=f(x[n],u[n])`, both semantic
modes compose the set-parameterised update with **all** `u[n] in U[n]`,
keeping the complete input box as a fresh independent parameter domain
at each tick. An auxiliary `a[n]=g(x[n],u[n])` uses the current complete
box; a later `u[n+1]` has a new coordinate. Preserve source constraints,
parameter correlations, the validated enclosure representation and the
original dynamics. Semantics may affect subdivision, stopping rules,
branch admissibility, time-horizon policies and the choice/contract of
a continuous integrator; it must **never silently fix, sample or optimise**
an inclusion input unless the model or an explicit separate API declares
that restriction. For a discrete map, UPPER/LOWER are not different
numerical right-hand sides. Do not infer an arbitrary INNER set claim
from the representation: the numerical guarantee and semantic
continuation contract must be documented separately and verified.

**Violation and correction.** v7.12 `15b898c` introduced
`bind_lower_input_witness()`, adding artificial midpoint equality
constraints to every LOWER input coordinate; it also introduced an
initial midpoint-feasibility gate and an unjustified ban on resuming
auxiliary-bearing sources. All these v7.12-only behaviors and their
midpoint-witness regression expectations are **removed**. The three
modified implementation/header/test files are rebuilt directly from
the full-input v7.11 base, changing only the semantics acceptance gate,
the UPPER-only subdivision policy, documentation and precise semantic
regressions. The original UPPER path and its C++ tests are retained.
The old v7.12 journal entry stays visible only as historical record;
its claimed LOWER input-selection contract is explicitly superseded.

**Cross-family audit against the current source (26 November 2026).**
- Ordinary validated discrete equations: `src/discrete_time_evolver.cpp`
  calls `apply_discrete_time_map_step(_system->function())` for both
  semantics; UPPER subdivides optionally, LOWER has its own optional
  termination policy. No input selection was found.
- Ordinary validated discrete inclusions: **v7.12 was a confirmed
  violation**; corrected here to use the unchanged whole-input map
  composition `compose(update_function(), join(state_function(), inputs))`
  for both semantics. Earlier v7.11 supports only UPPER.
- Validated hybrid discrete equations: `HybridDiscreteTimeStepImage::apply`
  applies the **same** equation map or transition-reset map after the
  branch is selected. The LOWER-only bounded parameter-subbox refinement
  in `hybrid_discrete_time_evolver.cpp` addresses uncertain **event
  guards and feasible branches**, not a free inclusion input: it retains
  the source numerical map and original constraints. Its optional
  midpoint is used only to *check* a feasible parameter source, not
  to bind a numerical input. Nevertheless hybrid LOWER remains partial
  and currently rejects inclusion dynamics explicitly.
- Validated hybrid discrete inclusions: full-box UPPER input-history
  propagation exists in `hybrid_discrete_time_step_image.cpp`, but
  public hybrid LOWER continues to reject inclusion source/target modes;
  no LOWER-specific input sampling was found there.
- Ordinary validated continuous equations: propagation uses configured
  numerical integration steps for both semantic paths, with distinct
  horizon/subdivision/termination rules. Ordinary continuous inclusions
  reject LOWER rather than selecting one free input.
- Hybrid validated continuous equations: existing LOWER has extra
  guard/event-admissibility checks and some paths call
  `apply_integration_step(..., Semantics::UPPER)` to construct a
  conservative structured tube. **Follow-up technical review needed**
  before treating those guards/tubes as a general LOWER contract; no
  explicit midpoint-input binding analogous to v7.12 was found. Hybrid
  continuous inclusions explicitly reject LOWER. This audit does not
  certify the full mathematical semantics of every continuous event.
- Approximate families: simulations may use caller-declared input
  policies, but that does not define the validated LOWER contract.

**Regression contract.** Tests now demonstrate the entire input box
survives LOWER at each tick: for x[n+1]=x[n]+u[n] with U=[-1,1] and
x[0]=0 the two-tick image must still cover [-2,2], with two free
independent INPUT axes and **zero artificial constraints**. The
multidimensional-input test retains its two original initial state
constraints and the bounds of the complete image; lower subdivision
cannot discard branches. Auxiliary-bearing LOWER retains the full
current and arrival input boxes, the a[n]=x[n]+u[n] correlation and
resumability, with no midpoint-equality constraints. Existing UPPER
expectations and legitimate rejection scenarios stay intact.

**Matrix:** ordinary/discrete/inclusion/Validated LOWER is a **new
set-propagation candidate** (partial pending the precise v7.13 CI).
This supersedes the v7.12 description of an existential selected-input
LOWER. Approximate 8/8; UPPER 6 full + 2 partial; LOWER 2 full +
3 partial candidates + 3 absent. This is a capability classification,
not a proof that numerical enclosures constitute an INNER image.
Only the three scoped inclusion files and the journal changed. All
five workflows (Ubuntu, macOS, Windows, Coverage, Sanitizers) must
pass on the exact corrective SHA before acceptance.

## v7.14 — Fix duplicate local variable in the v7.13 inclusion test

The macOS AppleClang 21 arm64 Debug job for v7.13
`9f23103e03fa37707e3dc76774b74d36d1faf42f` failed while compiling
`tests/test_discrete_time_inclusion_evolver.cpp` with `-Werror`:
`redefinition of 'lower_zero'` at line 287, first declared at line 261
within `test_rejection_is_explicit_even_at_zero_horizon`. This is a
**test-scoping/compilation error**, not an observed semantic failure of
the set-oriented LOWER propagation.

The second local, which holds the zero-horizon LOWER orbit for the
auxiliary-bearing inclusion, is renamed to `lower_auxiliary_zero`,
including its four expectations. The first `lower_zero` is unchanged.
Both tests remain enabled; no test assertions, production kernels,
semantics, input boxes, constraints, System, Integration or tutorials are
altered. In particular, the binding design rule from v7.13 holds:
`x[n+1]=f(x[n],u[n])` is applied to the full parametrized set with a
fresh full input box at every discrete tick under both UPPER and LOWER;
no lower-specific point or midpoint is selected.

**Verification:** require all five CI workflows (Ubuntu, macOS, Windows,
Coverage and Sanitizers) green on the precise v7.14 SHA, including
compilation with warnings as errors. No prior passing workflow can
stand in for this check. Matrix classification unchanged from v7.13.

## v7.15 — Hybrid discrete inclusion LOWER with whole-input progression and certified events

This version extends validated hybrid-discrete `orbit(..., Semantics::LOWER)`
to modes governed by bounded-input difference inclusions. **The dynamical
map and full set of input choices are identical to UPPER**; neither
input sampling nor midpoint/equality constraints are introduced.
The shared `HybridDiscreteTimeStepImage::apply` continues to compose
`f(x[n],u[n])` on the whole parameterized set and appends a fresh
independent full input box at each tick. For auxiliary outputs, the
current complete input box is available at zero horizon and after
resumption; destination input boxes are attached after a reset.

The shared discrete branch planner is retained. An inclusion source
can advance in LOWER only if the branch is unique and certified over
the **entire** source and input box. Uncertain or competing event
branches are rejected, not resolved by choosing a favourable input.
Existing bounded source-parameter refinement remains available for
equation modes only, and an equation may now reset into an inclusion
mode, preserving the target's full arrival input domain.

C++ regressions demonstrate full-input event-free inclusion progression
for two ticks, with independent inputs, zero synthetic constraints,
state bounds containing [-2,3], zero-horizon behavior and source
immutability. A forced equation-to-inclusion reset followed by full
input progression covers [9,12]. Auxiliary-bearing inclusion progression
exposes the full input box already at zero horizon, preserves
`a=x+u[n]` and covers [-3,7] at tick 2, with successful resumption
and no input-fixing constraints. A forced reset into an auxiliary-
bearing inclusion followed by progression covers [5,9]. A dedicated
urgent guard `a>=-2`, true for every possible input, certifies
`x'=x+a+3` over [2,6] while preserving the full input interval;
`a>=0`, uncertain across inputs, remains rejected. Previous UPPER
regressions and equation LOWER checks remain intact.

**Capability matrix candidate:** the discrete hybrid inclusion LOWER
cell advances from absent to partial (only globally certain branches),
so LOWER has 2 complete, 4 partial, 2 absent families, subject to CI.
No general INNER-image claim is made from outer numerical enclosures;
the semantic restriction on event choice is explicit. No changes to
System, Integration, UPPER stepping, approximate evolvers, or tutorials.

The base v7.14 Ubuntu and macOS workflows were green when this patch
was authored; the remaining v7.14 workflows were still in progress.
The entire v7.15 commit must independently pass Ubuntu, macOS,
Windows, Coverage and Sanitizers, including Debug and installed
tutorial, before acceptance.

## v7.16 — Repair duplicate local name in hybrid inclusion LOWER regression

The v7.15 `466c8ede` Ubuntu Clang 16 Debug build failed compiling
`tests/test_hybrid_discrete_time_evolver.cpp`:
`redefinition of 'bounds'` at line 435 (previous declaration line 414),
inside `test_public_reset_into_inclusion_with_auxiliary_and_continue`.
Both variables were in the same C++ function scope. This is a test source
**compilation** failure; no conclusion about the runtime or mathematics
of the newly introduced LOWER hybrid inclusion path follows from it.

Only the second `bounds` is renamed to `lower_bounds` and its two
existing assertions are updated to use the renamed variable. The
original UPPER bounds and all tests remain unchanged and enabled.
An additional source inspection counted `const auto` declarations
within all 19 hybrid evolver test methods and identified no other
duplicate method-local names.

No production kernel, numerical method, input domain, event planner,
guard classification, source/target System model, or contract changes.
In particular `Semantics::LOWER` continues to apply the same full
set-valued discrete inclusion map as `Semantics::UPPER`, without sampled
or fixed input choices. The v7.15 coverage of forced resets, guarded
auxiliaries, independent full input boxes and conservative rejection of
uncertain branches is preserved.

**Verification**: this is a corrective candidate; Ubuntu, macOS,
Windows, Coverage, and Sanitizers must all pass on the **exact v7.16
commit SHA**, including tests and installed tutorial, before acceptance.
Do not attribute the v7.15 failures solely to this issue until the
v7.16 test suite has run. Matrix entries are unchanged.

## v7.17 — Hybrid inclusion LOWER state-only guard refinement

**Precondition and unchanged contract.** Both validated UPPER and LOWER
advance discrete inclusions with the same map over the *whole* input
box at each integer tick. Independent current and historical INPUT
parameters and their ranges are retained, along with the input-
dependent auxiliary functions. Input sampling, midpoint selection,
equality constraints on inputs and any new numerical map are forbidden.

**Addressed incompleteness.** The v7.15/v7.16 hybrid inclusion LOWER
path admitted only a single uniformly certain branch on the entire
source. This rejected state-boundary guard cases that can be certified
on a strict subbox of the initial state *without narrowing any input*.
In `enqueue_parameter_halves()` a new optional variable-kind filter
allows `certified_lower_step()` to subdivide only coordinates tagged
`EnclosureVariableKind::INITIAL` when the source mode is an inclusion.
Every `INPUT` coordinate (current or historical), `PARAMETER` and
other protected kinds are excluded from this refinement. Equation
sources retain the exact previous axis-rotation behavior with no kind
filter. Candidate boxes preserve the complete input ranges, and the
same proof tests as before establish feasible source constraints and
exactly one certain branch over the *entire* candidate input domain.
The existing `constrain_lower_source()` only adds constraints on the
narrowed INITIAL state parameters; no input constraints are added.
The numerical `HybridDiscreteTimeStepImage::apply` remains unchanged
and is shared by UPPER and LOWER. This is a bounded sufficient search,
not a general decision procedure.

**Regression.** Added and registered
`test_lower_inclusion_refines_state_without_restricting_inputs`:
an auxiliary-bearing inclusion `a=x+u`, `x'=x+a`, with
`u in [-1,1]` and an urgent guard `x>=1/2` is initially
undecidable at its boundary. From initial state x in [1/2,1],
state-only refinement certifies the reset `x'=x+a+3`, preserves the
full input domain and contains both witnesses 3.75 and 5.75 for
x=7/8 and u=-1/+1; event history and source immutability are checked.
From x in [0,1/2], state-only refinement certifies event-free
progression, preserves both extreme input choices and contains -0.75
and 1.25 for x=1/8. Original full-input uncertain auxiliary guard
rejection, the already registered v7.15 cases, and UPPER tests remain.
All newly added test-local `const auto` names were checked for
duplicates across the 20 test methods before publication.

**Scope and acceptance.** Production change only in the hybrid LOWER
axis-selection helper and call site; one public C++ test method and
registration, plus this journal. No System, Integration, tutorial,
UPPER map, scheduling, input policy or mathematical semantics changes.
The v7.16 result is separately tracked and not inherited by this
commit. All five workflows (Ubuntu, macOS, Windows, Coverage and
Sanitizers) must pass on this v7.17 SHA before the new refinement is
validated. Matrix remains 8/8 approximate, UPPER 6+2 partial,
LOWER 2+4 partial+2 absent candidates; no cell is upgraded merely
because bounded refinement works.

## v7.18 — Canonical semantic contract pinned at the beginning of the design

This is a **documentation-only** change. The opening section of
`doc/system_families_design.md` is now the binding, non-negotiable
semantics contract for all validated evolvers and overrides historical
experimental descriptions in the journal. It explicitly distinguishes
rigorous outer enclosures of supersets (UPPER) from rigorous outer
enclosures of nonempty reachable subsets (LOWER), without mislabeling
LOWER as INNER. It fixes parity of the physical dynamics, full-input
set propagation at each discrete tick, proper use of continuous
integrators, LOWER subdivision/aggregation/provenance obligations and
hybrid event admissibility. It also makes exact-SHA five-workflow CI,
non-weakening semantic regression tests and explicit user authorization
for semantic changes mandatory.

**Important status warning:** the v7.17 code is unchanged by this
commit. Its test failure in `public_guarded_inclusion_upper_orbit`
and the older hybrid LOWER independent-subbox strategy remain open
issues requiring a separate, reviewed code correction. This version
does not pretend to repair that failure or advance the capability
matrix. Historical v7.05–v7.17 text that conflicts with the canonical
section is superseded as a design recommendation, but retained for
traceability. The earlier v6.53 witness audit remains consistent with
the lower-vs-inner distinction but must not be taken to authorize
arbitrary slicing of a later outer LOWER enclosure.

## v7.19 — Preserve LOWER reachable subsets by removing independent hybrid subbox selection

**Binding mathematical contract:** the canonical contract at the top
of this document governs every result. UPPER is a rigorous outer
enclosure of a superset of R; LOWER is a rigorous outer enclosure of
a nonempty subset of R, not an INNER approximation. A LOWER enclosure
must preserve at least one actually reachable point of the represented
subset. All temporal families apply the original dynamics.

**Confirmed error:** the v7.17 public hybrid regression
`public_guarded_inclusion_upper_orbit` failed because LOWER started
returning a selected state-refined child on an uncertain event rather
than rejecting the unresolved full-source branch. Child selection
had originated in v7.05 for equations and was extended to inclusions
in v7.17. From S_L nonempty and S_L subset E, splitting E into
two pieces does not prove the arbitrarily chosen child meets S_L.
This violates the established hybrid continuous LOWER principle.
The original negative regression is deliberately kept, not rewritten
to excuse independently selected subboxes.

**Production change:** eliminate the finite-depth subbox work queue,
axis splitting, per-box guard classification, synthetic source
parameter-box constraints and `CertifiedLowerStep` substitute source.
Both equation and inclusion LOWER now apply the SAME shared
`HybridDiscreteTimeStepImage::apply` numerical map as UPPER to the
ORIGINAL, UNSPLIT current enclosure, retaining complete independent
input boxes and source constraints. The initial source-feasibility
sufficient gate is checked once, not replaced by a per-tick search.
The shared branch planner must find one uniformly certain admissible
event/progression branch on the whole source. This is a conservative
SUFFICIENT supported transition policy, not the definition of LOWER
or a new physical dynamics. Event-free progression uses the normal
one-tick map without special event certification, while unresolved
branches are conservatively rejected.

**Semantic regression:** all 20 public test methods remain registered.
The existing uncertain-guard negative assertion remains unchanged.
Three formerly split-oriented tests now require rejection of an
ambiguous whole or boundary-inclusive initial enclosure, while
retaining positive mandatory checks for both reset and progression
on entire caller-provided guard-true/guard-false source boxes in
one and two dimensions. Auxiliaries, event identity/history, clock,
public image geometry, state constraints and initial-source
immutability remain under test. The inclusion case still checks
both extreme values of every relevant full [-1,1] INPUT interval:
no sampling, no midpoint equality, no invented subbox constraints.
Static preflight checked registration and duplicate const-auto
names, but is not a substitute for C++ compilation.

**Scope/status:** only hybrid LOWER orchestration, its public header,
public hybrid tests and this historical journal entry changed.
The top-of-document non-negotiable semantics are unmodified.
UPPER, System, Integration and numerical step-image code are untouched.
This is a candidate until the Ubuntu, macOS, Windows, Coverage and
Sanitizers workflows pass on this exact v7.19 commit. The capability
matrix does not advance on an unvalidated candidate.

## v7.20 — Safe permissive-event continuation in hybrid discrete validated LOWER

**Base:** v7.19 `4f261d55f46ce919dffbe714a0715df751a41a5b`.
Ubuntu, macOS, Windows, Coverage and Sanitizers all completed green
on that exact SHA. The canonical UPPER/LOWER definitions at the very
start of this document remain binding and unchanged.

**Mathematical reason for the minimal extension.** In a hybrid
system with a PERMISSIVE transition, an admissible progression is
still an actual behavior even if a reset could also be taken.
If the common discrete branch planner proves that the progression
invariant is TRUE on the complete current LOWER enclosure, then
applying the original unsplit map retains a genuinely reachable
subset. The LOWER enclosure can choose that admissible
uninterrupted subset; it does not claim to cover the optional
resets. This follows the existing validated hybrid-continuous
LOWER treatment, which explicitly ignores optional PERMISSIVE
jumps when the uninterrupted flow is valid.
No model-level priority between progression and reset is added:
UPPER still explores every allowed alternative.

**Implementation.** In `unique_certified_lower_branch`, keep the
existing shared, three-valued full-source planner. Accept a
certain PROGRESSION branch if every other listed alternative is
a PERMISSIVE transition; explicitly reject any competing
URGENT or unsupported event. When no such progression exists,
preserve the previous fallback of one uniquely certain
transition/progression. A possible guard never becomes TRUE or
FALSE just to produce an answer. Uncertain progression
invariants or urgent blockers are NOT ignored. Critically,
there is NO state/input/parameter splitting, no input selection,
and no changed numerical propagation kernel: both semantics
still apply the SAME full set-valued discrete map
`HybridDiscreteTimeStepImage::apply`, retaining all source
correlations and independent input boxes.

**Regression scope.** The existing public
`test_lower_unique_permissive_branch` retains its mandatory
reset on a forbidden progression invariant and its mandatory
progression with disabled guard; it additionally checks
progression when both PERMISSIVE reset and progression are
certain, when the optional guard is UNKNOWN but progression
is certain, and conservative rejection when both guard and
invariant are UNKNOWN. The existing
`test_lower_excludes_unproved_branches` retains zero-horizon,
input inclusion and invalid resume coverage, but its formerly
incorrect expectation of failure for a provably allowed
PERMISSIVE continuation is replaced by a positive semantic
assertion with the expected state image. An auxiliary-bearing
hybrid inclusion regression also checks optional guard
`a>=0` with safe invariant `a<=3`, intact full input
`u in [-1,1]`, no event, independent arrival input and
state image covering [-1,3]. Its otherwise identical URGENT
ambiguous guard still requires conservative rejection.
All 20 public test methods remain registered without duplicate
local `const auto` declarations.

**Scope/acceptance.** No changes to the automaton model,
numerical step image, UPPER planner, continuous hybrid
evolver, or the canonical UPPER/LOWER semantics. This is
a narrow supported case of LOWER hybrid transitions, not a
new semantics nor a capability-cell promotion. Full
Ubuntu, macOS, Windows, Coverage, Sanitizers CI must be
green on the exact v7.20 SHA before acceptance.

## v7.21 — Certain PERMISSIVE reset with undecidable progression invariant

**Base and acceptance state.** The v7.19 `4f261d55` had all five
Ubuntu, macOS, Windows, Coverage and Sanitizers workflows green.
v7.20 `4e460edc` was still in progress when v7.21 was prepared;
v7.20 CI results are not inherited and this version needs its own
five-workflow validation.

**Canonical proof obligation (unchanged).** Both UPPER and LOWER produce
rigorously outer numerical enclosures; LOWER encloses a nonempty
reachable subset. Neither changes the original discrete dynamical
map. A LOWER enclosure is never independently split, narrowed to a
favourable source child, or given a sampled input. Transition
admissibility is distinct from applying the map on the unsplit set.

**New narrowly supported hybrid case.** In the common planner a
PERMISSIVE transition with guard TRUE on the *entire* current source
is an admissible jump for the actual reachable subset, regardless
of whether the *separate progression invariant* is TRUE, FALSE or
UNKNOWN. When the progression invariant is UNKNOWN, v7.20 was
needlessly rejecting the valid reset because the planner included
a possible progression alongside the certain transition.
v7.21 accepts a certain PERMISSIVE transition if and only if
every other planned branch is merely PROGRESSION, and the target
mode has a supported equation/inclusion. It does not select
between multiple events or ignore any competing URGENT transition.
A certain progression with optional competing resets remains
supported by the v7.20 case. UNKNOWN guards remain UNKNOWN.
No priority is added to model semantics: choosing a proven
legal behavior yields a valid LOWER reachable subset, while
UPPER still considers all possible alternatives.

**Public regression.** Extended existing
`test_lower_unique_permissive_branch`: for x in [0,1],
`guard x>=-1` is uniformly true while independent progression
invariant `x<=1/2` is undecided; the reset x'=x+3 must keep the
complete x source, yield target mode and event provenance, and
cover [3,4] without synthetic constraints. Retained the
explicit negative case where BOTH the guard and the invariant
are undecidable. Extended the inclusion/auxiliary regression:
`a=x+u`, x in [0,1], full u in [-1,1],
`guard a>=-2` TRUE, invariant `a<=0` UNKNOWN.
The correlated reset `x'=x+a+3` must cover [2,6],
retain the entire u input and event history, and create no
artificial constraints. Existing v7.20 permissive progression,
urgent rejection, input-history and other regression checks
remain intact. Static source-level preflight found 20 registered
test methods and no duplicate method-local const-auto names.

**Scope.** Only one full-source branch admission case, public
header docs, two existing test methods and this journal entry
changed. No numerical kernel, integrator, model, UPPER
planner, test removals, source splitting or canonical
semantics change. All five CI workflows must be green on the
exact v7.21 commit before its behavior is accepted; no matrix
cell is promoted on an unvalidated candidate.

## v7.22 — Aggregate every uniformly certain hybrid LOWER branch without splitting the source

**Mathematical contract:** the immutable introductory section remains binding.
A LOWER result is a rigorous outer numerical enclosure of a nonempty
reachable subset, not an INNER image. If an unsplit source enclosure
contains such a subset and two hybrid branches are admissible
throughout that entire source, both propagate a genuine reachable
subset. Keeping both successor enclosures aggregates actual alternative
behaviors without independently splitting the numerical state/input
enclosure or fabricating an event priority. A branch with UNKNOWN
preconditions is not thereby declared certain.

**Production:** replace the former singled-out jump/progression
preference logic with the shared planner's complete list of
`certain` branches whose target dynamics are supported. At each tick,
validated LOWER applies the UNCHANGED `HybridDiscreteTimeStepImage::apply`
to the original complete source for every certain branch and preserves
each successor's location and event history. A frontier of enclosures
is advanced breadth-first. Reach/intermediate/final lists retain
all resulting certain branches. Uncertified branches are omitted;
if none reach the requested horizon, the call conservatively fails.
The existing maximum-node limit is applied per actual successor.
No new numerical integration method, free input selection, synthetic
parameter constraints, state splitting, urgency priority or changed
model dynamics is introduced.

**Public regression extensions (20 registered test methods):**
A PERMISSIVE guard and progress invariant both certainly true yield
BOTH the reset and uninterrupted evolution, with original state
intervals, separate event histories and two final enclosures.
A second two-tick model with generously true invariants produces
three final branches and five reachable nodes; a branching node
budget triggers as expected. Two simultaneously true URGENT events
yield two reset descendants with their individual event identities
and exact model-based bounds, and no false progression. A hybrid
inclusion `a=x+u`, `u in [-1,1]`, with optional guard `a>=-2`
and progression invariant `a<=3` yields two distinct successors:
full-input progression over [-1,3] and correlated reset over [2,6].
The original input coordinate range and arrival input binding
are checked separately, with no input-equality constraints.
Old rejections for genuinely uncertain transitions and all existing
UPPER, resume, invalid-source and other regressions remain.
Pre-commit static inspection found exactly 20 registered methods
and no duplicate local `const auto` names.

**Acceptance:** v7.20 was five-for-five green. v7.21 was green for
Ubuntu, macOS and Coverage but Windows and Sanitizers had not
finished at preparation time. This new behavior is a CANDIDATE
until all Ubuntu, macOS, Windows, Coverage, Sanitizers workflows
succeed on the exact v7.22 SHA. Do not promote any capability matrix
cell merely due to a new candidate. No changes to System, Integration,
the continuous evolvers, UPPER scheduler or canonical semantic rules.

## v7.23 — Structured validated LOWER for ordinary continuous inclusions, without changing dynamics

**Canonical mathematical contract:** UPPER encloses a superset of the
reachable set; LOWER rigorously encloses a NONEMPTY reachable subset,
not a geometric INNER subset. Both numerical enclosures are outer.
No model-specific LOWER right-hand side, input midpoint, or invented
physical input policy is allowed. The distinction for continuous time
is the supported validated integrator and the safe handling of
enclosures, not a different dynamics.

**Narrow implementation:** the existing selected-parametric
multifunction integration path already computes a rigorous
structured reach patch and structured endpoint patch, keeping
the independent free input coordinates per step. The ordinary
continuous evolver may now admit `Semantics::LOWER` only on that
selected structured integration path, using exactly the SAME
`compute_inclusion_step`/`compute_clocked_inclusion_step` and
`apply_integration_step` functions as UPPER. Explicit runtime
checks require both `has_reach_models()` and
`has_evolve_model()` before a LOWER step. The existing box-only
`FirstOrderMultifunctionIntegrator` remains UPPER-only: its
Cartesian box has no carried numerical trajectory-to-parameter
map and is explicitly prohibited by the shared integration layer
for LOWER. Constant physical-parameter LOWER remains unsupported
by the preexisting general parameter restriction; the new support
is currently for the selected structured solver with no such
physical parameters. The user chooses input-approximation
methods when constructing the integrator; these are its original
numerical approximation configuration, not a new unilateral
LOWER sampling policy.

**Existing unsound edge case corrected:** `_append_initial_set`
used to split the initial outer enclosure on radius without
checking semantics, even when LOWER held an existing witness
possibly absent from each arbitrary child. For LOWER it now
appends the complete original source directly; UPPER retains
the existing initial subdivision behavior unchanged. Subsequent
LOWER steps already avoid subdivision when radius is large.
This correction applies to both ordinary continuous equations
and inclusions, without touching their numerical kernels.

**Regression tests:** `test_structured_lower_inclusion_preserves_switching_inputs`
compares validated UPPER and LOWER on the SAME
`x'=u, y'=x` inclusion for two 1/16-length steps. It asserts
one final enclosure, four independent parameters (2 initial
state coordinates plus full u[0],u[1]), no synthetic constraints
and enclosure of the analytically reachable switching trajectory
`u[0]=+1,u[1]=-1`, `x(1/8)=0`,
`y(1/8)=1/256` under BOTH semantics.
`test_lower_preserves_unsplit_initial_enclosure` forces a
small max radius and enable_subdivisions, then checks that a
zero-horizon LOWER enclosure for x in [-1/4,1/4] retains
the WHOLE source as a single final enclosure. It also checks
that box-only integrator LOWER rejects at zero and positive
horizon rather than silently treating outer boxes as
trajectory certificates. All existing UPPER inclusion
tests and numerical bounds remain unchanged. Static preflight
confirmed all 13 test methods registered and no duplicate
method-local `const auto` declarations.

**Capability classification:** continuous/ordinary/inclusion/
Validated LOWER becomes a PARTIAL CANDIDATE (structured
selected-input method only), subject to CI and deeper
numerical-contract verification. The box integrator is not
promoted. This changes the candidate matrix LOWER to
2 complete + 5 partial + 1 absent (hybrid continuous
inclusion), with no changes in Approximate or UPPER.
The immutable rule at document start is unchanged. All
five CI workflows must pass on THIS exact SHA before
acceptance: Ubuntu, macOS, Windows, Coverage, Sanitizers.
The v7.22 CI has not been inherited.

## v7.24 — Correct the ordinary dynamics contract regression for structured inclusion LOWER

**Observed CI failure on v7.23** `3a546bfa`: in
`test_ordinary_dynamics_contract`, case
`inclusion_consumer_rejects_unsupported_arguments`, the old shared
assertion at line 104 expected that *both* explicit-time and
trajectory-constant-parameter inclusions throw on
`orbit(..., Semantics::LOWER)`. Ubuntu Clang 16 Debug (x86 and arm)
and macOS AppleClang Debug failed with “no exception thrown”.
The v7.23 capability extension intentionally made LOWER available
for a selected, structured continuous inclusion integrator, including
**explicit time**; the old negative assertion was not audited and
is therefore a **stale test expectation**, not evidence of a
production kernel defect. The second case, a physical
trajectory-constant parameter, is STILL intentionally rejected by
the evolver and must retain its negative coverage.

**Exact test correction, without weakening coverage:** keep the
existing independent model/interface dimensionality checks for
both variants. For `x'=t+u`, `u in [-1,1]`, source `x(0)=0`,
horizon `H=1/16`, verify that selected structured integration
under LOWER returns one final enclosure and a nonempty reach,
preserves at least one free input parameter and zero synthetic
constraints, and rigorously contains the analytical endpoint
interval `[-H+H^2/2, H+H^2/2]=[-31/512,33/512]`.
Also check that its zero-horizon LOWER orbit is valid. For
`x'=p+u` with `p in [0,1]`, retain mandatory
`std::invalid_argument` checks for LOWER at both positive
and zero horizon, since physical parameter LOWER is not yet
supported. The existing static type assertions, equation
negative cases and interval-function capability checks
remain unchanged. No production code is modified.

**Mathematical contract:** this is a test update matching the
already established *structured* selected integration capability;
it changes neither the model's dynamics nor the immutable
UPPER/LOWER definitions at the top of the document.
Validated numerical enclosures are outer under BOTH semantics;
LOWER represents a genuinely reachable nonempty subset, not a
geometric inner set. No input sampling or independent
enclosure splitting is introduced.

**Verification:** the prior v7.23 was CI-red in Debug. The
corrected v7.24 must pass ALL five workflows (Ubuntu, macOS,
Windows, Coverage, Sanitizers), including Debug tests,
on its own exact SHA before claiming validated support for
continuous inclusion LOWER. The capability matrix is still
2 LOWER complete + 5 partial candidates + 1 absent, and
must not be promoted solely by publishing this change.

## v7.25 — Structured hybrid continuous LOWER for unblocked inclusion flow

**Canonical semantics:** every numerical UPPER enclosure rigorously
contains a superset of the reached set; every numerical LOWER enclosure
rigorously contains a nonempty reached subset. Both numerical
enclosures are OUTER, not geometric INNER approximations. The source
is never split into independently claimed LOWER children, and the
physical differential inclusion and all admissible input boxes remain
unchanged.

**First supported hybrid-continuous inclusion LOWER case:** the
existing `SelectedParametricMultifunctionIntegrator` already
constructs a validated reach-tube patch and endpoint patch with the
independent admissible numerical inputs. `_evolve_inclusion`
now uses that SAME step for UPPER/LOWER. It is called with the actual
`Semantics`, checks that the structured reach and endpoint models
exist for LOWER, and applies the common
`EvolutionDetail::apply_integration_step` to the whole source
enclosure without sampling inputs or changing the dynamics.
The lower step takes the downward validated bound of the requested
time horizon; the UPPER step retains its former upward bound.

**Admissible hybrid behavior:** LOWER supports event-free modes and
modes whose transitions are exclusively PERMISSIVE. Choosing never
to take an optional jump yields a legal nonempty reachable subset,
including when that optional guard is UNKNOWN during continuous
progression. The ordinary shared `publish_outer_step` helper is
called with an EMPTY event map on this LOWER path: it propagates the
same unsplit time-flow enclosure, but cannot invent potential
UPPER guard-intersection jump candidates as LOWER witnesses.
UPPER continues passing all real transitions, unchanged.
Full-source initially certain optional jumps (if any) remain
separately governed by the established common hybrid starting-event
processing; no priority or source splitting is introduced.
URGENT/IMPACT/PROGRESS blockers and state INVARIANT constraints
remain explicitly unsupported by this new LOWER inclusion path
until witness-preserving event analysis is developed. Physical
trajectory-constant parameters are rejected consistently with
ordinary continuous inclusion LOWER. First-order box-only
integration is still UPPER-only, including zero horizon.
This is a conservative supported subset, NOT a new definition
of LOWER semantics.

**Public regressions:** two test methods added to
`tests/test_hybrid_inclusions.cpp`. The first runs
`x'=u, u in [-1,1]` from x in [-1/8,1/8] to time 1/8
with two semantic modes, checking same full [-1/4,1/4]
analytical bounds, correct time, original mode/event provenance
and at least four independent full INPUT domains under LOWER.
It checks that the old source remains unmodified and that
box-only LOWER rejects at both zero and positive horizons.
The second checks `x'=u` with a PERMISSIVE threshold
x>=1/32; LOWER advances without spuriously taking a jump and
still encloses [-1/8,1/8]. A corresponding URGENT guard
and a state INVARIANT are explicitly rejected, at zero and
positive horizons as applicable. No existing UPPER assertion
was deleted or weakened. Static preflight confirms 14
registered test methods and no duplicate method-local
`const auto` declarations.

**Matrix (candidate only):** Approximate 8/8, Validated UPPER 6
complete+2 partial, Validated LOWER 2 complete+6 partial+0 absent.
The hybrid continuous inclusion LOWER cell becomes PARTIAL
(candidate) for validated unblocked structured flows only,
subject to a full semantic and CI review. The five workflows
Ubuntu, macOS, Windows, Coverage and Sanitizers must all be
green on the exact v7.25 SHA before promotion. The canonical
contract at the beginning of this document is unchanged.
No changes to the numerical integrator implementation, UPPER
event processing, other evolvers, System, or tutorials.

## v7.26 — Hybrid inclusion LOWER under certified whole-tube invariants

**Source and legal semantics:** starts from v7.25 `923cce33`.
The immutable UPPER/LOWER contract at the document beginning
remains binding. Both semantics use the same actual differential
inclusion and the same selected structured numerical integrator.
LOWER encloses a nonempty reached subset, not an INNER image.
No sampling of inputs or independent LOWER subdivision is allowed.

**Sufficient validity proof:** an ordinary whole-source LOWER
witness transported by a structured validated inclusion reach
patch remains reachable in a hybrid mode carrying state
INVARIANT `g(x,t)<=0` whenever that guard is rigorously true
on the ENTIRE initial LOWER enclosure and on EVERY full reach
tube returned by the existing validated integrator. The
resulting trajectory family cannot leave the legal hybrid
state region. A mere nonempty intersection of an OUTER tube
with the invariant is not evidence of a reached LOWER subset:
that would risk discarding the original reachable witness.
This extends the existing v7.25 event-free/PERMISSIVE-only
case strictly by a check over the original sets, and introduces
no new flow equation, guard priority, integrator, or state
partition. The already-established hybrid continuous
invariant convention is `g<=0`.

**Production:** allow EventKind::INVARIANT together with
PERMISSIVE in inclusion LOWER modes; still reject all
URGENT, IMPACT and PROGRESS blockers at this stage.
Immediately after extracting transition predicates in
`_process_jumped_set`, and crucially BEFORE checking an
initial/terminal zero horizon or applying any invariant
restriction, verify that every state invariant has validated
upper bound <=0 over the ENTIRE source enclosure.
Otherwise reject with `std::invalid_argument`. In
`_evolve_inclusion`, construct a full `HybridEnclosure`
for each already computed structured reach tube, preserving
time/auxiliary/source arguments, and require each invariant
bound <=0 throughout each tube. If a bound is unknown or
positive, reject rather than creating a favourable numeric
subbox or choosing a special input. For certified whole
tubes the original common numerical map and
`publish_outer_step(..., {}, ...)` are unchanged; optional
jumps may still be ignored, and UPPER continues its
unmodified outer branch exploration.

**Public regression:** add one registered
`test_lower_inclusion_certified_whole_tube_invariants` in
`tests/test_hybrid_inclusions.cpp`. With `x'=u`,
full `u in [-1,1]`, initial x=-1, invariant x<=0
and H=1/8, all trajectories stay inside the invariant,
and LOWER must return one final source-mode enclosure
covering [-9/8,-7/8]. Test the 4 independent full input
intervals, unchanged event history and exact final time.
A boundary initial x=0 at zero horizon remains legal,
but must reject H=1/8 because positive inputs could
violate the invariant and no child can be selected.
Already-invalid x=1/4 is rejected even at horizon
zero. Existing UPPER invariant tests remain unchanged.
All 15 test methods are registered and there are no
duplicate `const auto` declarations in the static
source-level preflight.

**Status and matrix:** the hybrid continuous inclusion
LOWER cell remains PARTIAL. The v7.25 five-workflow CI
was still running when v7.26 was prepared; therefore
the v7.25 candidate is not presumed valid. v7.26 needs
its own green Ubuntu, macOS, Windows, Coverage and
Sanitizers workflows on its EXACT SHA, including Debug
and tutorial. No capability is promoted on publication.

## v7.27 — Inactive hybrid flow blockers preserve LOWER reachability

**Binding semantics:** the canonical contract at the start of this document
continues to govern. UPPER encloses a reachable superset; LOWER encloses a
nonempty genuinely reachable subset, with an OUTER numerical enclosure in
both cases. A numerical child of LOWER or a chosen input is never treated
as a free-standing witness. The original differential inclusion and
structured validated integration are IDENTICAL to the existing v7.25/v7.26
numerical path; no new integrator, input policy or event priority is used.

**Mathematical sufficiency:** for an uninterrupted flow, the mode's
INVARIANT `g<=0` must be validated on the ENTIRE source and ENTIRE
validated reach tube (already supported by v7.26). Equally, each blocking
URGENT, PROGRESS, or IMPACT guard must remain STRICTLY INACTIVE
(`upper_bound(guard)<0`) on the entire source and every reach tube.
The numerical outer tube covers all actual trajectories represented by
the original unsplit source, so a strictly negative guard can never be
reached during that step. No event can interrupt those trajectories.
PERMISSIVE jumps do not block and are freely ignored when LOWER chooses
uninterrupted progression. Conversely, a merely possible activation
must conservatively reject the progression: a guard-intersection of
an outer enclosure is not a LOWER crossing witness.

**Production implementation:** a single shared helper in
`src/hybrid_continuous_time_validated_evolver.cpp` evaluates these
sufficient predicate conditions. For PERMISSIVE, uninterrupted flow
is automatically allowed. For INVARIANT it requires upper<=0; for
blocking URGENT/PROGRESS/IMPACT it requires strictly upper<0. The
helper is applied BEFORE any initial source restriction, including
at zero time in `_process_jumped_set`, and to EVERY complete reach
tube in `_evolve_inclusion`. Missing/UNKNOWN guard evidence produces
`std::invalid_argument`, not an independently selected subbox.
The prior event-kind blanket rejection has been removed for these
provably safe cases only; UPPER's event planning is untouched.
The existing `publish_outer_step(..., {}, ...)` still publishes
only uninterrupted LOWER flow; this new extension does NOT claim
to implement LOWER first-crossing resets of an active blocker.

**New regression (16 registered hybrid inclusion tests):** for each of
URGENT, PROGRESS and IMPACT guarding x>=0, the differential inclusion
`x'=u, u in [-1,1]`, with initial x=-1, can legally flow for H=1/8;
LOWER must preserve the complete [-1,1] input boxes, event-free
provenance, correct final time, and the entire validated endpoint
range enclosing [-9/8,-7/8]. From x=-1/16, the initial guard is false
but the full trajectory tube can cross zero; LOWER must therefore
reject rather than claim the arbitrary uninterrupted behavior.
At x=0 the blocking event is active or undecidable already at time
zero; the presently conservative LOWER flow path rejects this case.
These are separate semantic cases, not modifications to earlier
v7.26 invariant tests. Static source inspection checked all 16
registered test methods and no duplicate function-wide `const auto`
names. The prior five workflows on v7.26 were 4/5 green with Windows
still running at preparation time; they do NOT validate this commit.

**Acceptance/matrix:** this is a new candidate partial-hybrid LOWER
capability, with no matrix upgrade. Require Ubuntu, macOS, Windows,
Coverage and Sanitizers green on the EXACT v7.27 SHA, and review the
existing continuous-hybrid event semantics before extending to active
URGENT/PROGRESS/IMPACT crossings. The binding contract and older
journal entries remain unchanged.

## v7.28 — Distinguish initial LOWER validity from future IMPACT tangency

**Observed failure and classification.** On v7.27 `7b08ca17`,
`test_hybrid_inclusion_impacts` reported 4 passed and 1 failed.
The failing case `tangency_with_structured_and_box_methods_and_both_backends`
expected a `std::invalid_argument` for LOWER at horizon **zero**,
for both the structured and box-only integrators. After v7.27,
that old blanket expectation is incorrect for a structured integrator:
the whole initial source has a strictly inactive IMPACT guard.
This is a stale capability/expectation mismatch, not sufficient
evidence of a kernel or numerical correctness bug. The box-only
integrator is still not supported under LOWER.

**Mathematics of the exact test model.** The source mode obeys
`x'=1, y'=u, u in [-1/8,1/8]`; its IMPACT guard is
`g(x,y)=y-x^2 >= 0`. Initially `x=-1/16,y=0`, so
`g=-1/256 < 0` on the **entire** source. At horizon zero,
no flow step or reset is needed, and the initial enclosure
rigorously encloses a genuinely reachable nonempty subset
regardless of later contact. For positive horizon H=1/8,
the family has possible crossings (including the tangency
at t=1/16 for u=0). The current LOWER path implements only
uninterrupted flow under fully inactive blockers and cannot
claim that whole flow is legal when the reach-tube guard is
not provably negative.

**Precise regression change.** Preserve all existing UPPER
tangency, target-mode reset, metadata, exact-contact time and
state assertions in all four structured/box x flat/tree
backend combinations. For the box-only method, continue to
require LOWER rejection even at zero horizon. For the
structured method, require one event-free LOWER final
enclosure at zero time, with original source location,
original exact state and no invented jump. Additionally
require conservative rejection of the structured LOWER
uninterrupted-flow attempt to H=1/8, where the impact
guard can become active. This latter check is a stated
implementation limitation, not a prohibition against
future certified LOWER impact transitions.

**Source and validation:** only
`tests/test_hybrid_inclusion_impacts.cpp` and this design
journal are changed. The immutable contract in the
document's opening section is unchanged. No production
kernel, original UPPER behavior, physical input choice,
test registration, matrix cell or guard interpretation
changes. Require Ubuntu/macOS/Windows/Coverage/Sanitizers
success on the exact v7.28 SHA, including Debug and installed
tutorial, before accepting the regression correction.

## v7.29 — Full-source URGENT immediate resets for continuous hybrid inclusion LOWER

**Semantic contract preserved:** validated LOWER returns rigorous OUTER
numerical enclosures of a NONEMPTY genuinely reached subset, never a
geometrical INNER approximation. The initial source remains the full
original enclosure with all independent input domains. Neither
priorities among events nor sampled input selections are invented.

**Distinguish a jump from the permission to progress:** the v7.27
full-tube blocker policy deliberately certifies uninterrupted
flow only while every URGENT/PROGRESS/IMPACT guard is strictly
inactive. However it incorrectly applied this FLOW precondition
also to CERTAIN URGENT IMMEDIATE RESETS. A nonempty source set with
uniformly positive URGENT guard guarantees an immediate jump,
whether or not ordinary progression is permitted. Treating
the two conditions as one was a conservative but unnecessary
capability limitation, not a licence to ignore the actual event
semantics.

**Production:** `_process_jumped_set` continues to check all
INVARIANTs and absent blockers over the full source, but now
also accepts an URGENT guard only when
`classify_initial_guard` certifies `FULL_SOURCE`
(strictly positive lower guard bound). Uncertain guards,
active PROGRESS, and active IMPACT are still rejected by
this sufficient source proof: a mere impact guard truth
does not prove its direction. The existing common
`_process_starting_events` already admits only uniformly
certain initial LOWER resets; its `guarded_clocked_lower`
safety check now also includes inclusion modes. Thus a
certain URGENT jump is emitted with full source and
event provenance while an ordinary no-jump flow is
suppressed if that blocker can activate. Multiple
certain URGENT events produce separate legal reset
branches; no arbitrary order is imposed. Numerical
integration and the later full-tube test remain
unchanged for all uninterrupted branches.

**Public regression:** new
`test_lower_inclusion_full_source_urgent_reset_without_flow`
with source `x'=u`, entire `u in [-1,1]`, initial
`x=1`, and TWO urgent transitions `x>=0` to distinct
stationary destinations with resets `x'=x+2`
and `x'=x+4`. LOWER with H=1/8 and max two
discrete transitions must reach BOTH target modes,
with one event each, source values 3 and 5, and
time exactly H. There must be NO spurious source
final. This checks explicit aggregation and does
not change the existing whole-tube input tests.
A guard only known on the boundary `x=0` remains
UNKNOWN and is rejected. An IMPACT starting with
a certainly true guard remains rejected without
a directional contact witness. UPPER and the
original v7.27/v7.28 tests remain untouched.
17 `test_hybrid_inclusions` methods are registered;
the new test passes static local-variable preflight.

**Acceptance:** v7.28 had green Ubuntu/macOS
CI at preparation time; other workflows were
still in progress. ALL five v7.29 workflows
Ubuntu, macOS, Windows, Coverage and Sanitizers
must be green on its exact SHA before claiming
this supported capability. The partial LOWER
cell and canonical opening design contract
remain unchanged.

## v7.30 — Whole-source PERMISSIVE inclusion LOWER resets beside PROGRESS

**Proof before implementation.** The original physical inclusion remains
`x'=f(x,u)` on its complete declared input domain. Let the incoming outer
enclosure `E_L` carry a nonempty reached subset `S_L`. If every source
invariant holds on all of `E_L` and a PERMISSIVE guard is strictly positive
on all of `E_L`, the original shared reset `r` is legal on all of `S_L`.
Consequently `r(S_L)` is nonempty, genuinely reached and enclosed by the
validated numerical reset of the whole source. PROGRESS can prohibit
continuous advancement; it does not prohibit this instantaneous optional
reset. No subbox, sampled input, new dynamics, event ordering or unknown
guard interpretation is needed. Existing state/input/time correlations
and source constraints must survive. The canonical opening contract is
unchanged.

**Small production change.** Before restricting a LOWER inclusion source,
`_process_jumped_set` detects whether an existing PERMISSIVE transition is
`FULL_SOURCE` according to the shared strict initial-guard classifier.
Only in that case may an active or undecidable PROGRESS predicate bypass
the uninterrupted-flow precondition. The existing common starting-event
handler emits all certified reset alternatives and suppresses flow if a
blocker is not certainly inactive. Invariant checks, unresolved URGENT
guards, IMPACT direction obligations, structured-integrator requirements
and physical-parameter limitations retain their previous treatment.
Neither the numerical integration kernel nor UPPER is changed.

**Public regressions.** Two new registered methods bring
`test_hybrid_inclusions` to 19 methods. In both flat and tree backends,
two certain PERMISSIVE resets beside `x>=0` PROGRESS retain both target
branches and no source final. The source is either `x=1` (active blocker)
or a carried LOWER enclosure produced by four full-input steps of `x'=u`,
`u in [-1,1]`, to `t=1/8` (undecidable blocker). With a one-jump budget,
both resets retain the original time, entire parameter domain, variable
kinds and event provenance. Numerical bounds check both the complete
reset image and its accuracy; source state/domain/constraints are not
mutated. Separate negative cases retain rejection for an uncertain
PERMISSIVE guard, an invariant not certified on the whole source,
an unresolved URGENT guard and an active IMPACT lacking direction.
Existing whole-tube, zero-horizon, urgent-reset and UPPER tests remain.

**Validation and formatting.** The three focused baseline tests passed on
v7.29. The new positive regression failed on the old production path with
the expected whole-source-continuation exception, then passed after the
change alongside both focused companion suites. As explicitly requested
by the user, `apply-format` and `check-format` run before committing;
formatting also normalizes accumulated drift in other root-repository
C++ files. No submodule source or gitlink is changed.
The full local Debug build with warnings as errors and enforced Python
tests passed all 49 CTest entries. Installation and both existing C++ and
Python tutorials against the installed library also passed. All 51
collateral C++ file changes exactly match clang-format applied to HEAD.

**Acceptance.** On 26 November 2026 all five functional v7.29 workflows
(Ubuntu, macOS, Windows, Coverage and Sanitizers) succeeded on exact SHA
`025d7b7f5fddb365cb901cb02f8dd306091de9b2`. The two separate Codecov
threshold checks failed; this is distinct from the successful Coverage
workflow. v7.30 remains a candidate until the same five functional
workflows succeed on its exact published SHA. The capability matrix
remains unchanged: Approximate 8 complete, UPPER 6 complete + 2 partial,
LOWER 2 complete + 6 partial. Python tutorial parity remains deferred to
the final phase; no work is added to the separate paving decoupling.

## v7.31 — Shared initial-reset policy beside PROGRESS

**Proof before implementation.** For a full original LOWER enclosure E
carrying nonempty reached S, whole-source invariants and a FULL_SOURCE
initial reset guard imply that the original reset r maps S to a nonempty
reached r(S) enclosed by the validated reset of E. This argument applies
to both policies that reset initially: URGENT and PERMISSIVE. PROGRESS
limits continuous flow, not an instantaneous certified reset. No event
priority, input selection or source subdivision is introduced.

**Implementation.** Extract the source proof into
`validate_lower_inclusion_initial_continuation`, keeping initial event
processing distinct from its sufficient continuation certificate. Use
the shared `resets_initially` transition policy instead of an ad hoc
PERMISSIVE-only condition. An active or unknown PROGRESS predicate can
therefore coexist with a FULL_SOURCE URGENT reset as well. Unproved
invariants, unresolved URGENT guards and IMPACT direction obligations
retain their existing rejection behavior. The common reset handler and
numerical integration remain unchanged; unsafe source flow is suppressed.

**Regression.** Expand the v7.30 positive regression over both initial
reset policies, retaining both flat/tree backends, active and uncertain
PROGRESS sources, two reset branches, complete carried input domains,
state/time/parameter correlations, event provenance, tight numerical
bounds and source immutability. The URGENT extension fails on v7.30 with
the expected continuation exception. Existing negative and architectural
tests remain registered without weakened expectations.

**Acceptance.** All five functional workflows passed on v7.30 exact SHA
`d069c087164a90cada89e648b8bab0129cb683d6`; the separate Codecov threshold
checks failed. v7.31 is a candidate until its five workflows pass on the
exact published SHA. The canonical contract and capability matrix remain
unchanged. Run apply-format/check-format and the full local suite before
publication.

**Local verification:** Debug build with warnings as errors and enforced Python tests, 49/49 CTest entries, both installed tutorials, apply-format and check-format all passed.

### Requested support checkpoint: native systems and public evolvers

This table records the current implementation scope, preserving the
existing conservative capability classifications. "Implemented" does not
extend a native model to undeclared time/parameter APIs or implicit DAEs.
"Partial" means that the stated subset works, not that the entire family
is absent. Validated LOWER remains an OUTER enclosure of a nonempty
reached subset under the immutable opening contract.

All named template evolvers in the table have `<ValidatedTag>`.

| Native family | Public validated evolver | UPPER | LOWER |
| --- | --- | --- | --- |
| Ordinary discrete equation | `DiscreteTimeEvolver` | Implemented | Implemented |
| Ordinary discrete bounded-input inclusion | `DiscreteTimeInclusionEvolver` (not a template) | Implemented | Partial: whole-input map, independent input history and correlated explicit outputs; no input sampling or source splitting |
| Ordinary continuous equation / explicit ODAE | `ContinuousTimeEvolver` | Implemented, including explicit time and fixed physical parameters | Implemented for the existing no-physical-parameter scope, including explicit time; fixed physical parameters still rejected |
| Ordinary continuous inclusion / explicit ODAI | `ContinuousTimeEvolver` | Implemented for both configured structured and box methods | Partial: structured selected-input method, full input domains and explicit time; no box-only method or fixed physical parameters |
| Hybrid discrete equation modes | `HybridDiscreteTimeEvolver` | Implemented | Partial: certified feasible initial source and whole-source certain progression/reset branches; aggregates certain alternatives, rejects unresolved continuations |
| Hybrid discrete bounded-input inclusion modes | `HybridDiscreteTimeEvolver` | Partial: full input propagation, correlated auxiliary guards/resets and inclusion targets are implemented | Partial: same full-input numerical map and certain-branch policy, preserving state/input/auxiliary provenance |
| Hybrid continuous equation modes / explicit ODAE | `HybridContinuousTimeEvolver` | Implemented, including time/parameters | Partial: time/parameters, certified immediate resets and supported certified crossings; not every geometry or competing-event case |
| Hybrid continuous inclusion modes / explicit ODAI | `HybridContinuousTimeEvolver` | Partial: both methods, time/parameters, invariant restrictions and conservative event/contact candidates | Partial: structured method, time, whole-source/whole-tube certified invariants and inactive blockers, certain immediate URGENT/PERMISSIVE resets; fixed physical parameters and general inclusion crossings remain unsupported |

**v7.31 candidate delta only:** the certain-immediate-reset case in the
last row now admits PROGRESS active or unknown beside either URGENT or
PERMISSIVE, according to the shared initial-reset policy. The PERMISSIVE
case is already covered by green v7.30 CI. The URGENT extension requires
all five new workflows on the v7.31 SHA. No whole cell is promoted.

**Time and parameters.** Continuous explicit time is an external clock,
not an extra user state with derivative one. Discrete models currently
declare state/auxiliary assignments and, for inclusions, bounded inputs;
they do not yet expose the continuous model's explicit clock and named
trajectory-constant physical-parameter declarations. Integer scheduler
indices do not imply general symbolic nonautonomous discrete dynamics.
Unifying ordinary continuous LOWER parameter support with the existing
hybrid equation path, then structured inclusion paths, is the next audit.

**Approximate evolvers.** All eight native equation/inclusion families
have an implemented approximate path (`DiscreteTimeEvolver`,
`ContinuousTimeEvolver`, `HybridDiscreteTimeEvolver` and
`HybridContinuousTimeEvolver`, each with `<ApproximateTag>`). Inclusion
inputs follow an explicit approximation/sampling policy. UPPER/LOWER are
validated semantics and do not apply to those trajectories. General
implicit DAEs/DAIs and final Python tutorial/binding parity are not
implied by this matrix. Paving decoupling stays outside this work.

## v7.32 — Persistent physical parameters in continuous validated LOWER

**Proof recorded before implementation.** The physical dynamics remain
`x'=f(x,t,p,u)`, with the same trajectory-constant tuple `p` in its full
declared domain P and every admissible input function u. Let the complete
outer LOWER source E carry a nonempty reached subset S. An existing named
physical-parameter binding is retained; an unbound initial source is
extended by all of the nonempty P using the same preparation as UPPER.
The structured integration step encloses the original physical flow for
the entire source and full input domains. Its argument map composes the
original correlated state and constant-parameter functions, so p keeps
the same value along each trajectory, across steps and resumed orbits.
Fresh input coordinates remain independent of the persistent physical
coordinates. The exact flow image of S (paired with P when initially
unbound) is nonempty, reached and enclosed by the numerical image of E.
The same argument covers equation auxiliary functions and pre-bound
parameter maps whose uniform uncertainty has persistent coordinates.
No source subdivision, input sampling, new feasibility gate or change
to the physical dynamics is required.

For hybrid modes, every original whole-source/whole-tube invariant and
blocking-event condition remains mandatory. A full-source certified
initial reset transports the same correlated physical tuple; the common
mode preparation checks compatibility with the destination declaration.
Thus allowing parameters does not license uncertain guard intersections,
IMPACT resets without direction or general inclusion crossings.

**Small implementation change.** Remove the blanket parameter/LOWER
rejection in ordinary continuous orbit preparation, covering equations
and structured inclusions, and the corresponding rejection in hybrid
inclusion mode validation. Reuse all integration, clock lifting,
parameter binding, auxiliary composition and reset kernels unchanged.
Existing parameter name/order/domain checks and preconditioning
exclusions remain. Box-only inclusion methods stay UPPER-only.

**Regressions, retaining original UPPER assertions.**

- Run the existing ordinary parameter/auxiliary analytic checks under
  both semantics, with two physical parameters, autonomous/explicit-time
  fields, all three initial-set entry forms and nonzero clock origin.
  Existing static coordinate/reconditioning checks remain enabled.
- Run the existing timed and autonomous inclusion extremal bounds,
  input-axis counts and physical-parameter sensitivity under LOWER too.
- New two-state ordinary inclusion regression integrates
  `x'=p+u+t, y'=x` (also without t), `p in [2,3]`, `u in [-1,1]`.
  Two resumed periods retain one physical axis and four independent full
  input axes, including a +1/-1 switching trajectory at both physical
  parameter endpoints. Analytic x/y/time values, zero-horizon resume and
  source immutability are checked. No internal time state is exposed.
- Extend the persistent-parameter rounding-error regression to LOWER
  over successive orbits, and retain explicit name/domain mismatch
  rejection for LOWER as well as UPPER.
- Extend hybrid inclusion method gates over both semantics and both
  backends, with/without explicit time and physical parameters. Structured
  LOWER works at zero/positive horizons; box-only LOWER still throws.
- New hybrid whole-tube regression checks invariants and PROGRESS guards
  depending on p, with/without explicit time and in flat/tree backends.
  An inactive blocker permits the full-input numerical image; a possible
  later crossing is rejected while its legal zero horizon is retained.
- New hybrid reset/continuation regression checks both URGENT and
  PERMISSIVE with joint parameter/time guards, active PROGRESS, a reset
  `x+=2p`, and a parameterized inclusion destination. It retains the
  original tuple, full independent input domains, event provenance,
  nonzero absolute time and tight correlated numerical bounds. A guard
  uncertain over p is rejected without selecting a parameter subbox.

The expanded existing tests fail against v7.31 with the two expected
parameter/LOWER admission exceptions; this is a deliberate capability
extension. No architectural test is weakened to conceal a failure.

### Support checkpoint after v7.32 (candidate until its exact CI passes)

Only the following parameter limits in the v7.31 matrix change:

| Validated evolver / family | UPPER | LOWER with physical parameters |
| --- | --- | --- |
| Ordinary continuous equation / explicit ODAE | Existing support unchanged | Added: autonomous and explicit time, persistent named tuples and correlated auxiliary outputs |
| Ordinary continuous structured inclusion / explicit ODAI | Existing support unchanged | Added: autonomous and explicit time, complete input domains and persistent tuples across resume |
| Hybrid continuous structured inclusion / explicit ODAI | Existing partial support unchanged | Added within the existing certified-flow/reset subset, including parameter/time predicates and compatible inclusion targets |

Continuous inclusion LOWER remains partial: the box-only method is not
promoted, and general hybrid crossing/contact reasoning remains open.
Hybrid equation LOWER, both discrete families and Approximate are
unchanged. The conservative family-level counts stay Approximate 8
complete, UPPER 6 complete + 2 partial, LOWER 2 complete + 6 partial.
This extension must pass Ubuntu, macOS, Windows, Coverage and Sanitizers
on the exact new SHA before its parameter capability is accepted.
The immutable opening contract, deferred Python tutorial parity and
separate paving decoupling remain unchanged.

**Architectural-test migration and verification.** The first full run
identified two stale capability expectations: the ordinary-dynamics C++
contract and Python consumer test still required parameterized LOWER to
throw. Replace those obsolete rejections with positive semantic checks,
retaining the original unsupported box/preconditioning cases. The C++
equation check validates `x(H)=p*(1-exp(-H))`; the inclusion check retains
the complete `[-1/16,1/8]` image at H=1/16. It checks that the only source
constraint is the original physical domain `0<=p<=1`, rather than an
input-selection equality. Python verifies successful structured LOWER,
resume/zero-horizon calls and the box-only rejection through the existing
bindings. No new bindings or tutorial features are introduced here.

On 27 November 2026 the final local Debug build with warnings as errors,
all 49 CTest entries including Python, installation, both installed
tutorials, apply-format and check-format passed. The six original focused
suites passed too. No submodule source/gitlink, numerical integration
kernel, UPPER behavior or canonical semantic contract is changed.
All five v7.31 functional workflows now succeeded on exact SHA
`87da52ae630e49981fa39b89a5aa66cafe9520b0`; its separate Codecov threshold
checks failed. Those results do not replace v7.32's own exact-SHA gate.

## v7.33 — Shared certified URGENT crossings for structured inclusions

**Proof recorded before changing the evolver.** Keep the original
`x'=f(x,t,p,u)`, every admissible input function, persistent p, and the
complete outer LOWER source E carrying a nonempty reached S. One complete
validated reach tube and its endpoint enclosure cover the uninterrupted
physical trajectories of that entire source/input family. If a continuous
URGENT guard g is strictly negative over E and nonnegative over the entire
endpoint, continuity gives a first zero within the step for each actual
trajectory represented by S. The complete endpoint time must also be
certified at or before the requested final horizon. Every invariant must hold on the whole tube,
and every other blocking guard must be strictly inactive on that tube.
This makes the first zero a legal hybrid reset point. The equality slice
`g=0` of the complete numerical tube contains the nonempty genuinely
reached crossing subset, and the common reset encloses its image.

This is the existing sufficient equation proof. It neither samples nor
freezes u, selects a source child, orders competing events, nor assumes
monotonicity or a Lie-derivative sign. Reach and endpoint approximations
can differ numerically; both must cover the entire original input family.
If any required certificate is unknown, no LOWER crossing is published.
Multiple reach tubes and IMPACT direction are outside this extension.

**Common handler.** Extract the existing clocked/parameterized equation
crossing proof and reset publication into
`EvolutionDetail::try_publish_lower_urgent_crossing` in the private common
hybrid outer-flow module. Equation orchestration calls the same handler;
structured inclusion orchestration tries it before requiring an entirely
uninterrupted tube. The numerical step, input approximation configuration,
UPPER path and initial source policy remain. Uncertified inclusion flows
retain their existing rejection behavior.
This separates crossing certification from integration and publication
without duplicating an equation/inclusion event engine. A failed attempt
publishes nothing; it never falls back to a possible UPPER jump.

**Public regression.** Three additional methods bring
`test_hybrid_inclusions` to 22 registered methods. The positive regression
uses `x'=1+u` or `x'=p+t+u`, `u in [-1/8,1/8]`, with/without explicit
time and `p in [1,2]`, under flat and tree backends. It requires an actual
URGENT boundary reset `x^+=x+1`, one event, the target mode, a time within
the certifying step, complete input/physical domains and an immutable
source. A joint guard `x+t>=p/64` retains the analytically reached witness
`p=1+1/256, u=0, t=1/128, x^+=1+t+t^2`; this is a test witness, not a
production sampling policy. The original code rejects the new positive
case with the expected whole-tube blocker exception.

Negative regressions retain rejection when the full endpoint fails to
prove a crossing (`x'=u`, `u in [-1,1]`), or another URGENT, PROGRESS,
INVARIANT or IMPACT condition cannot be excluded on the whole tube.
No competing event is prioritized. Existing equation crossings,
closed-endpoint contacts, inclusion tangencies, initial resets and
architectural tests remain enabled. General inclusion crossings and
IMPACT contacts are still partial; no whole matrix cell is promoted.
Acceptance requires all five functional workflows on the exact new SHA,
including the existing installed tutorials. The canonical opening
contract and deferred Python tutorial parity remain unchanged.

**Horizon audit and regression.** A certificate within a numerical step
is insufficient when source arrival times are uncertain: the step is
bounded from the earliest source time, while the actual LOWER subset
could arrive later and first cross beyond the requested horizon. The
shared handler now requires the whole endpoint absolute-time range to
lie within the final horizon before publishing a crossing. This is a
sufficient horizon certificate, not selection of early source times.
The new regression exercises both equations and inclusions, flat/tree,
with an outer arrival-time range [0,1/32], crossing x=1/32 and a short
horizon 3/64. It fails against the prior equation publication (one final
was returned without the necessary universal horizon proof). A larger
horizon 1/8 admits the same uncertain arrival source. Existing equation
closed-endpoint regression coverage is retained. This is a production
certificate omission exposed by a new test, not a reason to weaken an
existing expectation or redefine LOWER.

**Verification checkpoint (27 November 2026).** The final Debug build with
warnings as errors, all 49 CTest entries including Python, installation,
both installed tutorials, apply-format and check-format passed. All
22 inclusion methods remain registered. Existing equation crossings and
closed-endpoint tests pass after extraction and horizon certification.
The v7.32 parameter support is now accepted in its stated scope: Ubuntu,
macOS, Windows, Coverage, Sanitizers and Format all succeeded on exact SHA
`0201d8442a69e24ad4532635423587329a852f97`. v7.33 remains a new candidate
until its own five functional workflows succeed. No whole family cell is
promoted; numerical kernels, free input domains, physical dynamics and
the canonical contract remain unchanged.

## v7.34 — Certified PERMISSIVE endpoint resets beside uninterrupted LOWER flow

**Proof recorded before implementation.** Keep the physical dynamics
`x'=f(x,t,p,u)`, constant interval parameters p, every admissible input
function u and the complete original outer LOWER enclosure E carrying
nonempty reached S. Whole-tube invariants and strictly inactive
URGENT/PROGRESS/IMPACT blockers prove legal uninterrupted evolution,
including a nonempty actual endpoint subset S_h inside E_h. If a
PERMISSIVE guard is FULL_SOURCE (strictly positive) over all of E_h,
the original reset r is enabled on all of S_h. Require the whole endpoint
time range within the final horizon and every target invariant on the
whole reset image r(E_h). Then r(S_h) is nonempty, reached and enclosed
by that image. Publish this branch in addition to the no-jump behavior.
An unknown guard or target invariant yields no optional reset and never
removes the certified no-jump branch.

The configured numerical endpoint is an admissible PERMISSIVE event
time, not an input sampling policy. Independent full input domains,+persistent parameter maps, state/time/auxiliary correlations and source
constraints survive. Each reset copies the same complete endpoint; no
event priority, source child, target-invariant slice or chained same-step
reset is selected. The shared continuous prefix is credited once.

**Separate common handlers.** Move the unchanged uninterrupted-flow
predicate to `lower_flow_predicate_certified` in the private outer-flow
module and centralize its application to every whole reach tube in
`lower_flow_step_certified`. Initial inclusion checks and both ordinary
outer numerical producers now use the common policy. Preserve the
existing failure distinction: unsupported inclusion flow throws, while
the clocked/parameterized equation path reports no continuation.

After the no-jump result is assembled, the independent
`publish_lower_permissive_endpoint_resets` handler appends all endpoint
resets certified by the shared initial-guard classifier. It validates
target invariants over the complete reset image before publication,
without restricting the image to a favorable target subset. It retains
the no-jump contributions and continuous progress credit, and clears a
terminal stop reason if reset successors are added at the closed horizon.
The URGENT crossing handler stays distinct because its certified jump
replaces ordinary source continuation. UPPER and integration kernels
remain unchanged.

**Regressions.** Four additional public methods bring
`test_hybrid_inclusions` to 26 registered methods. Tests cover both flat
and tree backends, ordinary outer equation and structured inclusion
producers, autonomous fields, explicit time and constant interval
parameters. At one step's closed horizon, two certain PERMISSIVE guards
must produce both reset destinations plus the no-jump source final.
Analytic full-input bounds, original parameter domains/kinds, reset
images, auxiliary values, source immutability, time and event provenance
are checked. The original v7.33 path fails the new three-branch check
because it publishes only the no-jump branch; the negative cases already
pass before implementation.

Additional controls retain the no-jump result with an unknown endpoint
guard or unknown/violated target invariant; a reset at a late endpoint
cannot be inferred from an uncertain arrival-time enclosure beyond the
requested horizon. A later optional endpoint reset does not license
flow through an earlier PROGRESS blocker. A two-step regression checks
continued source flow and resets at successive endpoints, including
target continuation to the final horizon, independent input-history
counts, dwell times and one event per reset branch. All earlier tests,
including URGENT endpoint and horizon regressions, remain enabled.

**Capability boundary.** This extends the existing partial hybrid LOWER
support only when the entire prefix and endpoint reset are certified.
It does not add arbitrary within-step PERMISSIVE localization, IMPACT
direction, or independent input/source selection. No whole matrix cell
is promoted. All five functional workflows must pass on the exact new
SHA before acceptance. The canonical opening contract, deferred Python
tutorial parity and separate paving decoupling remain unchanged.

**Verification and next priority (27 November 2026).** The original
three-branch regression failed on v7.33 and all four new public methods
now pass. The final Debug build with warnings as errors, all 49 CTest
entries including Python, installation, both installed tutorials,
apply-format and check-format passed. The seven producer forms include
joint state/time/parameter endpoint guards and correlated target outputs.
All 26 inclusion test methods remain registered. Ubuntu, macOS, Windows,
Coverage, Sanitizers and Format all succeeded for v7.33 exact SHA
`01eae7f98bd611d442f57cf64459adee37ad860e`; this does not substitute for
the new v7.34 five-workflow gate.

At the user's request, after completing this LOWER increment the next
development priority is closing validated UPPER support. Audit the two
remaining partial native inclusion cells (hybrid discrete and hybrid
continuous), distinguish stale capability descriptions from real gaps,
retain architectural tests and add precise coverage regressions before
promoting a cell. General LOWER completion and final Python tutorial
parity are deferred while this UPPER priority is addressed.

## v7.35 — Complete native UPPER coverage audit and discrete initial outputs

**Priority and scope.** After v7.34's certified LOWER endpoint resets,
the requested priority is closing the two remaining native UPPER cells:
hybrid discrete inclusion and hybrid continuous inclusion. The audited
models are the native explicit update/auxiliary equations, full bounded
input inclusions and declared continuous interval-Taylor fields. This
does not introduce general implicit DAE/DAI solving, a symbolic discrete
clock or discrete named trajectory-constant parameter declarations.
Continuous explicit time and trajectory-constant interval parameters
retain their existing native contracts. Numerical success requires a
valid configured solver/model combination and its validity domain.
In particular, the structured inclusion path requires a free numerical
input when using an input-function representation; choose the box method
otherwise. Preconditioned integration is not supported for inclusion
approximations or time/parameter-dependent hybrid fields. These explicit
configuration limits remain errors, not missing event branches.

**Coverage proof recorded before implementation.** Preserve the opening
contract `R subset S_U subset E_U` for the complete requested finite
time/index and jump horizon. For discrete evolution, prepare the full
current input box and its correlated algebraic outputs. Retain every
possibly enabled transition and every possibly admissible progression;
UNKNOWN never discards a branch. Each branch copies the same source,
restricts only necessary closed pre-state predicates, and applies the
original update/reset. Introduce independent full destination inputs
only after this image. Induction on the integer index and complete
breadth-first fan-out covers every real branch, with parent/event
provenance retained. Project every horizon leaf; keep blocked leaves
separate. Node exhaustion throws and cannot publish an incomplete UPPER.

For continuous evolution, each validated numerical tube and endpoint
encloses the complete original flow family. Every real uninterrupted
pre-event point satisfies the source invariants and inactive blockers,
so their necessary closed restrictions retain it. Every actual first
URGENT/IMPACT contact belongs to the retained boundary; every legal
PERMISSIVE jump belongs to the retained guard-active tube slice.
Enumerating every possible candidate and applying the original
correlated reset therefore covers all actual alternatives. Spurious
contacts are permitted by UPPER; first-entry or crossing-direction
certainty is not required to retain an outer candidate. Initial events,
closed horizon restriction, compatible constant-parameter bindings and
both backend fan-out policies preserve the induction. The documented
HybridTime stopping rule remains reaching either the time or jump bound.
Numerical/resource failure aborts instead of returning partial coverage.

**Audit result and concrete correction.** The historical partial labels
for the two inclusion cells outlived the removal of their capability
gates (v6.97 and later). The present declared-model event/target paths
already implement the UPPER policy above. Their composed regressions
pass without changing event geometry, integration or physical dynamics.
One genuine initial-output gap was found in the discrete equation paths:
ordinary state-only seeds were rejected when the equation declared
auxiliaries, while hybrid zero-horizon roots omitted those outputs.
New tests reproduce both failures on the previous implementation.

The ordinary discrete enclosure factory and orbit entry now compose
the declared auxiliary map with the entire initial state. The hybrid
discrete UPPER scheduler does the same before publishing its root,
including zero horizon and resumed enclosures. Work on copies, retain
state-space validation and reject incompatible nonempty output labels.
This is a validated outer image of the declared algebraic outputs, not
a constraint selecting favorable initial states. The ordinary shared
entry also prepares LOWER outputs, without changing its propagation or
event policy. Correct stale deterministic-only and not-yet-public
descriptions of the hybrid discrete scheduler/evolver.

**Composed regressions.** Five new registered methods bring the ordinary
discrete suite to 6 methods, hybrid discrete to 22 and hybrid inclusion
to 28. The discrete inclusion regression follows two ticks through
correlated auxiliary guards and competing URGENT/PERMISSIVE resets into
inclusion targets: all nine nodes and five distinct horizon branches,
analytic full-input state maps/bounds, fresh full arrival inputs,
auxiliary residuals, time/dwell/event/parent provenance, orbit projection,
source immutability and hard node overflow. Initial-output tests cover
factory, raw enclosure and constrained-set entry, zero horizon, resume,
correlated values and incompatible labels in ordinary/hybrid equations.

Continuous regressions retain competing URGENT and optional PERMISSIVE
destinations with explicit time, constant interval parameters, source
and target invariants, correlated target outputs, and equation/inclusion
targets. A separate interval-Taylor-field regression retains urgent and
optional destinations. Both run structured/box methods and flat/tree
backends, verify analytic reached examples, absolute horizon and event
provenance. These examples guard the structural coverage proof; they do
not replace it with sampling. All earlier architectural tests, closed
horizon and IMPACT/tangency/competing-blocker tests remain mandatory.

**Current support checkpoint and acceptance gate.** All eight native
approximate paths remain implemented. For validated UPPER the candidate
matrix below closes the two historical partial cells only after Ubuntu,
macOS, Windows, Coverage and Sanitizers all succeed on this exact new
commit. Until then retain the accepted 6 complete + 2 partial count.
Validated LOWER remains 2 complete + 6 partial; no LOWER cell is promoted.
All template evolvers below use ValidatedTag unless stated otherwise.

| Native family | Public validated evolver | UPPER after exact-SHA acceptance | LOWER |
| --- | --- | --- | --- |
| Ordinary discrete equation | DiscreteTimeEvolver | Implemented | Implemented |
| Ordinary discrete bounded-input inclusion | DiscreteTimeInclusionEvolver (not a template) | Implemented | Partial: whole-input map, independent history and correlated outputs |
| Ordinary continuous equation / explicit ODAE | ContinuousTimeEvolver | Implemented, including time/constant interval parameters | Implemented, including time/constant interval parameters |
| Ordinary continuous inclusion / explicit ODAI | ContinuousTimeEvolver | Implemented, structured and box methods | Partial: structured method, full inputs, time/constant interval parameters |
| Hybrid discrete equation modes | HybridDiscreteTimeEvolver | Implemented | Partial: certified whole-source progression/reset alternatives |
| Hybrid discrete bounded-input inclusion modes | HybridDiscreteTimeEvolver | Implemented: full branching, auxiliary guards/resets and equation/inclusion targets | Partial: full-input maps and certified whole-source alternatives |
| Hybrid continuous equation modes / explicit ODAE | HybridContinuousTimeEvolver | Implemented, including time/constant interval parameters | Partial: certified flow, immediate/endpoint resets and supported crossings |
| Hybrid continuous inclusion modes / explicit ODAI | HybridContinuousTimeEvolver | Implemented: both methods, time/constant interval parameters, all native event/contact candidates and targets | Partial: structured certified flow, immediate/endpoint resets and supported URGENT crossings; general crossings/IMPACT remain open |

**Validation checkpoint (27 November 2026).** The five new methods and
existing impact suite pass after the initial-output correction. The
complete Debug/warnings-as-errors build, all 49 CTest entries including
Python, installation, installed C++/Python tutorials and
apply-format/check-format passed locally. The five exact-SHA functional
workflows are required before matrix acceptance. The preceding v7.34 is accepted:
Ubuntu, macOS, Windows, Coverage, Sanitizers and Format all succeeded on
`53e21a32e249cf1d825afb864918d56d7b329759`. The canonical opening section
is unchanged. General LOWER completion and final Python tutorial parity
stay deferred during UPPER closure; paving remains out of scope.
