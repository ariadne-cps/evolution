# Hybrid transition handlers

The hybrid evolver separates numerical capabilities, event policy, crossing
analysis, geometric guard operations and orchestration. The broader
[four-family roadmap](system_families_design.md) places this work before the
remaining continuous capabilities, discrete hybrid systems and simulators.
Event policies, guard handlers and crossing analysis are separate private
components. Public crossing data distinguish geometry from typed numerical
evidence while retaining the numerical rules and constraint order of the
equation path.

## Responsibilities

| Component | Owns | Does not own |
| --- | --- | --- |
| Capability validation | Accepted system/solver/semantics combinations and coordinate layouts | Crossing geometry or workload scheduling |
| Crossing analysis | Validated bounds, direction, curvature, crossing and critical-time information | Reset execution or publication |
| Event policy | Blocking, activation, initial-event treatment and impact direction | Numerical root solving |
| Geometric handlers | Flow restrictions and jump candidates with their elapsed-time maps | Other events, resets, horizons or backends |
| Orchestrator | Initial events, timing, competing guards, horizon, reset and local step results | Detailed geometric formulas |

Handlers are private functions selected by explicit policies. There is no
class for every product of dynamics category, semantics, event and crossing.
The evolver remains the sole owner of orchestration; handlers create no Orbit,
workload or backend and publish no results themselves. Implementations remain
in `.cpp` files.

## Event policy

These entries preserve the existing algorithm, including its separate initial
event phase. The initial impact treatment retains the current assumption that
an impact is not processed immediately after another event; revisiting that
assumption requires a separate numerical change and regression.

| Event | Blocks flow | Can reset during flow | Initial activation | Initial reset | Extra direction condition |
| --- | --- | --- | --- | --- | --- |
| INVARIANT | yes | no | separate invariant restriction | no | no |
| PROGRESS | yes | no | restrict flow if active | no | no |
| PERMISSIVE | no | yes | inspect guard | yes | no |
| URGENT | yes | yes | inspect guard and restrict flow | yes | no |
| IMPACT | yes | yes | deferred by the existing algorithm | no | guard Lie derivative |

`is_blocking` and `is_activating` consume this same policy. Optional jumps do
not remove the continuous branch. Initial urgent coverage remains certified
only by the existing strict whole-source activation test. A sole successor
does not prove full-source coverage.

## Two distinct geometric operations

`apply_flow_guard` restricts an existing reach, endpoint or jump branch to the
part allowed by a blocking event. `make_jump_candidate` constructs a possible
jump before reset. They share crossing information but have different
mathematical conditions: a flow restriction uses the admissible side, an
urgent entry uses the boundary, and a permissive jump uses the active side.

Jump candidates return their enclosure and elapsed-time map together. The
orchestrator must not infer that map from the number of parameters. Initial
activation candidates have zero elapsed time; transverse candidates carry
the composed crossing-time map; tube candidates carry the added time axis.
This also corrects the structural domain mismatch when an initial or decreasing
candidate is subsequently restricted by another blocking event: its time map
must use the source domain rather than a nonexistent extra temporal coordinate.

| Crossing information | Flow restriction retained in this refactor | Urgent entry retained in this refactor |
| --- | --- | --- |
| First-entry time certificate (increasing or concave geometry) | elapsed time no later than crossing time | evaluate flow at crossing time |
| INCREASING or CONVEX without a time certificate | guard on admissible side | tube restricted to guard boundary |
| Concave guard-maximum certificate | split using critical time and guard maximum | boundary and entry direction |
| CONCAVE without a certificate | upper: existing two temporal samples; lower: discard continuation | boundary and entry direction |
| UNRESOLVED | upper: existing two temporal samples; lower: discard continuation | upper: boundary and direction; lower: discard candidate |
| NEGATIVE, POSITIVE or DECREASING | initial phase has already handled activation | restrict initial source to active side |
| Invalid geometry or incompatible certificate | rejected at construction | rejected at construction |

Permissive activation does not require a crossing classification in the current
algorithm. Impact applies its directional condition after the common entry
handler. The table describes existing numerical behavior; shared code between
UPPER and LOWER is not itself a proof of inner reachability.

The maximum-based flow handler keeps the named hitting, missing and past subsets local.
It returns additional branches to a single list traversal. Primary branches
retain their original positions; additional branches follow in the same order
as before. No handler holds a reference into a vector-backed list while that
list is extended.

## Orchestration order

1. Validate capabilities and process invariants and initially active events.
2. Compute the numerical step, crossing information and timing.
3. Construct reach and endpoint using the independent endpoint certificate
   when the complete duration is taken.
4. Apply all blocking guards to reach, and non-deferred guards to endpoint.
5. Copy the endpoint branches used for future continuation.
6. Apply guards with a first-entry certificate, and increasing guards, after that copy.
7. Record horizon-limited reach, final and intermediate contributions.
8. For each activating event, construct its jump candidate and elapsed map,
   enforce the horizon and apply every other blocking event.
9. Apply the reset with target labels and auxiliary outputs; return the local
   `HybridEvolutionStepResult` to the backend.

The distinction between a reported endpoint and a continuation checkpoint is
intentional. Guard processing order, simultaneous events and the closed time
horizon must survive the extraction. Handlers must preserve validated
`definitely`/`possibly` decisions and original numerical tolerances.

## Crossing information and future capabilities

`CrossingData` stores `CrossingGeometry` independently of a variant containing
no evidence, `CrossingTimeCertificate` or `GuardMaximumCertificate`. The old
`CrossingKind` and its TRANSVERSE/GRAZING/DEGENERATE tags are removed without
aliases. UNRESOLVED means that sign and curvature were not established; it
makes no claim of actual geometric degeneracy.

A time certificate owns the first-entry map. A maximum certificate owns the
critical-time map, the guard range and state bounds at the maximum. Constructors
check finite, nonempty domains/ranges and consistent chart dimensions; the
numerical producer remains responsible for certifying their mathematical
meaning. Certificate data are immutable through the API. A maximum certificate
forces concave geometry; a time certificate requires increasing, convex or
concave geometry. Geometry-only data never claim a numerical certificate.

The private `analyze_crossings` function receives the configured solver and
an autonomous state-only flow. `_compute_crossings` remains the evolver strategy
hook and delegates to that analyzer. The extraction preserves all predicates,
solver settings, the 15/16 reduced search interval and fallback choices:

| Numerical case | Geometry retained | Evidence/result |
| --- | --- | --- |
| Strictly increasing, root solve succeeds | INCREASING | First-entry time |
| Strictly increasing, EquationSolverException | INCREASING | None |
| Strictly increasing, unexpected runtime error | — | Propagate original error |
| Strictly decreasing | DECREASING | None |
| Strictly convex | CONVEX | None |
| Concave, certified maximum strictly negative | — | Omit event |
| Concave, positive maximum and first root found | CONCAVE | First-entry time |
| Concave, positive maximum but first root solve fails | CONCAVE | Guard maximum |
| Concave, maximum guard range contains zero | CONCAVE | Guard maximum |
| Concave, critical solve fails (equation or runtime error) | CONCAVE | None |
| Sign and curvature unresolved | UNRESOLVED | None |

The nested first-entry runtime-error path retains its existing outer fallback
to concave geometry without evidence. A guard-maximum certificate therefore
does not mean tangency: its maximum may be strictly positive after failed root
solving. Critical times in the past remain valid information. Timing, creep
and deferred endpoint guards test certificate presence; a concave first-entry
certificate receives the same treatment as the former TRANSVERSE result.
Backends retain geometry and evidence kind as `CrossingSummary`, without
retaining numerical patches in transition provenance.

ODE classification uses a single-valued vector field and Lie derivatives.
For ODAI with measurable inputs, second derivatives and a smooth crossing-time
map cannot be assumed. An outer structured tube does not certify a family of
realizable trajectories for LOWER. New handlers must declare the numerical
evidence and semantics they require, and reject unavailable capabilities.

## Work plan and acceptance criteria

- [x] Extract the private event policy and independent flow/jump handlers.
- [x] Return elapsed-time maps explicitly and migrate all callers without
      compatibility aliases for removed hooks.
- [x] Centralize the before/after-continuation restriction phase and preserve
      branch order, competing events, reset history and horizon handling.
- [x] Validate existing mathematical suites and add focused cases for distinct
      jump/flow behavior, upper/lower unresolved geometry and grazing splits.
- [x] Separate geometric properties from typed crossing evidence, including
      explicit unavailable/unsupported states.
- [x] Extract crossing analysis from the evolver and audit its solver fallbacks.
- [x] Decompose timing, creep and unwinding while preserving their strategies.
- [x] Replace parallel endpoint/continuation lists with explicit paired branch
      results and make result publication a small orchestration phase.
- [x] Document and verify the supported capability matrix for each handler;
      extend ODAI events only with corresponding mathematical certificates.

The first six tasks form the implemented behavior-preserving extraction. Changes
to inner guarantees, impact assumptions, approximation rules or unresolved-case
handling are separate numerical work. Run the existing hybrid evolution,
contract, backend and tree regressions, then the full Evolution suite. Use
analytic witnesses and exclusion bounds for new mathematical checks; do not
replace them with tests which merely repeat policy-table entries.

Format is currently non-blocking by decision; maintain the functional
compiler/test/coverage/sanitizer checks before expanding numerical capabilities. Keep simulator work
after the required hybrid capabilities. Full Python/C++ tutorial parity and
the bindings it needs remain the final Evolution task.

## Verified extraction

The first four tasks above are implemented in `src/hybrid_transition_handlers.hpp`
and `.cpp`. The protected `_apply_guard` and `_apply_guard_step` hooks are removed;
all repository callers use the private operations directly. There are no naming
aliases or compatibility wrappers for those removed hooks. Crossing analysis
and typed evidence are implemented in the dedicated crossing components. Timing decomposition, paired continuation results and the capability
audit have since been implemented. ODAI UPPER events use conservative outer
guard intersections; first-entry/trajectory capabilities and LOWER integration
remain unavailable without further mathematical certificates.

The Debug build treats warnings as errors. Validation passes all 30 Evolution
CTest targets, including Python, and the 13 System C++ targets in the consuming
build. The hybrid contract suite contains 13 cases. Four new analytic cases
check transverse versus permissive elapsed maps in both semantics, impact
direction and zero-time source candidates with competing guards, unresolved
upper/lower restrictions, and stable grazing branch order with explicit feasible
trajectory witnesses. Existing hybrid, backend, tree and inclusion tests remain
active. The new crossing suite has five analytic cases covering certificate
invariants, increasing roots and decreasing direction, concave maxima and
strictly negative exclusion, convex/unresolved geometry, equation/runtime solver
fallbacks, and first-entry provenance for concave geometry. Analytic residuals
check root and critical-time maps; the positive-maximum fallback explicitly
excludes an interpretation as tangency. The public header is included first.

The complete consuming build passes 43/43 CTest targets. A fresh install supports
the independently built C++ tutorial and the existing Python tutorial.
`apply-format` and `check-format` pass for all 83 tracked C++ files before commit.
Full Python/C++ tutorial parity remains deferred to the final task.

## Timing, creep and unwinding extraction

`src/hybrid_timing.hpp/.cpp` owns horizon selection, guard creep and clock
unwinding. The existing evolver `_estimate_timing` hook delegates with the
configured function factory and solver, preserving the step strategy.

Horizon selection retains the exact final-time combinatorial choice and
halved time-scaling rule. The guard-creep eligibility test excludes permissive
events and checks whether the full-step endpoint is empty. Solver-based creep
preserves guard ordering, the 0.9375 default, clobbering, existing exception
fallback and removal only of guards with successfully computed implicit creep
times. Clock unwinding runs only after composing the time maps, and only with
no remaining crossing, nonzero clock spread, a constant step and
`BEFORE_FINAL_TIME`; both previous clock-narrowing formulas are unchanged.

No result publication, backend selection, reset or workload scheduling moved
into timing. The existing hybrid timing, creep, solver-budget, crossing,
backend and tree suites are the regression criteria. This extraction does
not change numerical guarantees or add ODAI event capabilities. Paired
continuation branches and result publication remain the next task.

## Paired endpoint branches and continuous-result publication

`src/hybrid_step_results.hpp/.cpp` introduces `HybridEndpointBranch`,
pairing the reported endpoint with its continuation checkpoint. The pair is
created after non-deferred blocking restrictions and before deferred
first-entry/increasing restrictions. Deferred restrictions affect the reported
member only. If a geometric handler splits a reported branch, the additional
branches inherit the original, unmodified checkpoint; original primary
branches precede the additional branches, matching guard handler order.
No parallel list-position assumption remains.

`publish_step_contributions` is the sole continuous-result publication
phase of the hybrid equation step. It retains the prior exact-final-time
case (final contributions before reach), the non-final case (intermediate
and continuation per branch before reach), the unbounded continuation
checkpoint when a reported endpoint must be horizon-bounded, and the
closed reach-horizon final slices. Jump construction and reset still follow
continuous publication and remain owned by the evolver. The private
component publishes only into a local `HybridEvolutionStepResult`; backends
remain responsible for orbit publication and scheduling. Numerical
crossing/solver policy, LOWER limitations and ODAI capability are unchanged.

## Validated handler capability matrix

The live selectors `flow_guard_capability` and `jump_guard_capability`
in `src/hybrid_transition_handlers.cpp` are the source of the dispatch
matrix, rather than a second, independently maintained table. They return
the chosen private method, required typed crossing evidence, whether a
timed flow is used, and whether a guard Lie derivative is required.
`select_flow_handler` and `select_entry_handler` dispatch through those
results. All complete event-analysis calls still require the existing
single autonomous state-only flow patch, independent endpoint, and
supported chart; this is enforced by `_evolution_step`. This matrix
documents the *present ODAE implementation*, not generic ODAI support.

| Crossing certificate or geometry | Flow restriction (UPPER / LOWER) | Blocking jump (UPPER / LOWER) | Required numerical information |
| --- | --- | --- | --- |
| Certified first-entry (INCREASING, CONVEX or CONCAVE) | FIRST_ENTRY_BOUND / same | FIRST_ENTRY / same | Certified entry-time patch on original state chart, structured flow |
| Certified concave maximum | MAXIMUM_SPLIT / same | DIRECTED_BOUNDARY / same | Critical-time patch, maximum guard range/state bounds; Lie derivative for directed entry |
| INCREASING or CONVEX without entry certificate | GUARD_SIDE / same | BOUNDARY / same | State guard, validated local tube; no root-time claim |
| CONCAVE without certificate | TWO_SAMPLES / DISCARD | DIRECTED_BOUNDARY / same | Upper flow samples at two fixed times; blocking jump requires guard Lie derivative |
| UNRESOLVED | TWO_SAMPLES / DISCARD | DIRECTED_BOUNDARY / DISCARD | Upper: two flow samples or directed boundary; lower conservatively discards |
| NEGATIVE, POSITIVE, DECREASING | KEEP / same | SOURCE_ACTIVATION / same | Initial-phase treatment and source guard; no invented future entry time |

Event-policy effects are independent of that crossing matrix:

| Event | Flow handler | Jump handler | Additional requirement |
| --- | --- | --- | --- |
| INVARIANT | Blocking flow handler | Not applicable | No reset |
| PROGRESS | Blocking flow handler | Not applicable | No reset |
| PERMISSIVE | No blocking flow restriction | TUBE_ACTIVATION in both semantics, independently of crossing geometry | Tube with elapsed-time coordinate and active-side guard |
| URGENT | Blocking flow handler | Blocking entry from matrix | Event reset, competing blockers, closed time horizon |
| IMPACT | Blocking flow handler | Blocking entry from matrix | Additional impact-direction Lie derivative even with a first-entry certificate |

The `uses_timed_flow` and `requires_lie_derivative` flags describe
*operations performed*, not a proof that an inner approximation exists.
In particular, `LOWER` is an existing ODAE algorithmic path, not a
general guarantee for inclusion solutions. A positive concave maximum
does not imply tangency, and `UNRESOLVED` does not imply degeneracy.

**Dynamical capability gate (current):**

| Mode | UPPER | LOWER | Events / guards / invariants |
| --- | --- | --- | --- |
| Autonomous state-only ODAE, no fixed parameters | Existing validated local equation method | Existing equation method (no new inner claim) | Existing event policy and crossing classifier, subject to certified local flow layout |
| Autonomous ODAI, no fixed parameters | Explicit structured or box method; necessary blocking restrictions and conservative urgent/permissive/impact jump candidates | Rejected, including at time zero | State invariants and PROGRESS/URGENT/IMPACT blockers; UPPER urgent/permissive active-side and impact-boundary candidates |
| Explicit-time or fixed-parameter ODAE/ODAI modes | Rejected by current hybrid mode validation | Rejected | Not reached |

ODAI with arbitrary measurable input functions lacks, at present, a
certified trajectory-level guard derivative and first-entry/maximum
analysis compatible with this equation classifier. Neither an outer box
nor an independent endpoint certificate supplies these properties.
Conservative UPPER ODAI transitions and blockers can use outer reach/guard
intersection without such certificates. IMPACT uses the full closed boundary
without directional exclusion, as detailed below. Certified first entry and
LOWER transitions require an explicit compatible numerical capability and a
separate proof of the semantics it supports.

The focused handler contracts exercise geometric witnesses, root-time
maps, ambiguous sampled/discard behavior, impact direction, branch order,
and incompatible certificate rejection. The inclusion tests exercise state-only
event restrictions and conservative UPPER jumps with both explicit integration choices; unsupported LOWER, time
and fixed-parameter capabilities are rejected before integration.

## Initial ODAI capability: state invariants under UPPER

The private `src/hybrid_inclusion_invariants.hpp/.cpp` component initially
admitted only `EventKind::INVARIANT`. Its UPPER capability is extended below
to conservative urgent and permissive reset candidates. The existing explicit
structured or box integration method remains required, and LOWER remains
unsupported. No input trajectory or first-contact time is certified.

An invariant `g(x) <= 0` is a *necessary pointwise condition* on
every valid solution. The validated outer inclusion reach tube and the
independent endpoint are each intersected with that state constraint.
At the initial, stopping and post-reset phases the same state-only
invariant restriction is reused, without querying a nonexistent
single-valued vector field. The original state, time, dwell, auxiliary
mapping and event-history representations remain unchanged.

**Soundness boundary:** intersecting a validated outer tube with
`g(x) <= 0` retains an outer approximation of solutions that stay
in the invariant. It does **not** prove that each point in the
intersection is dynamically realizable, that a whole trajectory
respects the invariant between snapshots, or that the first hitting
time is known. This necessary-condition filter does not, on its own, justify source
coverage, first-entry constraints, or any LOWER transition. Under UPPER,
however, a guard intersection of the outer reach can generate additional
possible (including spurious) transitions without such certificates.

Analytic regression cases exercise an autonomous inclusion `x' = u`,
`u in [-1,1]`, under `x <= 0`, across structured and box
integration, including both reach and endpoint state constraints.
A closed zero-horizon source slice and rejection of ODAI LOWER semantics
are checked independently.

The next *certified* event extension requires an explicit trajectory
or crossing certificate from the inclusion solver, with separate
justification for each event kind and semantics. UPPER conservative jumps
can instead use the validated reach/guard intersection, as detailed below.

## Conservative ODAI URGENT/PERMISSIVE transitions in UPPER

The key distinction is approximation semantics, not whether every
candidate point actually takes a transition. In UPPER it is sufficient
to retain all feasible transition states/times in the approximation;
additional spurious jumps are permitted. LOWER requires a suitable
all-points realizability or coverage guarantee on a selected subenclosure
before claiming that subenclosure can take the transition.

For autonomous inclusion modes, the already validated outer reach tube
contains all realizable state/time pairs in the step. A candidate is
constructed by intersecting this tube with the guard `g(x) >= 0` and
the closed time horizon, then applying the existing target reset. The
candidate retains its time chart and input-independent previous event
history; the resulting successor retains the event kind and history.
This construction is conservative: **the time is not a certified first
entry**, and no source-wide coverage is inferred. It can retain
post-first-contact candidates, including spurious jumps. The ordinary
continuous branch is retained conservatively to avoid silently dropping
realizable inclusion trajectories for urgent events.

The inclusion transition extraction does not call
`dynamic_function` or `lie_derivative`: those operations would
silently treat a set-valued field as a single trajectory field.
The transition handler operates only on validated outer reach
enclosures. The same initial-event and invariant filtering phases apply.
`EventKind::PROGRESS` is supported as a blocking state constraint
on the validated outer reach and endpoint, without reset and without a
first-hit certificate. The initial/zero-horizon source is not
retroactively constrained by PROGRESS. IMPACT support is described separately
below: it uses only the necessary boundary condition and does not claim a
directional certificate. ODAI LOWER remains explicitly rejected.

The analytic regression covers `x'=u`, `u in [-1,1]`,
guard `x>=1/32`, and reset `x'=0`, with both URGENT and
PERMISSIVE kinds and both existing structured and box integrators.
The feasible witness `u=1` activates the guard at `t=1/32`.
The test checks that target final states and source reach are
retained without claiming all candidate transitions are realizable.
At zero horizon, the stopping-time phase records the invariant-valid
terminal source slice before initial event processing and does not schedule
a transition. The zero-horizon regression checks termination, not an
urgent reset. LOWER/PROGRESS/IMPACT are not expanded by this change.

Do not reuse the ODAE first-entry, concavity or directional impact machinery
for ODAI without a compatible numerical certificate from the integrator.

## Single-source initial guard activation and deferred collective LOWER verification

The private `classify_initial_guard` policy distinguishes `EXCLUDED`,
`UNRESOLVED`, and `FULL_SOURCE` against the entire *invariant-valid*
source. `FULL_SOURCE` is established by universal guard satisfaction
on that source, never by satisfying the guard only after intersecting or
selecting a favorable subset. `EXCLUDED` requires a certified empty
guard intersection. `UNRESOLVED` includes both genuinely partial and
numerically inconclusive activation.

The initial-event phase carries `Semantics` explicitly to the common
policy: UPPER may create a possible, potentially spurious reset for
`UNRESOLVED`, but LOWER creates that initial reset only for
`FULL_SOURCE`. The event policy still determines whether the guard
blocks flow, whether a reset is optional, and whether the event is
inspected on entry. Existing complete-source urgent provenance is
preserved; permissive transitions remain optional and therefore do
not assert urgent full-source resolution.

No LOWER subdivision is enabled. Numerical splitting is conceptually
compatible with LOWER *only when every part of the complete partition
is retained*. One cannot certify a LOWER transition by selecting
the guard-active pieces and discarding the others. A future tree-level
collective guard decision could certify full activation when **every**
part of a complete partition is certified active; that aggregate check
does not exist yet, and is not simulated by independently certifying
a few child enclosures.

The new backend contract tests use constant ODAE flow, full guard
activation on `[1/4,1/2]`, and partial activation on
`[-1/4,1/4]`, for both URGENT and PERMISSIVE transitions.
LOWER must preserve a fully active reset, reject a partial reset
without discarding the whole computational state, and UPPER must
retain a conservative possible reset with unresolved provenance.

This change addresses *initial* activation only. Time-interior ODAE
LOWER handlers and the guarded-reach UPPER ODAI path retain their
separate numerical policies. The inclusion bridge still rejects
LOWER modes: the present UPPER integration methods do not furnish a
LOWER reachable-set representation merely because a guard is fully
active. No claim of such support is made here.

## Single-source activation decision contracts

A focused analytic contract distinguishes a guard that is impossible on
the source, a guard active on only a subset, and a guard satisfied by
every point of the source. On the closed parameter interval
`[-1/4,1/4]`, the constraints `x>=1`, `x>=0` and `x>=-1`
exercise `EXCLUDED`, `UNRESOLVED` and `FULL_SOURCE`, respectively.
UPPER is permitted for unresolved candidates; LOWER only for whole-source
activation. Both reject excluded candidates. Invalid semantics must not
silently select either branch. These tests do not claim inner ODAI
reachability or a collective full-partition check.

The ODAE first-entry certificates used during a continuous step must
not be confused with the simple initial-source predicate: a certified
parameter-dependent entry time may establish valid per-trajectory
crossings even where a fixed-time reach slice activates only a subset.
Changing such numerical handlers requires separate proof and tests.

## UPPER ODAI PROGRESS blockers (no transition)

A PROGRESS guard `g(x)>=0` is an upper-blocking condition: a valid
continuous branch must stay on `g(x)<=0`. `apply_inclusion_flow_blockers`
adds this necessary pointwise condition to every certified inclusion
reach tube and independently calculated endpoint, alongside invariants.
No reset or target mode is queried for PROGRESS. As with UPPER invariant
filtering, the guarded outer tube may contain spurious trajectories:
the intersection is a necessary set condition, not a first-hitting-time
certificate. An urgent or permissive jump candidate formed from the
restricted reach must also respect PROGRESS blockers.

Unlike invariants, a PROGRESS guard is not imposed on the initially
terminated source at a closed zero-time horizon; the initial phase
retains its existing ODAE-style blocking policy. This does not enable
LOWER inclusion modes or directional IMPACT. Focused tests use
`x'=u, u in [-1,1]` and `x>=1/32` as PROGRESS, checking both
structured and box integrators, constrained source reach, lack of
spurious reset provenance, and preservation of the t=0 terminal set.

## Concurrent UPPER ODAI blocking guards

A validated inclusion outer reach and its independent endpoint are now
filtered by every blocking guard, following the same event-kind policy
used by ODAE: INVARIANT, PROGRESS and URGENT restrict the continuous
branch to the non-active side; PERMISSIVE does not. This does not
identify a first entry, nor justify LOWER reachability. Closed
boundaries are retained for an urgent candidate, so simultaneous
URGENT events with a common threshold are not ordered by declaration
or silently dropped. Jump candidates are constructed separately from
the guarded outer reach before resets, so competing blockers apply
without requiring a single-valued Lie derivative.

The regression uses `x'=u`, `u in [-1,1]`, two URGENT
transitions with `x>=1/32` to distinct zero-dynamics targets, and
checks that both target branches survive for the structured and box
methods. A separate PERMISSIVE case checks that optional jumps do
not prevent the source branch from reaching the horizon. Neither
test claims the two simultaneous jumps can both occur along one
concrete execution. The implementation retains their conservative
UPPER alternatives and does not add a collective LOWER guard check.

## Conservative UPPER ODAI IMPACT boundary candidates

An admissible impact under the existing blocking event policy lies on the
closed guard boundary `g(x)=0`. Intersecting a validated outer reach with that
boundary and the closed time horizon contains every feasible impact state/time.
As with the urgent/permissive path, it may also contain spurious candidates.
This necessary-condition construction is sufficient for UPPER possibility;
it does not certify an actual impact, crossing direction, first entry, smooth
event-time map or whole-source coverage.

`inclusion_upper_transitions` now extracts IMPACT resets using the common event
policy. The reach and independent endpoint respect all blocking guards,
including IMPACT, while each impact candidate explicitly adds `new_guard`
(`g=0`) before reset. This explicit equality also preserves the helper's
boundary contract independently of prior flow filtering. No single-valued
`dynamic_function` or Lie derivative is queried for inclusion modes. The
ordinary ODAE directional handlers retain their existing numerical rules.

No directional exclusion is applied to these UPPER candidates: measurable
inputs need not supply a unique velocity or derivative at contact. In
particular, a tangency with zero guard derivative is retained. Inward contacts,
repeated boundary contacts and post-first-contact candidates can remain
spurious; their presence does not assert realizability. The ordinary finite
event budget still limits reset sequences. LOWER remains rejected, even at
zero horizon. The terminal zero-horizon phase continues to retain the
invariant-valid source without scheduling a reset.

Competing invariant, progress, urgent and impact conditions remain on the
pre-reset enclosure; permissive events do not block. Closed simultaneous
boundaries preserve separate impact/urgent alternatives without imposing an
event priority. Candidates retain the original state/time chart until reset,
then the existing reset operation changes state/auxiliary labels and clears
dwell time. Transition provenance records IMPACT with no crossing certificate
and unresolved full-source resolution.

Five analytic regression cases check:

- The exact witness `x(t)=t-2`, `y(t)=0` touches `g=-x^2` at absolute `t=2`
  with zero derivative. Every post-reset constraint admits the witness;
  off-boundary points are excluded. Reset dimension, auxiliary output, clock,
  dwell and event provenance are checked.
- A horizon before that contact excludes it; the closed horizon at contact
  retains it.
- A preceding PROGRESS blocker excludes a later contact, while simultaneous
  IMPACT and URGENT boundaries retain both reset alternatives.
- For `x'=1`, `y'=u`, `u in [-1/8,1/8]`, the genuine `u=0` tangency survives
  both structured and box integration on the flat and tree backends. Zero
  horizon and unsupported LOWER behavior remain explicit.
- An ODE reset enters that inclusion at a parameter-dependent arrival time;
  the subsequent contact retains event order, feasible reset states and the
  full analytic contact-time interval, including the final horizon.

Arrival and local time can be correlated: a raw time-range box need not itself
lie below the horizon even when the constrained enclosure does. The regression
proves the bound from the retained constraints, using `c<=u` and
`c-time+T>=u` to establish `time<=T`, rather than discarding that correlation or
relaxing the horizon. Exact analytic contact witnesses separately check the
closed boundary and exclusion of future/off-boundary points.


Validation of this extension uses the macOS Debug build with warnings as errors:
44/44 integrated CTest targets (31 Evolution including Python and 13 System C++),
the installed standalone C++ tutorial and the existing Python tutorial.
`apply-format` and `check-format` pass on all 90 tracked C++ files, also correcting
formatting accumulated in the preceding hybrid refactoring commits. The remote
functional CI of the new commit must be checked independently of these local
results and the five successful functional workflows of predecessor `85c0781`.
The four-family roadmap and full Python/C++ tutorial parity remain open.
