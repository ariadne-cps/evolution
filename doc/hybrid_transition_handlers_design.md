# Hybrid transition handlers

The hybrid evolver separates numerical capabilities, event policy, crossing
analysis, geometric guard operations and orchestration. Refactoring these
responsibilities precedes further ODAI event support. Event policies, guard handlers and crossing analysis are separate private
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
- [ ] Replace parallel endpoint/continuation lists with explicit paired branch
      results and make result publication a small orchestration phase.
- [ ] Document and verify the supported capability matrix for each handler;
      extend ODAI events only with corresponding mathematical certificates.

The first six tasks form the implemented behavior-preserving extraction. Changes
to inner guarantees, impact assumptions, approximation rules or unresolved-case
handling are separate numerical work. Run the existing hybrid evolution,
contract, backend and tree regressions, then the full Evolution suite. Use
analytic witnesses and exclusion bounds for new mathematical checks; do not
replace them with tests which merely repeat policy-table entries.

Run `apply-format` and `check-format` before committing. Keep simulator work
after the required hybrid capabilities. Full Python/C++ tutorial parity and
the bindings it needs remain the final Evolution task.

## Verified extraction

The first four tasks above are implemented in `src/hybrid_transition_handlers.hpp`
and `.cpp`. The protected `_apply_guard` and `_apply_guard_step` hooks are removed;
all repository callers use the private operations directly. There are no naming
aliases or compatibility wrappers for those removed hooks. Crossing analysis
and typed evidence are implemented in the dedicated crossing components. Timing decomposition, paired continuation results and the capability
audit remain the unchecked tasks above, before new ODAI event capabilities.

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
