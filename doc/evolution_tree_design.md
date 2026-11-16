# Evolution tree design

The hybrid evolver separates numerical propagation from execution and result
storage. Flat and Tree implement the same step/result contract. Tree has one
checkpoint scheduler with configurable admission, grouping, evaluation timing,
dispatch, numerical recombination and acceptance. All modes use the same
DynamicWorkload for sequential and concurrent execution.

The generic append/prune genealogy belongs to
`ariadne-cps/tree/doc/derivation_tree_design.md`. Hybrid trajectory semantics and
all scheduling decisions belong to Evolution.

## Numerical work and Orbit contributions

`HybridEvolutionWorkItem` contains a HybridEnclosure and a processing phase:
PROCESS_STARTING_EVENTS or EVOLVE_CONTINUOUSLY. The enclosure retains its own
numerical configuration. The evolver's callback returns a complete
`HybridEvolutionStepResult` containing:

- explicitly classified reach, intermediate and final contributions;
- successors classified as continuation, transition, subdivision or reconditioning;
- transition event, kind, crossing information and whole-source resolution evidence;
- a local stop reason.

An intermediate contribution does not schedule work implicitly. A final
contribution can coexist with successors. Exhausting work does not invent a
final contribution. The backend publishes a successful result before any of its
successors executes. A failed callback leaves incomplete history and prohibits
Orbit release; outstanding workers are joined before propagating failure.

`HybridEvolverBase::set_evolution_backend_factory()` clones immutable factory
settings. Each orbit() call creates independent run state, supplies semantics,
runs the same numerical callback and releases a conventional Orbit. Numerical
propagation and termination rules are independent of the selected backend.

## Flat backend

`FlatHybridEvolutionBackendFactory{}` is the evolver default. It has no
backend-specific configuration. All successors enter DynamicWorkload
immediately; contributions are accumulated directly into Orbit under a result
lock. It records no genealogy and performs no recombination.

Both backends use ThreadManager's concurrency setting: zero executes directly,
positive concurrency uses available workers. Concurrency is not duplicated in
backend configuration.

## Trajectory history

`HybridEvolutionTree` is
`DerivationTree<HybridEvolutionTrajectory, HybridEvolutionTrajectoryMetadata>`.
Each node owns an ordered vector of entries. An entry stores a work item, its
optional same-trajectory predecessor and its optional committed result.
`HybridEvolutionCheckpointId` identifies a node and entry index.

Pure continuous continuations and reconditioning append entries in the same
node. A sole urgent/impact transition continues that node only with certified
FULL_SOURCE resolution and whole-trajectory eligibility. Permissive/unresolved
transitions and distinct alternatives create child trajectories. Every
subdivision creates explicit child trajectories, under every admission policy;
scheduling no longer changes how splitting is represented.

Node metadata records the origin entry/successor, local completion, pending
entry count, whole-trajectory eligibility, recombined checkpoint provenance and
an optional recombination time barrier. A complete node can have incomplete
descendants. Completion is published only after its results and children.

There is one retirement representation: the trajectory's retired_entries.
Processed entries retain their original results; replaced pending entries retain
provenance and have no fabricated result. There is one replacement-provenance
representation: recombined_checkpoints, including for sibling-only merges.
Replacements are attached below their nearest common ancestor. Their position
in the genealogy does not override admission dependencies.

## One checkpoint scheduler

Every policy shares:

- ready: eligible checkpoints, each backed by a workload token;
- waiting: unprocessed checkpoints withheld from execution;
- in-flight count: claimed numerical computations that have not published.

There is no second parent-first execution branch or separate per-parent
ready/waiting data structure. The scheduler claims and snapshots work under its
lock, performs numerical evolution without the lock, then publishes results,
creates successor records and applies admission under the same lock.
Recombination currently runs under that publication lock as well. Numerical
propagation can run in parallel; concurrent merge construction is not implemented.

A batch is complete only when ready is empty and the in-flight count is zero.
Retired workload tokens can outlive their checkpoints and become harmless
no-ops. A policy cannot select unknown, retired or duplicate pending work; it
must admit a nonempty selection if the batch is complete and waiting remains.

### Admission

`HybridCheckpointAdmissionPolicy` separates admission from genealogy and dispatch.
It decides whether a successor is immediately ready, selects eligible waiting
IDs after publication and chooses the queue used for recombination. select_next
receives batch_complete and can return no work at intermediate publications.

| Policy | Immediately ready | Waiting release | Merge queue |
| --- | --- | --- | --- |
| ParentCheckpointAdmissionPolicy | Continuations, reconditioning, subdivisions and certified same-node transitions | When the producing continuous family has finished | WAITING |
| ImmediateCheckpointAdmissionPolicy | Every successor | No ordinary waiting | READY |
| LocationCheckpointAdmissionPolicy | Continuous successors in the current location | At the batch barrier, all waiting IDs in the first waiting location | WAITING |

Parent admission is the default. A continuous family includes a trajectory and
all descendants born by subdivision, excluding its discrete descendants.
Discrete successors wait for all pieces of that family. Independent completed
families can release descendants while other families still run, so parent
admission does not impose a global level barrier. For a recombined replacement,
all source-family dependencies are followed through checkpoint provenance;
attachment below a completed common ancestor cannot release it prematurely.

Location admission runs one location at a time. Every transition waits,
including self-loops and certified transitions preserving trajectory identity.
At a barrier all compatible waiting checkpoints can be recombined across
parents, then all checkpoints of the first waiting location become ready.
Waiting for other locations can persist through multiple batches and absorb
later compatible arrivals. There is no physical-time synchronisation or clock
shift. Arrival order decides the next location, not a configurable priority.

### Grouping

`HybridRecombinationGroupingPolicy` is the single grouping interface. It receives
explicit pending checkpoint IDs rather than deriving eligibility from topology.

- NoRecombinationPolicy disables numerical recombination.
- AllCompatibleRecombinationPolicy groups all compatible pending checkpoints,
  without time-spread or cardinality limits. Disjoint physical time intervals
  are allowed; the numerical recombiner covers their clock values.
- WindowRecombinationPolicy(max_group_size, max_time_spread, scope) bounds group
  cardinality and the hull of time ranges. ANY allows nonadjacent checkpoints
  but requires overlap with the group anchor. SIBLINGS selects consecutive
  compatible sibling seeds and permits disjoint times within the time span.
- ProximityRecombinationPolicy(maximum_center_distance, max_group_size,
  maximum_enclosure_radius, scalings) groups spatially nearby checkpoints across
  arbitrary parents and disjoint clock intervals. Its normalized infinity-norm
  diameter bounds every pair of centres, rather than linking a chain of close
  neighbours. State bounds are evaluated over the current reduced parameter
  domain; the optional preview radius bounds their unscaled bounding hull.
  The numerical context must also specify a radius limit when the actual
  recombined output must satisfy it: a small source hull does not guarantee a
  small validated output.

AllCompatible and Window share HybridRecombinationScope: ANY or SIBLINGS.
SIBLINGS requires the same immediate genealogy parent and successor meaning
(kind, event and event kind). Replacement provenance preserves that meaning
when an early merge is considered again. ANY imposes no parent restriction.

Proximity uses an immutable balanced spatial index built for each selection,
evaluates each candidate's geometry once, and greedily selects disjoint groups
of at most the requested size. It does not promise an optimal partition.
Unprocessed replacement seeds are excluded until they execute: repeatedly
merging a held replacement can accumulate residual parameters without any
reconditioning or physical evolution. Original singleton candidates remain
available for later arrivals. The numerical recombiner still preserves all
source constraints and clocks; grouping only controls which union is attempted.

The backend independently requires equal locations, processing phases, state
and auxiliary spaces, exact previous-event histories and structurally identical
symbolic auxiliary maps. Neither common ancestry nor geometric overlap proves
compatibility or the numerical covering contract.

### Evaluation timing

`HybridRecombinationTriggerPolicy` receives a publication/admission boundary,
the pending queue and the global committed-step count. The count has the same
meaning under every admission policy.

- BeforeAdmissionRecombinationTriggerPolicy is the default: before releasing
  eligible waiting work, or after each publication for immediate ready admission.
  It replaces separate parent-completion and location-synchronisation triggers.
- AfterEachStepRecombinationTriggerPolicy inspects the selected queue after
  every committed numerical result, allowing early waiting merges.
- EveryNStepsRecombinationTriggerPolicy(n) inspects every n global steps and
  flushes waiting work before admission. Immediate admission does not turn the
  flush into an inspection on every step.

At most one evaluation is made per publication. After recombination, eligibility
is recalculated before admitting replacements. Early waiting replacements stay
waiting until their admission conditions hold.

### Dispatch and acceptance

`HybridCheckpointDispatchPolicy` orders ready work independently of admission.
ImmediateCheckpointDispatchPolicy preserves FIFO. BoundedTimeLagCheckpointDispatchPolicy
moves a candidate to the tail while an earlier-time checkpoint of a different
trajectory is ready within max_time_lag, at most max_deferrals times. This
preserves the checkpoint and its physical time and does not wait indefinitely
for partners. Under location admission it only affects order inside that location;
it does not change when waiting sets become merge candidates.

`HybridRecombinationAcceptancePolicy` assesses a numerically valid result.
AcceptAllRecombinationsPolicy accepts every sound effective reduction.
RadiusTrendRecombinationAcceptancePolicy additionally requires a previously
executed continuous continuation with positive time progress, bounded observed
radius growth and bounded output amplification. It is an empirical quality gate,
not a contraction proof. It supports READY only: choosing it with enabled
WAITING recombination is rejected by both backend and factory. With no-op
grouping the unused acceptance gate has no effect.

## Temporal progress toward the run horizon

Flat and Tree share HybridEvolutionProgressTracker and expose a consistent
HybridEvolutionProgress snapshot through backend.progress(). Let t0 be the
initial enclosure's lower absolute-time bound, T the requested final time, and
D=max(0,T-t0). For the currently represented pending branches:

```text
C = sum(clamp(lower_time(checkpoint) - t0, 0, D)) + completed_terminal_branches * D
R = sum(max(0, T - lower_time(checkpoint)))
progress = C / (C + R)
```

C is credited_time, not accumulated numerical computation cost. Every new branch
inherits its elapsed prefix: splitting or a permissive transition at time 120
with t0=0 and T=200 creates another branch already 60% through the horizon,
not a new branch starting at zero. All sums include ready, waiting and in-flight
checkpoints. Claiming work, changing admission state or reordering the queue does
not alter the estimate. The bounds are conservative clock estimates, without
synchronising clocks or averaging a time interval's endpoints.

Publication removes the source's credit and remainder and registers every
successor's corresponding credit and remainder. An early-terminal checkpoint
has no future evolution to perform and contributes one complete horizon of
completion credit. A committed merge removes its source estimates and adds
its replacements' estimates; retired sources are not counted as terminal branches.
Rejected merges and retired no-op workload tokens leave the accounting unchanged.

For example, one set at time 120 with horizon 200 gives C=120, R=80 and 60%.
Two split pieces at the same time give C=240, R=160 and still 60%. Merging them
back into one gives 60% again. Accumulating hundreds of permissive branches
born at or after time 120 cannot by itself reset the indicator to zero.
The percentage can still change when topology changes the relative weights of
branches at different times, or when a replacement widens its clock bounds.
No monotonic smoothing hides these changes.

The numerical callback's continuous_evolution_time remains a separate diagnostic.
processed_time sums that value once per published computation, independently of
its number of Orbit contributions or successors. The built-in evolver reports
the upper range of its chosen parameter-dependent flow duration, clipped to the
remaining horizon. Structural steps contribute zero. This cumulative cost is
not used as the progress numerator: counting shared prefixes once in the numerator
but every new branch's future in the denominator would collapse progress during
progressive permissive transitions. Invalid cost metadata marks processed_time
unavailable (NaN) without resetting valid clock-based temporal progress.

These estimates are not CPU-time estimates, geometric probabilities or certified
completion bounds. Without a finite horizon or a finite initial/pending time
bound, remaining_time is unavailable and no temporal percentage is estimated
until successful exhaustion. Zero-length horizons avoid division by zero.
The displayed percentage is capped at 99 while work remains; complete is set
only after successful workload exhaustion, independent of final Orbit contributions.

Counters are updated incrementally under the publication mutex, without traversing
the full trajectory history. Snapshots are safe while numerical workers execute.
Factory-created Tree runs receive the same horizon as Flat. A directly constructed
Tree backend can call set_time_horizon() before run(); omission leaves the finite-time
estimate unavailable.

## Structural-operation diagnostics

Recombination logs selection, a begin marker before numerical construction, and
accepted, rejected or failed outcomes. Events include the source count, location,
source time span, maximum source parameter count and numerical elapsed seconds;
accepted merges also report output count, radius and clock range. No-compatible-group
and no-time-progress skips are explicit. These events use the normal structural
logging level, visible at verbosity 3 through the standard evolver entry point.

ADMIT reports the next cohort's location and clock, ready/waiting/active checkpoint
counts and its temporal-progress estimate. In particular, Workload's generic
w/p/c counters count workload tokens, not all hybrid checkpoints. Tree waiting
checkpoints have no runnable token until admission. Thus w=0,p=1 can mean that
the last publisher is serially grouping or recombining many waiting sets under
the publication lock; it does not establish that those sets have already merged.
A begin marker without an outcome identifies numerical construction still in
progress (or a subsequent interrupted execution), whereas an accepted marker
confirms the replacement was committed before the next ADMIT event.

## Numerical construction and acceptance contract

Only UPPER semantics permits enabled recombination. A proposal must have at
least one output, strictly fewer outputs than inputs, no discarded sources,
complete source-index provenance and a satisfied radius request. Failure or
quality rejection leaves all original work runnable. The backend relies on the
recombiner's validated upper-covering contract; it does not infer coverage from
outer projection overlap or source-index bookkeeping alone.

EnclosureRecombinationContext supplies an optional maximum output radius
explicitly. It is not automatically copied from maximum_enclosure_radius in the
evolver. ReferenceEnclosureRecombiner constructs at most one output and does not
partition it to fit a requested radius. Without a recombination radius limit,
the evolver can subsequently subdivide a large replacement.

READY replacements inherit a time-progress barrier preventing immediate
split/recombine cycles at the same physical time. WAITING replacements can absorb
later arrivals without advancing time, since they cannot execute or resplit
before admission. All merged clock ranges must remain finite.

## Configuration examples

One constructor is shared by the concrete backend and its factory. Its arguments
are admission, grouping, trigger, recombiner, context, dispatch and acceptance;
the concrete backend additionally takes the initial enclosure first. Policies
are cloned. Defaults select parent admission and no recombination.

```cpp
// Parent admission, conventional Orbit, no recombination.
evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory{});

// One active location; merge all compatible waiting checkpoints at its barrier.
evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory(
    LocationCheckpointAdmissionPolicy{}, AllCompatibleRecombinationPolicy{}));

// Immediate admission with periodic time-window grouping of ready checkpoints.
evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory(
    ImmediateCheckpointAdmissionPolicy{}, WindowRecombinationPolicy(8u, 0.125),
    EveryNStepsRecombinationTriggerPolicy(8u)));

// One active location; compact waiting work periodically in spatial groups.
// Position units are scaled separately from other state variables.
const Scalings scales(1.0_x, {{RealVariable("x"), 12.0_x}, {RealVariable("y"), 12.0_x}});
evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory(
    LocationCheckpointAdmissionPolicy{}, ProximityRecombinationPolicy(1.0, 8u, 12.0, scales),
    EveryNStepsRecombinationTriggerPolicy(128u), ReferenceEnclosureRecombiner{},
    EnclosureRecombinationContext(12.0_x)));

// Parent admission, bounded consecutive sibling groups before release.
evolver.set_evolution_backend_factory(TreeHybridEvolutionBackendFactory(
    ParentCheckpointAdmissionPolicy{},
    WindowRecombinationPolicy(8u, 0.125, HybridRecombinationScope::SIBLINGS)));
```

## Verification and deferred work

Tests compare sorted complete descriptions and exact collection counts against
Flat without recombination, for urgent/permissive transitions, splitting,
invariants and termination limits. Enabled recombination tests cover numerical
injection, rejection without work loss, provenance, cross-parent grouping,
resplitting progress, early waiting merges, parent dependencies after a merge,
independent family release, location batches and simultaneous numerical workers.
Periodic triggers and admission boundaries are tested separately.
Spatial grouping tests cover normalized group diameter, cardinality, hull radius,
reduced constraint domains, semantic compatibility, unrelated parents, a dispersed
frontier and periodic compaction without repeatedly merging waiting replacements.

The backend supports forward evolution. Backtracking, numerical-configuration
alternatives, explicit recombinable-family keys, numerical radius-oriented
partitioning and concurrent merge construction remain future
work. General insertion/moving of genealogy nodes is outside the generic tree's
append/prune contract. The normal evolver return value remains a conventional
Orbit; exposing history alongside it is a separate API decision.

### Continuous enclosure recombination

`EnclosureRecombiner` is a `Handle<const EnclosureRecombinerInterface>`, owning a
polymorphic clone when constructed from an interface reference. Handle copies
share an immutable strategy; numerical work belongs to each invocation. The
interface transforms a collection into an upper covering collection without
changing its sources. `ReferenceEnclosureRecombiner` implements the existing
nonlinear-reference algorithm and produces at most one output. Grouping and
partitioning strategies can supply other concrete implementations later.

Each call receives an `EnclosureRecombinationContext`. Its optional maximum
enclosure radius is a finite, nonnegative per-run value supplied explicitly
when configuring the backend. It is not automatically copied from the evolver's
maximum_enclosure_radius and is not mutable strategy state.
An omitted value means unrestricted recombination. The result contains:

- `outputs`: records containing an enclosure and its original `source_indices`;
- `discarded_sources`: indices removed only after proving emptiness;
- `radius_limit_satisfied`: whether every output is certified within the request,
  measured by the existing `Enclosure::radius()`, or no limit was requested.

Indices refer to the original input list, before contraction or empty-source
removal. A future partitioning strategy may reference the same source from more
than one output. Provenance is recorded during construction; it is not inferred
later from bounding-box overlap. The result's joint upper-coverage contract is
independent of whether a radius request succeeds. The reference implementation
reports an unmet limit while retaining the covering output; it does not split
or discard nonempty sources to force compliance.

Empty batch input produces an empty result. All-proven-empty input produces no
outputs and reports all input indices as discarded. The existing free operation
`recombine(const List<Enclosure>&)` remains a single-output convenience wrapper
around the reference implementation, retaining its nonempty-list precondition
and returning a typed empty enclosure when no output survives. These operations
do not change trajectory identity or modify a backend. The tree backend uses
the independent grouping, admission and trigger policies described above. The existing
C++ benchmark still recombines post-reset checkpoint enclosures externally,
so its measurements do not quantify the integrated scheduler overhead.

Contract source copies with their own constraints and discard only proven-empty
ones. A singleton surviving source is returned without reparameterisation.
Require initialised models, matching state dimensions,
bounded reduced domains and structurally identical symbolic auxiliary mappings.
The first surviving source supplies the numerical configuration and reference.

Map a common unit parameter box onto each source's reduced domain with explicit
affine maps. Pad unused coordinates when parameter counts differ. This defines
the parameter correspondence; no physical identity between unrelated source
parameters is assumed. Compose state, absolute time, dwell time and every source
predicate with that source's own affine map. Let the resulting joint state/clock
maps be `F_i(u)` and take `P(u)=F_0(u)` as the nonlinear reference.

Bound each residual `F_i(u)-P(u)` on that source's feasible parameter set. Use
its predicates both to bound expressions directly by
`g(u) in C => f(u) in C+range(f-g)` and to contract the graph of the image with
`ConstraintSolver`. This uses constraints even when their coordinate-domain box
cannot shrink. The reference residual is zero. Hull the validated residual bounds
over the alternatives and append an independent parameter for each nonzero
component. The result is `P(u)+r`, with `r` in the residual box. There are at most
`state_dimension+2` new parameters, independent of the number of sources.

For each reference predicate, bound its values separately over every source's
feasible normalised parameters and hull those bounds. Retain that predicate with
the widened bounds, then contract the resulting enclosure. Source predicate
lists need not have matching lengths or ordering. They are never concatenated
as a conjunction across alternatives. Auxiliary variables are derived from the
common auxiliary state mapping, rather than independently hulled.

For any feasible source parameter, its normalised parameter remains permitted by
each widened reference predicate, and its state/clock difference from `P` lies
in the residual box. This establishes the upper inclusion. The construction
relaxes the disjunction and may lose correlations between residual components;
it does not claim an exact union or preservation of every state-space equality.
Results can depend on reference order and cannot be used as lower enclosures.
An oriented residual representation or additional transported predicates may
later improve tightness; they are not part of this implementation.

#### Constructive upper-coverage proof and numerical obligations

The following is a **semantic inclusion argument**, not a statement about
the inclusion of independently computed interval projections. It applies to
the *joint* state/absolute-time/dwell-time set represented by each input
`Enclosure` under upper semantics. Write a source as

\[
 S_i = \{ F_i(p) : p\in D_i,\ c_{ij}(p)\in C_{ij}\ (j=1,\ldots,q_i)\},
\]

where \(F_i=(x_i,t_i,d_i)\) includes all state coordinates and both clocks.
Each source's functions and predicates are validated function patches and
may include their own uniform error terms. The feasible parameter set is
interpreted in that validated enclosure semantics, not as the bare
parameter box \(D_i\).

**1. Exact affine parameter witnesses.** The implementation first contracts
each source copy and excludes it only if emptiness is *proved*. Let
\(D_i=\prod_k[\ell_{ik},h_{ik}]\) be its remaining bounded parameter box
and \(m=\max_i\dim D_i\). On \(U=[0,1]^m\), define

\[
 \phi_i(u)_k=\ell_{ik}(1-u_k)+h_{ik}u_k
 \quad (k<\dim D_i).
\]

Every original feasible parameter \(p\in D_i\) has at least one
\(u\in U\) with \(\phi_i(u)=p\), including degenerate intervals and
zero-dimensional boxes; the remaining coordinates are unused. This
is a *mathematical* affine map. Its implementation uses
`FloatDPBounds` coefficients to enclose the exact value.
The over-rounded model of \(\phi_i\) may report a range extending beyond
\(D_i\) at a boundary; that does not place the exact \(\phi_i(u)\)
outside \(D_i\). Consequently `unchecked_compose` is used
to avoid a false domain-error precondition, **with the obligation** that
validated model composition still encloses evaluations along the exact
admissible map. The unchecked operation itself does not prove that
obligation. Let \(\bar F_i(u)=F_i(\phi_i(u))\), and similarly
\(\bar c_{ij}(u)\).

**2. Sound residual and reference-predicate images.** Set
\(P(u)=\bar F_0(u)\), where index zero denotes the first surviving source.
For each source \(i\), use its *own* feasible normalised parameters
\(E_i=\{u\in U:\bar c_{ij}(u)\in C_{ij}\}\). The validated
`constrained_image_bounds` bounds the joint image of

\[
 H_i(u)=\big(\bar F_i(u)-P(u),\
                 \bar c_{01}(u),\ldots,\bar c_{0q_0}(u)\big),
 \qquad u\in E_i .
\]

For the reference, the first difference is identically zero. Direct
validated function-image bounds are intersected with implications
from each *source* constraint
\(\bar c_{ij}(u)\in C_{ij}\), using the identity
\(h(u)=\bar c_{ij}(u)+(h(u)-\bar c_{ij}(u))\).
A graph contractor then propagates the source constraints and
\(z=H_i(u)\); it may only remove points proved infeasible. The returned
box \(B_i\) therefore includes every value \(H_i(u)\) for feasible
\(u\in E_i\). Conversion via `cast_exact_box` copies the
outward-computed endpoints exactly; it does not narrow them.

**3. Explicit witness in the output.** Hull the validated boxes:
\(B=\operatorname{hull}_i B_i\). Let \(R\) be the first
\(\dim(F_i)\) coordinates of \(B\), and \(K_j\) its coordinates
for the reference predicates. The output is parametrised by

\[
 G(u,r)=P(u)+r,\quad
 u\in U,\quad r\in R,\quad
 \bar c_{0j}(u)\in K_j\quad (j=1,\ldots,q_0).
\]

Zero residual components are omitted only when their entire validated
range is exactly \([0,0]\); the remaining residual coordinates are
independent. For **any** \(v\in S_i\), choose one feasible source
parameter \(p\) representing it, a corresponding \(u\) with
\(\phi_i(u)=p\), and the residual
\(r=\bar F_i(u)-P(u)\) associated with that representation
(choosing compatible values of the validated patch models).
Since \(H_i(u)\in B_i\subseteq B\), this choice is admitted by
every output residual bound and every widened reference predicate.
It satisfies \(G(u,r)=v\) in *all* state and clock components at once.
The output's final validated `reduce()` contraction can
only remove infeasible output parameters, so it cannot remove this
feasible witness. Hence, conditional on the numerical obligations above,

\[
 \boxed{\ \bigcup_i S_i\ \subseteq\ S_{\mathrm{recombined}}\ } .
\]

The same state inclusion preserves any auxiliary values derived from
the **identical** symbolic auxiliary mapping required for all sources.
The radius request is advisory: a failed radius test does not truncate
the output.

**4. Scope and verification responsibilities.** The proof depends on
(a) validated source models and composition including the exact affine
maps, (b) outward image arithmetic and sound hull/intersection,
(c) sound graph/parameter-domain contractors and proved-empty
classification, and (d) exact copying of outward bounds. It does **not**
depend on projection-hull nesting, numerical speedup, guard equality
being satisfied by every *newly added* point, or source parameter maps
representing the same physical parameters. A different recombination
strategy requires its own argument; this proof covers the reference
upper-enclosure algorithm, not hybrid trajectory scheduling or lower
semantics.

The implementation uses `ConstraintSolver::reduce` for both
graph and final parameter contraction; its contract is validated
feasible-point preservation. `cast_exact_box` preserves
raw outward endpoints in `geometry/box.hpp`. Scaled patch
`unchecked_compose` delegates to validated Taylor-model
composition in `function/scaled_function_patch.tpl.hpp`;
its explicit domain compatibility check is skipped, and this is the
most important contract to guard with regression tests on rounded
non-dyadic parameter endpoints. The existing tests for different
constraints, source dimensions, clocks and guard equality exercise the
other steps. These checks support the implementation contracts but do
not replace their mathematical assumptions.

The C++ suite checks whole nonlinear curves through symbolic identities and
residual-domain inclusion, non-dyadic affine normalisation with differently
sized (including degenerate) source parameter domains and clock maps,
disjoint constrained sources, contraction without
coordinate-box shrinkage, differing predicate ordering and parameter dimensions,
clock dependence and subsequent nonlinear evolution, rounding and uniform errors,
point and empty sources, singleton/duplicate behaviour and auxiliary compatibility.
It also verifies interface dispatch with a separate implementation, polymorphic
clone lifetime and handle copying, original source indices after empty removal,
and radius reporting without coverage loss. Python exposes the interface, handle,
reference implementation, context and result records as well as the free operation.
Its checks cover source-list conversion, clone ownership, provenance, radius
requests and the distinct empty-list contracts of the batch and free operations.

### Representative recombination tests and benchmarks

Hybrid trajectory-level recombination is not integrated into the evolver.
The external C++ benchmark already recombines post-reset enclosures generated
by actual hybrid evolution, including guard, time and progress constraints,
instead of synthesising a list of translated boxes.
Correctness tests and performance measurements must use the same scenario
definitions. Recombination is initially an upper-semantics operation.

Both baseline scenarios have two locations, source and target, one transition
with guard `x+y >= 0`, and the identity reset `(x,y) -> (x,y)`. The vector fields
are different and nonlinear:

```text
source: dx/dt = 1,       dy/dt = x^2
target: dx/dt = 1+y^2,   dy/dt = 1
```

Both flows have curved trajectories. Different entry points on the guard give
different target trajectories: in the target, `d(x+y)/dt = 2+y^2 > 0`, so a
trajectory cannot meet the guard boundary twice. The identity reset preserves
entry geometry; it does not make the two vector fields interchangeable.

These nonlinear fields also admit explicit flow maps. From `(a,b)` in the
source, after elapsed time `t`,

```text
x(t) = a+t
y(t) = b+((a+t)^3-a^3)/3.
```

From an entry point `(p,q)` in the target, after elapsed time `delta`,

```text
x(delta) = p+(1+q^2)*delta+q*delta^2+delta^3/3
y(delta) = q+delta.
```

Use these formulas as independent mathematical references, evaluated with
validated arithmetic for numerical checks. They are references for verification,
not replacements for the production evolver when generating benchmark inputs.

#### Progressive urgent crossing

Start at absolute time zero with `a in [-2,-1]`, `b in [-1/2,1/2]` and make the
transition urgent. The whole initial set is strictly outside the guard. Along
the source flow, `d(x+y)/dt = 1+x^2 > 0`, so each initial point has a unique
crossing time `tau`. Writing `p=a+tau`, it satisfies

```text
p+p^3/3 = a^3/3-b.
```

The crossing equation is strictly monotone in `p`; it can be bounded independently
by validated scalar root isolation. All crossings occur by time `5/2`, since
the initial guard value is at least `-5/2` and its derivative is at least one.
A generation horizon of three therefore includes the entire crossing.

With integration step `h` much smaller than the spread of crossing times, each
step produces only the portion whose crossing belongs to that step. The guard
and timing constraints must cut the feasible source parameters accordingly.
Do not construct these portions manually or replace them by unconstrained boxes.
Immediately after the reset the reference coordinates are

```text
(x,y,time,dwell) = (p,-p,tau,0).
```

The zero dwell and the spatial guard equality must hold jointly with the crossing
time dependence. Contraction should be exercised on the generated constraints.
Some portions, particularly near the extremes of crossing time, permit strict
coordinate-domain reduction. Do not require every portion to shrink its box:
a nonlinear or oblique constraint can cut out points without changing any
coordinate projection.

#### Repeated permissive transition of the whole set

Start at absolute time zero with `a in [1,2]`, `b in [-1/2,1/2]`, and make the
same transition permissive. The entire initial set is enabled and remains
enabled along the source flow. The no-jump continuation remains in the source;
the whole population can also jump at every later time. Generate jump
enclosures over a jump-time window `tau in [0,1]`, including any immediate
transition at zero. Their reference coordinates immediately after reset are

```text
(x,y,time,dwell) = (a+tau,b+((a+tau)^3-a^3)/3,tau,0).
```

The numerical representation still has parameter and time constraints, but the
guard must not select a different proper subset of initial points in each time
window. This distinguishes repeated optional jumps from progressive urgent
crossing. Filter generated checkpoints by target location and transition
provenance; the source's final no-jump contribution is not a recombination input.

#### Collection, contraction and continued evolution

Collect actual post-reset checkpoints at trajectory birth, before target
evolution. A generation run may stop on the transition event and use its target
final enclosures as this dataset. Continuing the dataset is a separate run with
that event removed from the terminating-event criterion and with a sufficient
discrete-step budget. Do not clear event history merely to bypass termination.

For each dataset, retain an untouched source collection for comparison. Contract
each source with its own constraints before measuring reduced geometry or
constructing a common parameterisation. Constraints belonging to different
sources are alternatives; concatenating their lists would impose an incorrect
intersection. If a recombination relaxes their union, measure that relaxation
as part of the enclosure inflation.

Compare independent continuation of all sources with continuation of the
recombined enclosure to a common absolute horizon: for example, four for the
urgent case and two for the permissive case. For a checkpoint born at `tau`,
substitute `delta=T-tau` into the target flow above. Check the resulting state,
absolute time and dwell jointly. The original guard equality holds at urgent
entry, not at later target states; copying it as a current-state equality would
incorrectly cut away valid target evolution.

Correctness tests cover certified coverage of the selected source union,
constraints and clock dependencies, proven-empty sources, singleton inputs,
and continued-evolution coverage. Analytic bounds and selected reference points
provide useful regression checks; sampled points alone are not a certificate
of inclusion. The grouping distance is a separate criterion and cannot certify
coverage. Derive centres and widths from maps evaluated over reduced domains:
the current `Enclosure::centre()` uses the unreduced state bounding box and does
not automatically account for `reduce()`.

Benchmark the following phases separately:

1. Hybrid generation of the source dataset, outside recombination timing.
2. Source constraint contraction, starting from fresh copies for each sample.
3. Construction of the recombined enclosure and its final contraction.
4. Continued evolution, both independently and after recombination.

Also report the total contraction-plus-construction cost; otherwise an algorithm
could appear faster by moving work into dataset preparation. Verification and
diagnostic reporting remain outside the measured phases. Do not hardwire logging
redirection into the benchmark or tests.

For window-width comparisons, fix the maximum step and window midpoint
and vary the physical birth-time interval width. Record the selected
checkpoint count and the costs of contraction, recombination and both
continuations separately. Algebraic whole-domain witnesses on
guard-constrained sources exercise joint state/time coverage, while
per-source validated analytic overlap checks on generated sources
remain necessary rather than sufficient for general inclusion.

For step-resolution comparisons, fix an absolute post-reset birth-time window,
impose it as a constraint on each generated enclosure, and include every
feasible checkpoint intersecting that window. The source count may then
change with the step; do not also fix group cardinality in this experiment.
Use the separate adjacent-cardinality sweep to measure the grouping effect.

Start with maximum integration steps `1/16`, `1/64` and `1/256`; record the actual
accepted steps and generated checkpoint count rather than assuming fixed steps.
Vary initial-set extent and local error tolerance independently. Avoid incidental
radius-driven subdivision in the baseline, and record any subdivision or
reconditioning that does occur. Subsequent cases can deliberately exercise them.
Use compatible nearby groups with fixed distance scales; do not force the entire
dataset into one enclosure to obtain a desired reduction factor.

Alongside timing, record source/group/output counts, parameters, constraints,
domain contraction, maximum normalised centre distance, validated state/time
bounds, and the cost and enclosure growth of continued evolution. Volumes alone
are unsuitable for urgent post-reset sets, which lie on a guard surface. Neither
small runtime nor a reduced enclosure count is sufficient if uncertainty growth
makes subsequent evolution more expensive or less informative.
