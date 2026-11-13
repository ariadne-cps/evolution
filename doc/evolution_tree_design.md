# Evolution tree design proposal

Status: the flat and serial trajectory-tree backends are implemented.
Continuous enclosure recombination is implemented; hybrid trajectory
recombination, backtracking and per-work numerical alternatives remain design proposals.

## Starting functions

1. Each node represents one trajectory containing an ordered list of sets.
2. Parent-child relations distinguish trajectories born from another trajectory.
3. Each node records whether its own trajectory is complete or still needs
   evolution.

The tree describes genealogy; the list inside a node describes progression of
that trajectory. These requirements replace the previous checkpoint-per-node
proposal. A new set does not by itself create a new tree node.

The generic structure belongs in
`ariadne-cps/tree/doc/derivation_tree_design.md`. Entry semantics, trajectory birth,
completion criteria and scheduling belong to Evolution.

## Trajectory object and metadata

```text
Trajectory
    ordered list of EvolutionEntry

EvolutionEntry
    computational checkpoint, or explicitly grouped numerical pieces
    Orbit contributions tagged reach / intermediate / final
    shared reference to const NumericalConfiguration
    information about the step producing the entry, when applicable

TrajectoryMetadata
    completion: NEEDS_EVOLUTION or COMPLETE
    origin entry in the parent trajectory, except at the root
    reason for birth
    event / branch information, when applicable
    optional application-specific labels
```

The generic owner supplies stable node identity, parent and children; those
links are not duplicated in Evolution metadata. A node ID can identify its
trajectory initially. A separate trajectory-ID namespace is unnecessary merely
to inspect genealogy; introduce one later only with an explicit relationship to
node identity.

Entries follow the logical progression of the trajectory, not a sort by time-
interval endpoints. A discrete step may preserve continuous time, and validated
time intervals may overlap. An entry must distinguish the checkpoint used for
further computation from the sets contributed to the Orbit; these are not
necessarily the same enclosure. The backend proposal below records both.

Keep configuration with the entry to which it applies. Unchanged configurations
can be shared immutably. Later changes must not alter earlier entries. Examples
include integration step limits, local error targets, parameter-count limits
and numerical discretization/reconditioning choices.

## Continuing versus creating a trajectory

| Situation | Effect |
| --- | --- |
| Pure continuous step | Append the resulting entry to the current node |
| Fully resolved transition with one continuation and unchanged trajectory identity | Append the transition result to that same node |
| A distinct trajectory is generated | Create a child with its own ordered list |
| A parent continues after a child's birth | Append to the parent; retain the child's original source reference |
| Another distinct trajectory is generated later | Add another child referring to its own source entry |

For one permissive/unresolved event, the no-jump continuation stays in the
original node and a distinct jump trajectory is a child. Creating that child
does not complete the parent. An internal tree node may still need evolution;
active trajectories therefore cannot be identified just by taking tree leaves.

A child starts with its own initial set, for example the reset image of the
relevant source portion. It need not copy the parent's entire list. Earlier
provenance is available through the parent and origin reference, and the child's
first set need not equal its parent's source set.

Record the particular parent entry responsible for birth, in addition to the
reason/event. Parent identity alone is insufficient because its list can keep
growing. An entry index is sufficient for an append-only initial history.
Backtracking may truncate and regrow history, but may not insert entries between
existing entries or move nodes. References to truncated history need explicit
invalidation, for example through a history generation; an index alone must not
silently resolve to a replacement entry at the same position.

## Completion is local to a trajectory

`NEEDS_EVOLUTION` means more work remains for the trajectory stored in that node.
`COMPLETE` means that local evolution has finished under the current termination
criteria: it will generate no more sets or child trajectories in normal forward
processing. This can follow an evolution horizon, exhausted/empty continuation,
or another explicit stopping condition.

Completion does not mean that descendants are complete. A complete parent can
have children that have not yet been evolved. Its list remains available as
history. Explicit future backtracking may reopen work, but that requires a
separate validity policy rather than silently appending to a completed history.

Do not publish completion until all results and generated children from that
trajectory's local work have been committed. Otherwise the child count observed
at completion would be premature. Worker/scheduler state, such as whether a
trajectory is currently running, is distinct from trajectory completion.

## Finish a parent before descending

The proposed initial scheduling policy favors completing an ancestor trajectory
before evolving the trajectories generated from it:

1. Select a node that needs evolution.
2. Continue its trajectory, appending entries and recording generated children.
3. When local work is finished, mark the node complete.
4. Inspect all children it generated; their count and origins are now known.
5. Decide whether to keep those child trajectories or propose recombination into
   one or more successor trajectories.
6. Descend to the resulting trajectories that still need evolution.

This constrains parent-versus-descendant processing; it does not yet choose a
global breadth-first or depth-first order among unrelated subtrees. The generic
tree stores the necessary information without implementing this scheduler.

This completion boundary is useful because children can be collected before
their own evolution is scheduled. It makes sibling-group recombination a first
policy to investigate, without requiring an initial algorithm for merging
already evolved arbitrary branches.

## Example: repeated guard crossings

```text
A: [A0, A1, A2, A3, ...]     COMPLETE
    +-- B: [B0]   origin A1, event e; NEEDS_EVOLUTION
    +-- C: [C0]   origin A2, event e; NEEDS_EVOLUTION
    +-- D: [D0]   origin A3, event e; NEEDS_EVOLUTION
```

A generated B, C and D while continuing its own list. Once A is complete, the
scheduler knows all three children before descending. It can inspect candidate
groups such as {B,C,D} or {B,C} and {D}, and decide whether a sound replacement
by one or more trajectories is useful.

As clarified, different portions crossing a guard during different integration
steps remain distinct trajectories, potentially in a common recombinable family.
Common parent, reason/event and source entries provide initial provenance; they
do not establish mathematical compatibility by themselves. A precise family
key and recombination criterion remain later design work.

## Recombination constraints and deferred operations

Parent completion establishes that the generated group is known; it does not
prove that its children overlap or may be recombined. A future policy must check
location, clocks, discrete context, numerical-alternative compatibility and
sound coverage of all selected source sets. Bounding-box overlap or common
ancestry is insufficient.

Recombination can produce one or more replacement nodes/trajectories. Preserve
which original children contributed to each replacement, rather than deleting
their origin information. Distinguish the collection originally generated by
the complete parent from the active collection after recombination. Structural
placement of replacements and absorbed-node handling are not finalized here.

The earlier enumeration of `SPLIT`, `UNBUNDLE`, `RECOMBINE` and
`ALTERNATIVE_CONFIGURATIONS` remains context, but its checkpoint-per-node
interpretation is superseded. Before adding those operations, specify how
subset and numerical-alternative branches relate to the trajectory list. Do
not flatten parallel pieces into an apparently sequential history. Additional
internal representation or explicit relations may be needed.

Stable IDs remain required: appending entries or adding children does not
renumber existing nodes. Backtracking can still invalidate results derived from
edited entries. Persistence, concurrent publication, general multi-input
dependency graphs and reclamation are outside the initial scope.

## Initial operations and verification

The wrapper needs creation of an initial trajectory, ordered append, creation
of a derived trajectory at an identified source entry, completion marking/query,
and access to lists and genealogy. The generic owner does not decide when a
trajectory is born or complete.

Test that continuous steps extend one list without increasing node count;
configuration history is retained; children identify the right parent and
source entry; a parent can continue after generating children; completion is
local and published after all generated children; and a complete parent can
have unevolved children. A scheduler test should confirm that descendants are
deferred until their ancestor is complete and the child-group decision has
been made. Recombination selection/construction is subsequent work.

## Hybrid evolver backends

### Implemented first stage

`include/evolution/hybrid_evolution_backend.hpp` defines work phases, successor
provenance, explicit Orbit contributions, local stop reasons and the backend and
factory interfaces. `FlatHybridEvolutionBackendFactory` is the default.
`HybridEvolverBase::set_evolution_backend_factory()` clones factory settings;
each call to `orbit()` creates independent run state, invokes the same numerical
callback, and obtains the conventional Orbit through `release_orbit()`.

The numerical callback returns `HybridEvolutionStepResult`. The five previously
coupled methods now write only to that local result; they neither dispatch work
nor publish to a shared Orbit. The flat backend reuses `DynamicWorkload`, commits
all contributions under one lock before scheduling successors, and rejects
release before successful completion or a second release. Exceptions propagate
through the executor after its outstanding work joins.

Work items currently contain the enclosure (including its existing enclosure
configuration) and processing phase. The evolver's integration/solver settings
continue to use the existing uniform configuration for this first stage;
per-work immutable numerical settings and alternatives are not implemented.
Successors record event, event kind, optional crossing kind and source-coverage
evidence. Whole-source resolution is never inferred from the successor count.

`tests/test_hybrid_evolution_backend.cpp` verifies output/checkpoint separation,
publication before successor dispatch in serial and concurrent modes, lifecycle
and failure handling, and transition/subdivision provenance. An independent LIFO
backend in that test is selected through the public factory setter and compares
all three output collections with the default backend for urgent and permissive
systems, including repeated calls after the original factory is destroyed.
The existing hybrid evolution tests continue to exercise the numerical behavior.

### Implemented serial trajectory-tree backend

`TreeHybridEvolutionBackendFactory` selects the tree backend through the same
factory setter and callback/result contract. The ordinary `orbit()` return type
is unchanged. A concrete `TreeHybridEvolutionBackend` can also be run directly;
`tree()` provides const inspection before/after execution and after Orbit
release. The standard `orbit()` call owns its backend only for the duration of
the run; returning the detailed history from that API remains separate work.

Evolution directly pins Tree at `e1372f9c7d5f1a6f1db11dc74b58053e1d8ad9e7`.
Its CMake target is established before transitive Tree consumers, so all compiled
consumers and the installed header bundle use that public API. Compared with
the older transitive Tree revision, BinaryTreeNode changes are documentation
only. Foundation, Threading and Build pins agree with the existing dependencies.

The representation is:

```text
DerivationTree<HybridEvolutionTrajectory, HybridEvolutionTrajectoryMetadata>
    trajectory.entries: ordered vector of HybridEvolutionEntry
        work: input checkpoint and processing phase
        origin: optional same-trajectory (entry index, successor ordinal)
        result: optional complete HybridEvolutionStepResult
    metadata.complete: local trajectory completion
    metadata.origin: optional parent (entry index, successor ordinal)
```

The source result retains successor reason/event information, so an origin can
resolve it without duplicating parent links or event fields. New entries are
appended. A missing entry result denotes an unprocessed/failed computation.
Subdivision pieces have a common predecessor and distinct successor ordinals;
they remain parallel pieces, not successive time steps just because their
records are adjacent. Their continuations and reconditioned representations
retain the same-node predecessor link. This does not yet assign independent
subset/alternative ID namespaces.

The scheduler drains every local piece of the current trajectory before marking
that node complete and releasing its children. Unrelated ready trajectories are
processed breadth first. Children generated at different source entries remain
distinct trajectories, including when they share parent/event provenance. No
recombination is performed. The completed-parent boundary is the insertion
point for a later domain-specific recombination policy.

A reset normally derives a child. It continues the same node only when it is
the sole successor, is urgent/impact, carries a FULL_SOURCE certificate, and
the trajectory has not been divided by subdivision or prior urgent branches.
Once portions have separated, even the last fully crossing piece becomes a
child: it must remain distinct from pieces crossing at earlier steps, and all
these children must be known before their evolution begins.
The numerical kernel currently certifies strictly active immediate urgent
transitions by proving the invariant-valid source lies wholly inside the guard.
Other cases remain UNRESOLVED, including permissive events and flow crossings
without an explicit coverage proof. Transversality or an empty other branch
alone is insufficient. The flat backend ignores this genealogy classification
and preserves the same numerical output contributions.

No global Orbit accumulator duplicates the tree history. Statistics count
published contributions; `release_orbit()` gathers their explicit reach,
intermediate and final roles after the complete run. It never turns a leaf or
an empty completed trajectory into an implicit final set. Output collections
may have a different traversal order from flat execution, but preserve the
same enclosures without recombination. History remains inspectable after the
one permitted Orbit release.

A callback exception leaves incomplete metadata and unprocessed entries; it
propagates and prohibits Orbit release or another run. The backend is serial,
including when ThreadManager is configured for parallel work. Editing history,
pruning/reopening a run, numerical alternatives and parallel tree publication
are not exposed by this initial implementation.

`tests/test_hybrid_evolution_tree.cpp` covers parent-first dispatch, distinct
children from repeated event births, explicit origins, whole-source versus
unresolved/permissive identity, parallel split/reconditioned pieces, exhausted
empty trajectories, failure histories and exact flat/tree output comparisons
for upper and lower semantics.

The remaining subsections retain the architectural rationale and deferred
operations; the implemented stages above define the current API.

### Coupling before the refactoring

The reference implementation examined is Evolution commit
`f34ff6e7e72abe2c2527bf7bc308cb4893010de9`.
`HybridEvolverBase` at that revision owns both the execution mechanism and result
publication:

- `orbit()` constructs a `DynamicWorkload` and a `SynchronisedOrbit`.
- A work item is `Pair<HybridEnclosure, Bool>`. The Boolean requests starting
  event/termination processing; it is true for the original input as well as
  after a reset, so it is not a trajectory identity or a completion flag.
- `_process_working_set`, `_process_jumped_set`, `_process_starting_events`,
  `_evolution_step` and `_apply_evolution_step` pass the workload and Orbit
  accumulator through numerical processing.
- Successors are submitted with `workload.append`, and results separately with
  `adjoin_reach`, `adjoin_intermediate` and `adjoin_final`.

Changing just the queue element to a tree node ID would leave the evolver
responsible for scheduling and flattening. It would also lose information
needed to distinguish continuation, transition and subdivision.

### Proposed boundary

Keep `HybridEvolverInterface::orbit(...) -> Orbit<HybridEnclosure>` unchanged.
Introduce a run-owned `HybridEvolutionBackendInterface` with a flat backend and
a tree backend. The evolver stores a factory/prototype selecting the backend;
each `orbit()` call creates independent mutable execution state. Backend state
must not be shared accidentally between calls on a const evolver.

The numerical evolver computes one local outcome. The backend owns pending work,
selection of the next work, history, completion publication, and production of
the conventional Orbit. It may later use a recombination policy at a defined
completion boundary. Numerical guards, resets, invariant restrictions,
termination predicates and upper/lower semantics remain in the evolver.

Conceptually, both implementations use the same calls:

```cpp
auto backend = backend_factory.create(initial);
backend->run([&](const HybridEvolutionWorkItem& item) {
    return _advance(item, termination, semantics);
});
return backend->release_orbit();
```

`run` owns dispatch, including any concurrent execution and joining before
returning. `_advance` must neither enqueue work nor mutate a shared Orbit.
It can accumulate a step-local result through an internal collector passed to
the existing numerical helpers. Publishing the result is a backend operation
after the callback has returned successfully.

Move step-local results into backend storage where the enclosure/container APIs
permit it. The abstraction does not require copying a parent's entire history
when scheduling one work item or creating a child trajectory.

Use a runtime interface at this boundary, rather than templating
`HybridEvolverBase`, its derived evolvers and its public interface on the
backend. Dispatch occurs per work item, not in scalar numerical operations.
The generic `DerivationTree<Object, Metadata>` remains unchanged, with its two
template parameters and no hybrid scheduling logic.

### Work item and step result

An initial common work item contains:

```text
HybridEvolutionWorkItem
    enclosure to process
    phase: PROCESS_STARTING_EVENTS or EVOLVE_CONTINUOUSLY
    shared immutable numerical configuration
```

Scheduling tickets and tree node/entry IDs belong to the backend. The numerical
callback does not need to know which concrete container owns its work item.
An entry-processing phase is independent of the reason a trajectory was born.

A step result must contain:

```text
HybridEvolutionStepResult
    Orbit contributions, explicitly tagged REACH / INTERMEDIATE / FINAL
    successors, each with enclosure, processing phase and configuration
    successor provenance: continuation, transition, subdivision,
        reconditioning, or numerical alternative
    event/crossing/resolution information where applicable
    local stop reason if no continuation remains
```

Keep associations between contributions and successor pieces where a guard
creates several pieces. A local result is a record of one computation, not a
claim that every collection in it is one sequential trajectory.

Do not infer trajectory birth merely from the number of successors. A numerical
subdivision and competing transitions have different meanings even when both
produce two sets. A permissive transition retains the no-jump continuation and
derives a jump trajectory. A resolved single transition can continue the same
trajectory only when that classification is justified for the whole affected
set. `CrossingKind::TRANSVERSE` alone is not that justification. Different source
portions crossing in different steps remain distinct derived trajectories in a
potentially common recombinable family.

Orbit contributions and successors must remain separate. In the current
`_apply_evolution_step`, `evolve_sets` can have urgent-guard/time constraints
which are intentionally absent from `next_working_sets`. An intermediate set
therefore cannot simply be reused as the next computational checkpoint.
Similarly, a reach set is a swept set, not an endpoint. Starting events can
publish a pre-reset reach set and a post-reset intermediate set together.

### Completion and scheduling

The flat backend initially reproduces the existing workload behavior and
accumulates the tagged contributions directly into an ordinary Orbit. The tree
backend records results in trajectory entries and can keep newly derived
children blocked while their parent is still generating continuations.

Returning from one `_advance` invocation completes one work item, not necessarily
its trajectory. A trajectory becomes COMPLETE only after all its continuations,
including numerical pieces belonging to it, have been exhausted, all in-flight
work has returned, and all its generated children and contributions have been
published. Only then may the backend inspect the complete child group, apply a
future recombination decision, and release children for evolution.

Track computational status (ready/running/blocked/retired) separately from
trajectory completion. Numerical subdivision may require several pending pieces
inside one trajectory; those pieces must not be appended as if they were
successive time steps. If subdivisions/alternatives are enabled, trajectory
entries must explicitly represent grouped pieces and their predecessor links.
Detailed subset/alternative identity rules remain separate design work and are
not solved by silently turning every work item into a trajectory node.

The backend may reuse DynamicWorkload as its executor. It must submit only work
which its own policy has made eligible: merely choosing a queue order cannot
prevent a concurrent executor from starting a child too early. State publication
and ready-work selection require synchronization. Numerical processing happens
outside that synchronization, with worker-owned work snapshots; no mutable
DerivationTree references or children spans are held by workers across commits.
Begin the tree backend serially before enabling parallel commit/dispatch.

An empty executor queue is not a complete run if the backend still has blocked
or in-flight work. Completion of a run requires all applicable work to be
finished and no deferred child-release/recombination decision to remain.
Exceptions must join outstanding work and propagate rather than return a
partially computed Orbit as a complete result.

### Orbit compatibility

Recombination is not a prerequisite for returning an Orbit from the tree backend.
Every trajectory can be evolved independently and its tagged contributions
projected into the conventional collections. With recombination absent, the flat
and tree backends must produce exactly equal cardinalities separately for reach,
intermediate and final, preserving multiplicity as well as the original initial
enclosure. Tree-node count, completion and genealogy do not replace these output
roles. The scheduling order may differ; the classified multiset of enclosure
contributions must agree.

The regression in `tests/test_hybrid_evolution_tree.cpp` checks these three
cardinalities explicitly, compares each full collection without removing
duplicates, and checks both initial enclosures against the original input. It
covers urgent/permissive events, upper/lower semantics, continuous and discrete
limits, initially active jumps, terminating events, subdivision and invariant
restrictions including excluded inputs.

Preserve all four existing roles: `initial()`, `reach()`, `intermediate()` and
`final()`. The initial enclosure is the original input, not a recombined root or
the collection of all trajectory seeds.

Tree leaves are not final sets, and COMPLETE does not imply a final contribution.
For example, an empty continuation completes local evolution without adding a
final set. Conversely, the existing time-straddling logic may emit a final slice
while retaining further work. Preserve explicit classifications made by the
evolver; neither topology nor a generic terminal flag can reconstruct them.

For the tree backend, materialize Orbit from valid tagged contributions at the
end. Avoid an append-only global Orbit as a second authoritative representation:
pruning/backtracking would otherwise leave obsolete contributions in it. Logical
retirement of pending recombination sources prevents further evolution, but does
not by itself invalidate previously validated history. Actual backtracking must
invalidate contributions derived from the truncated history as well as work
tickets; pruning generic nodes alone does neither.

Without recombination or different numerical choices, flat/tree results should
match the same classified enclosure contributions, disregarding scheduling order
where concurrent ordering was already unspecified. With recombination enabled,
the public result contract is preserved, but enclosure counts and partitions
may legitimately differ. Soundness must be established for the requested
semantics: an enclosing hull may be acceptable for an upper approximation but
is not a generic sound recombination for lower semantics.

`release_orbit()` moves the flat accumulator or gathers the tree's valid results.
Keeping a queryable tree after `orbit()` returns would need a separate detailed
result/session API; changing the existing return type is not required for this
first refactoring.

### Configuration and recombination ownership

Keep numerical configuration with the work/entry it governs, sharing immutable
snapshots when unchanged. Future alternatives must read the item's configuration,
not mutate the evolver's shared global configuration. Integrator/solver state
must follow the same rule; assigning a configuration pointer to a work item is
insufficient if helpers still read mutable evolver state. Initially all work may
use one per-run snapshot, preserving current numerical behavior.

The backend identifies candidate groups through provenance and controls when a
decision may run. A domain-specific policy must validate and construct any
recombined enclosure; the generic Tree container does not make this numerical
decision. Parent completion alone is not a recombination certificate.

Replacements are newly added trajectories with explicit source provenance.
Pending absorbed sources may be logically retired in Evolution metadata.
There is no reparenting or insertion between existing generic nodes.
`DerivationTree::prune(id)` removes every strict descendant, so it cannot erase
an arbitrary selected subset of sibling sources; do not assume that operation
exists. Preserve provenance before invalidating any referenced history.

### Continuous enclosure recombination

`recombine(const List<Enclosure>&)` constructs an upper enclosure of a nonempty
source list without changing its sources. It does not select nearby groups,
change trajectory identity, or modify a backend. Hybrid recombination is a later
extension, and the benchmark is deferred until hybrid traces can be generated.

Contract source copies with their own constraints and discard only proven-empty
ones. A singleton surviving source is returned without reparameterisation; an
all-empty list returns an empty source copy. Require matching state dimensions,
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

The C++ suite checks whole nonlinear curves through symbolic identities and
residual-domain inclusion, disjoint constrained sources, contraction without
coordinate-box shrinkage, differing predicate ordering and parameter dimensions,
clock dependence and subsequent nonlinear evolution, rounding and uniform errors,
point and empty sources, singleton/duplicate behaviour and auxiliary compatibility.
Python exposes the same operation and exercises actual source-list conversion
and rejection of an empty list.

### Representative recombination tests and benchmarks

Hybrid recombination and its benchmarks remain to be implemented. Use
enclosures generated by actual hybrid evolution, including its guard, time and
progress constraints, rather than synthesising a list of translated boxes.
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

### Incremental migration and verification

1. Replace the anonymous Boolean with a named processing phase and introduce
   the step-local result. Route all current result/queue writes through this
   record, initially retaining the existing executor and numerical choices.
2. Move execution and accumulation into the flat backend. Adapt the five coupled
   methods and make `_log_summary` read backend counters rather than an Orbit
   accumulator. Keep `orbit()` overloads and default behavior unchanged.
3. Add the serial tree backend with the same callback/result contract, append
   history, explicit trajectory birth and completion, and parent-first dispatch.
   Initially disable recombination so numerical behavior can be compared.
4. Specify subset/alternative records before claiming support for those tree
   cases; add completion accounting for all pieces. Introduce recombination at
   the completed-parent boundary only after its semantic checks are specified.
5. Add parallel tree dispatch and backtracking only with synchronized publication,
   stale-ticket rejection, configuration isolation and invalidation of derived
   contributions. These are not prerequisites for the initial backend boundary.

Verify the flat backend against existing orbit tests, then both backends on the
same systems with recombination disabled. Cover immediate/permissive/urgent
events, competing jumps, repeated guard crossings across steps, empty exhausted
continuations, invariant rejection, final-time slices, discrete termination,
reconditioning and subdivision. Check classified coverage as well as counts;
use deterministic serial runs where exact collection comparison is appropriate.
Scheduler tests must show that a child is recorded before its parent completes
but is not dispatched until the completed child group has been handled, and that
one completed split piece cannot prematurely complete its trajectory.
