# Evolution tree design proposal

Status: design draft; no implementation or finalized public API is claimed.

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
    set / enclosure
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
time intervals may overlap. Before implementation, settle whether an entry
stores an evolved checkpoint, a reach segment, or a record containing both.
This does not change the node-as-trajectory contract.

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
growing. An entry index is sufficient for append-only progression. Backtracking
prunes a future continuation and regrows from the retained checkpoint; it does
not insert intermediate nodes or move existing subtrees. Truncating trajectory
entries must also invalidate child trajectories derived from removed entries.
The generic `prune(id)` boundary retains that node and removes its descendants.
Lifetime and truncation of source-entry references remain consumer concerns.

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

The structural model permits growth and pruning followed by regrowth. Existing
nodes cannot be inserted between other nodes or reparented. Pruned IDs must
become invalid and must not be reused; retained nodes keep their identities.
The generic container implements `prune(id)` by retaining that checkpoint and
removing all strict descendants, ready for regrowth through `add_child`.
It preserves the retained payload and metadata; truncating the trajectory list
and reopening a completion flag remain Evolution responsibilities. Persistence, concurrent publication and general
multi-input dependency graphs remain outside the initial scope.

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
