# Evolution tree design proposal

Status: design draft. This document records the requested hybrid-evolution
semantics and proposed use of the generic structure; it does not describe an
implemented `EvolutionTree`.

## Motivation and document ownership

A wide initial enclosure can meet an urgent guard over many integration steps.
The evolver then creates many post-transition enclosures from different portions
of the source bundle. We need provenance that identifies related branches and
supports controlled recombination of active enclosures.

The hybrid-specific design lives here. The generic owner, stable node handles,
derivations, structural parent relation and multi-input dependency history are
specified in `ariadne-cps/tree/doc/derivationtree_design.md`. Evolution consumes
Tree; Tree must not acquire dependencies on hybrid systems or enclosures.

An important clarified requirement is that portions crossing the same urgent
guard during different integration steps remain distinct trajectories, but may
belong to a common recombinable family. Equality of trajectory IDs must therefore
not be the only way to find recombination candidates.

## Node object, metadata and configuration

An evolution node represents a checkpoint from which work may be continued.
Its object contains the enclosure/work state and a shared immutable numerical
configuration. Store the configuration once; the enclosure and its configuration
must not have two competing authoritative copies.

Conceptually:

```text
EvolutionObject
    enclosure / work state
    shared reference to const NumericalConfiguration

EvolutionMetadata
    TrajectoryId
    SubsetId
    AlternativeId
    TrajectoryFamilyId
    incoming branch description, when needed
```

The generic owner supplies `NodeId`, parent, children, derivation references and
lifecycle state. These are not duplicated in `EvolutionMetadata`. A derivation's
metadata supplies the reason for its outputs and the events/decisions associated
with the operation. The input and output records together describe a step;
storing checkpoints does not replace storage of the continuous reach segments
needed for the eventual orbit.

Configuration examples include integration step limits, local error targets,
function-domain parameter limits, and numerical discretization/reconditioning
choices. Distinguish planned numerical operations from actual automaton events.
An unchanged configuration is inherited by sharing its immutable object. A
changed configuration produces a new object while existing checkpoints keep the
old one. Global mutable evolver settings cannot silently override a node's
recorded configuration.

## Meaning and stability of identifiers

| Identifier | Meaning |
| --- | --- |
| Node | One immutable historical checkpoint, unique within its tree |
| Trajectory | A computational evolution strand that can continue across steps |
| Subset | A label for subdivision of a strand's represented set |
| Alternative | A label for a numerical configuration alternative |
| Trajectory family | A compatibility/provenance group containing potentially recombinable strands |

Allocate fresh identifiers independently of tree depth, sibling position and
storage addresses. Never reuse them after backtracking. Strongly typed ID
wrappers should prevent mixing the different namespaces.

When an operation creates several outputs, designate one as the identity-
retaining continuation. It inherits the parent's applicable trajectory, subset
and alternative IDs. Other outputs receive fresh values only in the dimensions
changed by the operation. Every output still receives a fresh `NodeId`.
Selection of the retaining continuation is explicit and stable; it is not
inferred from whichever child happens to be visited first.

IDs are labels, not a complete mathematical description of ancestry or coverage.
In particular, retaining a `SubsetId` through repeated splits does not mean that
the later subset equals the earlier one. Derivation provenance records that
relationship. Stable labels avoid renumbering after insertion; they do not
eliminate invalidation when an edit changes the meaning of later computations.

## Operation semantics

This table resolves the requested distinction between a single fully resolved
successor and branching successors. Exact C++ enum names remain provisional.

| Reason | Outputs | Trajectory IDs | Subset/alternative IDs |
| --- | --- | --- | --- |
| `CONTINUOUS_STEP` | One | Inherit | Inherit |
| `URGENT_TRANSITION` / `DISCRETE_STEP`, fully resolved, one successor | One | Inherit | Inherit |
| Fully resolved step with several possible discrete successors | Several | Retaining successor inherits; others fresh | Inherit unless independently refined |
| `SPLIT` | Two or more | Inherit for every child | Retaining child inherits subset; others fresh; alternatives unchanged |
| `UNBUNDLE` | One | Fresh | Preserve provenance; no time advance |
| `RECOMBINE` | One from several active sources | Retain primary source's ID as proposed policy | Retain primary labels and record all absorbed coverage |
| `PERMISSIVE_TRANSITION` / `UNRESOLVED_TRANSITION` | No-jump plus jump successor for one event | No-jump inherits; jump receives fresh ID | Inherit unless independently refined |
| `ALTERNATIVE_CONFIGURATIONS` | Several | Inherit | Subset unchanged; retaining alternative inherits, others fresh |

`SPLIT`, `UNBUNDLE` and `ALTERNATIVE_CONFIGURATIONS` do not advance continuous
or discrete time. A continuous step advances continuous time only. A resolved
transition records its actual continuous-time progress and discrete event;
`DISCRETE_STEP` can have zero continuous-time progress.

Splitting can be recursive: a split child can itself be split until a requested
criterion holds. Those are explicit derivations, not a flattened list that loses
the splitting history. A split must soundly cover its source; children need not
be geometrically disjoint when they are validated enclosures.

For `UNBUNDLE`, the proposed history representation creates a unary successor
with a new trajectory ID, rather than mutating the old checkpoint. This is an
intentional exception to ordinary unary identity inheritance. Genealogical family
membership may be retained, but must not automatically reverse a deliberate
decision to treat a subset independently.

The two-output transition case describes one unresolved/permissive event. Several
simultaneously possible events require one no-jump continuation where admissible
and a jump output per admissible event, or explicit staged derivations preserving
the same alternatives. Do not impose a global binary-child restriction.

## Families for repeated guard crossings

Different integration slices of a wide source bundle can generate distinct jump
trajectories. They must remain distinguishable while being discoverable as a
family. Do not use a fresh source `NodeId` or integration-step number as the sole
family key, since that would separate the very slices we want to collect.

A proposed family registry uses the source semantic family, the selected event
and target location, and relevant discrete-decision context to locate a target
family. Repeated slices of that same transition can refer to the same target
family while receiving fresh trajectory IDs. Distinct event histories, including
additional loop iterations, are not automatically identified. Configuration
alternative compatibility is checked separately.

A no-jump continuation keeps its pre-transition family. A resolved single
successor can retain its trajectory ID while updating its family to reflect the
new discrete context. A family is a candidate-search group, not a certificate
that every pair of its members may be merged.

The exact family-key representation and the default effect of `UNBUNDLE` on
automatic candidate selection remain design decisions. The requirement to keep
different slice trajectories distinct is fixed.

## Recombination eligibility and soundness

Search among active frontier nodes, including nodes in different structural
subtrees. Do not restrict search to siblings or equal trajectory IDs. Start with
family membership, then check the hybrid and numerical context.

Before accepting a group, require:

1. Compatible location, state space and discrete evolution context.
2. Compatible continuous-time representation, pending event obligations and
   termination/discrete-step budgets.
3. An explicitly compatible numerical-alternative context.
4. A sound result enclosure for the complete union of selected sources.
5. No concurrently committed continuation of a consumed frontier source.

Equal family IDs or intersecting bounding boxes are insufficient. Bounding-box
overlap is only a cheap candidate filter; it may be caused by lost correlations.
Geometric overlap is not itself necessary for a sound outer union, although it
may be useful for deciding whether recombination improves cost or precision.

If inputs are synchronized at a common time, the result must contain their union
at that time. If they carry different time ranges, either evolve them to a valid
common checkpoint or construct an enclosure of the joint time-state sets that
preserves a sound time-state relationship. Do not erase the clock or associate
every input state with an arbitrary new common time. An inclusion argument must
cover every source and its still-relevant constraints.

For upper reachability semantics the essential contract is

```text
union(represented_time_state_set(source_i))
    subset_of represented_time_state_set(result).
```

No lower-semantics guarantee follows from taking an outer hull. Recombination
must declare which semantics its enclosure construction supports.

The proposed result retains the primary source's trajectory label even when
other sources have different trajectory IDs. All source IDs and coverage
provenance remain recorded. The output represents their combined work; its
retained label is not a claim that the other trajectories never existed.

## Subsets versus numerical alternatives

Split children jointly represent coverage that must be retained until proved
empty, terminated or subsumed by a sound replacement. Numerical alternatives
are different computations of the same obligation: choosing one can supersede
the others if that computation soundly covers the obligation.

Do not treat alternatives as additional physical behavior, or discard a split
piece merely because another alternative was selected elsewhere. Nested split
and alternative provenance is required to determine which coverage is replaced.
Equal current labels alone do not establish this. The first implementation may
restrict recombination to a single compatible alternative context.

## Committing recombination and backtracking

Compute and validate a proposed result before changing the tree. Then use one
multi-input derivation to create the result under a selected primary source,
absorb the other source nodes, and update the active frontier atomically.
Absorbed nodes remain inspectable but are no longer scheduled. Other unaffected
pieces of the family remain active. Do not discard already computed reach
segments merely because future checkpoint work has been combined.

Backtracking and insertion keep existing IDs stable. They also invalidate every
result that depends on changed work, including recombinations whose primary
parent lies in another subtree. Recompute those results as fresh nodes unless
their validity is established independently. Scheduler entries must recheck node
state and revision before publishing work after an absorption or invalidation.

## Guard-width diagnostic and later policy

The relevant comparison is crossing-time spread versus integration step, rather
than spatial width versus a time step. For a transverse, monotone crossing with
guard function `g` and a positive bound `m` on `abs(L_f g)` over the relevant
flow, `width(g)/m` is a useful crossing-time-spread diagnostic. Treat it as a
bound only when the required flow-domain and monotonicity hypotheses have been
validated. Tangencies and zero-containing derivative bounds need another rule.

Optional metadata may record the guard, crossing-time range, time progress,
enclosure size and recombination quality/cost estimates. Thresholds, candidate
indexing, synchronization, hull construction and scheduling policy belong to
Evolution and are not part of generic Tree.

## Verification and next decisions

Domain tests should cover identity inheritance for every operation; recursive
splits; repeated urgent-guard slices with distinct trajectory IDs and a shared
family; preservation of the no-jump continuation; alternative selection without
coverage loss; recombination across subtrees; time-state inclusion; incompatible
locations/clocks/budgets; and backtracking through an absorbed dependency.

Settle family equivalence, unbundling policy, the output-label convention and
numerical-alternative compatibility before implementation. The generic tree can
then be implemented independently of a particular enclosure or merge operator;
the hybrid wrapper enforces this document's domain invariants.
