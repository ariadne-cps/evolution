# Ordinary discrete-time equations: witness-based LOWER audit

Checkpoint v6.53 (24 November 2026); repository `ariadne-cps/evolution`,
branch `multifunction`.

## Distinguish the four notions used in Ariadne

- **LOWER:** every published bounded enclosure `E` has at least one
  actually reachable witness, i.e. `E ∩ R != ∅`. The enclosure is
  **not** required to lie inside `R`.
- **INNER:** every represented point is actually reachable:
  `I ⊆ R`. This is strictly stronger than the project-specific LOWER
  witness contract and is not claimed by the validated evolver.
- **UPPER:** the numerical enclosure conservatively includes all
  evolution possibilities for the source and branch it represents.
  Correct images of selected branches alone need not be exhaustive.
- **OUTER:** all reachable points are covered by the union of the
  computed UPPER enclosures: `R ⊆ ⋃ U_i`. This additionally requires
  *complete exploration* of all applicable evolution alternatives.

For any bounded LOWER enclosure `E` with reachable witness `r`,
every point `x∈E` satisfies
`distance(x,R) <= distance(x,r) <= diameter(E)`.
There is no need to compute a distance in code; the diameter supplies
the bound automatically. If `E` has no reachable witness, the
diameter argument does not establish LOWER correctness.

## Why the deterministic discrete map needs no separate INNER kernel

Assume a **known nonempty initial enclosure** `E_0` with a real
initial witness `r_0`, and a deterministic map `f` defined along
its evolution. Suppose `apply_discrete_time_map_step` always returns
a sound UPPER image `E_(n+1)` of `f(E_n)`. Inductively,
`r_n∈E_n` implies `r_(n+1)=f(r_n)∈f(E_n)⊆E_(n+1)`.
Thus a single unsplit enclosure retains a reachable witness at
each tick, even when it is geometrically much larger than the true
reachable set. This proves the LOWER witness obligation on this
supported subset; it does **not** prove INNER.

The current `DiscreteTimeEvolver<ValidatedTag>::_evolution`
subdivides only under `Semantics::UPPER`, where preserving the
complete cover is appropriate. In LOWER, it propagates the unsplit
enclosure or, when configured for premature termination and the
radius budget is exceeded, publishes **no successor** on that path.
That distinction is mathematically necessary: splitting a
witness-bearing enclosure does not prove that *every* child contains
a reachable witness.

For an ordinary deterministic `DifferenceAlgebraicEquation`, there are no hybrid
guarded resets or alternative transitions. By propagating all map
images, a finite-horizon UPPER computation can additionally be OUTER
if it exhaustively covers its intended initial set and completes
the horizon. Its OUTER status must not be inferred from the word
UPPER alone in systems with branching, truncation or incomplete
exploration.

## Limits of this audit

The argument assumes an initial source with a **known witness**.
An arbitrary constrained initial enclosure can be empty, and the
current `Enclosure::is_empty()` may return INDETERMINATE for
feasibility. The public `orbit` overload does not demand an explicit
witness/non-emptiness certificate from every constrained caller.
Do not mark that general case as proven LOWER without a separate
input precondition or constructive feasibility certificate.

We are auditing only ordinary discrete equations. For hybrid LOWER,
each emitted branch needs a witness. A guard certified fully active
preserves the existing witness through the reset; a merely possible
guard does not. LOWER cannot silently split a source or retain a
branch solely because a UPPER guard test was unresolved.

## Regression coverage and acceptance

The new `test_discrete_time_evolver::test_lower_witness_and_upper_subdivision`
checks an analytic witness for the affine map `x[n+1]=x[n]+1`
from `x0=0` through two ticks, verifies LOWER never subdivides
even if UPPER subdivisions are enabled, and checks a nonlinear
`x[n+1]=x[n]^2` witness at `x0=0`. Existing tests cover
zero-horizon behavior and LOWER premature termination. Witness
coverage is grounded in these explicit analytic trajectories rather
than a fabricated inner enclosure or a measured distance.

These tests are submitted but not considered passed until the
functional GitHub CI completes for this exact commit; in addition
to the ordinary discrete suite, existing continuous and discrete
hybrid tests must stay green. Format remains outside the gate.

## Capacity matrix delta

Ordinary discrete equations / Validated LOWER:
**audit established for known nonempty initial enclosures**, provided
the deterministic numerical image is sound and no additional
branching/partitioning is introduced. The unconstrained full API
contract for arbitrary potentially empty initial enclosures remains
an explicit open assumption. This audit introduces no new
Approximate evolver or inclusion semantics.

## v6.54 — Compile fix, no semantic changes

The v6.53 Debug CI failed to compile `test_discrete_time_evolver.cpp`
because the single-element brace arguments to the `DifferenceAlgebraicEquation`
constructor were ambiguous between the copy constructor and
`DifferenceAlgebraicEquation(const List<PrimedRealAssignment>&)`. Both the affine and
quadratic regression models now explicitly construct
`List<PrimedRealAssignment>`. The mathematical witness argument, the
LOWER non-splitting requirement, and the test expectations are unchanged.
The v6.53 failure was a test compilation bug, not a failure of LOWER
semantics. v6.54 remains pending CI on its own commit.

## v6.57 — Constructive witness gate for constrained LOWER sources

The ordinary discrete validated public `orbit` entry now checks a
**sufficient** certificate before publishing any LOWER enclosure.
For an initial `Enclosure` with no parameter constraints, a nonempty
parameter domain supplies the witness. For constrained sources,
`has_certified_initial_witness` checks the midpoint and two
parameter-domain vertices (all-lower, all-upper) using validated
interval evaluation of every parameter constraint model. A candidate
is certified only if its whole validated function range is contained
in the permitted constraint interval. Passing therefore constructs a
feasible parameter value; it does not merely fail to prove emptiness.

If no tested candidate certifies all constraints, LOWER rejects the
request with `std::invalid_argument` **even at zero horizon**.
This is a conservative failure of certification, not a claim that
the underlying mathematical initial set is empty. UPPER semantics
and its exploration are unchanged.

For an accepted source, the deterministic `DifferenceAlgebraicEquation`
map transports a known witness at each tick, since numerical images
enclose the map image and LOWER never subdivides. Its bounded enclosure
diameter automatically bounds the distance to a reachable witness.
The entry check does not compute a distance and does not claim INNER.

The new discrete regression tests cover a strictly feasible
constraint at the midpoint, a feasible constraint away from the
midpoint, a provably impossible constrained source rejected at
both zero and positive horizons, and the original analytic witness
propagation. Existing ordinary continuous, hybrid continuous and
hybrid discrete test suites remain mandatory in CI.

**Matrix:** Ordinary discrete equation / Validated LOWER remains
partial in the sense that feasible arbitrary constrained initial
sets may be rejected when none of the finitely tested candidate
points provides a proof. This checkpoint strengthens the **soundness
contract for accepted inputs**, not completeness of admissible inputs.
No extra numerical or paradigm capability is claimed. The exact new
commit's functional CI must pass before this version is accepted.
