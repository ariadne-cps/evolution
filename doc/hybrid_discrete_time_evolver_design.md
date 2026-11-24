# Public hybrid discrete-time validated UPPER evolver

Design checkpoint v6.45, 24 November 2026, Evolution branch `multifunction`.

## Contract
`HybridDiscreteTimeEvolver<ValidatedTag>` is a public specialization
of `EvolverInterface<ValidatedTag, DiscreteTimeHybridAutomatonInterface,
HybridEnclosure, HybridDiscreteTimeTerminationCriterion>`.
The termination criterion carries exact arbitrary-precision initial and
terminal indices, preventing lossy inference of a tick from a validated
real-time enclosure. It rejects inverted horizons; the scheduler
certifies that every source and successor enclosure has a singleton clock
coinciding with the declared index.

The evolver clones the polymorphic System model. Copy and `clone()`
deep-copy the model and the configuration. The positive
`configuration().maximum_nodes()` budget defaults to 10000; exceeding
it fails rather than yielding an incomplete UPPER enclosure family.

## Results
- `orbit(..., Semantics::UPPER)` returns `Orbit<HybridEnclosure>`;
  only leaves at the requested horizon are in `orbit.final()`.
- `orbit_with_blocked(...)` also returns premature `BLOCKED` leaves
  separately, even if the root is blocked.
- `explore(...)` returns the complete parent/event/index tree and
  both terminal kinds; Orbit does not encode per-branch provenance.
- `Semantics::LOWER` is rejected unconditionally, regardless of
  certified preconditions. No LOWER semantics are inferred from UPPER.
- No difference inclusions are sampled; unsupported images fail.
  Guard-specific subset restriction, subdivisions, reconditioning,
  priorities, zero-time jump chains and implicit merging remain absent.
- No approximate hybrid discrete specialization or new Python binding.
  Complete binding/tutorial parity remains the final phase.

## Validation
A dedicated C++ test checks header/traits/interface, orbit horizons,
permissive alternatives, event provenance, blocked roots, nonzero initial
ticks, clone/configuration independence, budget and invalid inputs, LOWER
rejection and inclusion rejection. The tests are registered in CMake.
Actual CI outcomes must be checked for this commit before declaring it
functionally validated. Format remains deliberately excluded until the
redesign is complete. System and paving's `decoupling` branch are not
modified.

Next: check CI; continue with sound UPPER image refinement and quantified
input handling before introducing inclusion modes. Implement the
approximate family separately, and leave Python/tutorial parity last.

## v6.46: pre-state certified branch refinement

For `execute_certified` only, each branch is imaged from an independent
copy of the source with validated parameter constraints encoding the
**closed** pre-state predicate region. A transition receives its guard,
while a progression receives the intersection of every event's
progression invariant (the complement of each urgent guard, or the
explicit permissive invariant). The indicators use System's
`make_constraint_function` with sorted algebraic substitution and are
composed with the original source state map, preserving dependencies.
A nonnegative indicator condition includes the boundary even when the
symbolic predicate is strict, so this operation is a sound outer
relaxation, not a certified equality with the exact admissible subset.

The older callback-based `execute` method still produces full-source
images: the arbitrary callback does not supply a numerical restriction
certificate. No LOWER claim, branch priority, zero-time transition,
implicit fusion or difference-inclusion sampling is added. The
regression verifies permissive split guards including a dependent
auxiliary, the urgent guard's complement, independent constraints for
each alternative, no mutation of the original source, and unchanged
callback-based semantics. CI on the new commit is pending until checked.

## v6.47: shared source-guard certificates across both hybrid time models

The private `hybrid_guard_predicate` component now contains validated
`classify_source_guard` and `restrict_source_guard`. The classification
checks emptiness of the complete source, then whole-source strict positivity,
and finally emptiness of a candidate restricted to the **closed** positive
guard side. It returns EXCLUDED, UNRESOLVED or FULL_SOURCE and never
turns a partial UPPER candidate into a LOWER reachability certificate.

The continuous hybrid `classify_initial_guard` delegates to this same
component, preserving its existing enum and UPPER/LOWER policies; source/tube
activation and immediate jump candidates use the common restriction.
Crossing geometry, integration, and impact processing remain continuous-only.

The discrete hybrid evaluator constructs an indicator with canonical algebraic
substitution and delegates classification to the same component. The certified
tick executor applies the same restriction on independent pre-state copies;
the external-evidence callback path remains full-source. Inclusion sampling
and LOWER evolution remain rejected.

Tests for **both** cases are mandatory: `test_hybrid_step_contracts` checks
continuous guard status, source immutability, equivalent restriction, and
UPPER/LOWER gating, alongside the existing continuous evolution suites;
`test_hybrid_discrete_time_predicate_evaluator` checks all three statuses
with correlated auxiliary outputs and retains urgent/permissive branch tests,
alongside the public discrete evolver and scheduler suites. A CI pass on
the actual new commit is required before this version is considered validated.
Format remains outside the functional acceptance gate.
