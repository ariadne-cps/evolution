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
