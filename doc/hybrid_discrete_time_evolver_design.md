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

## v6.48: validated event-free bounded-input inclusion UPPER progression

For a DifferenceAlgebraicInclusion in an event-free mode without auxiliary
outputs, each tick appends a fresh INPUT parameter axis for each bounded
input and composes the validated difference update with the existing
parameterized state and the new independent input coordinates. No sample
is selected; previous input axes remain independent and the original
enclosure stays unchanged. This produces an UPPER set-valued image.

Transitions involving inclusion sources or targets and input-dependent
auxiliary outputs remain unsupported, requiring the missing quantified
guard and auxiliary representation. LOWER remains explicitly rejected.

Tests cover image construction, independent inputs across two ticks,
the full x=0, u in [-1,1] two-tick extent [-2,2], scheduler and Orbit
integration, the public validated evolver, immutability and rejection
of unsupported auxiliary/transition paths. Continuous guard and evolution
regressions must continue to pass alongside all discrete suites.
Acceptance depends on CI for the actual new commit, not on prior commits.

## v6.49: scheduler capability gate repaired

The v6.48 image/executor subset had a mismatched scheduler precondition:
`HybridDiscreteTimeUpperScheduler::evolve` still rejected all inclusion
modes before reaching the new event-free kernel. The initial capability
gate now accepts equation modes or only inclusion modes with no outgoing
events and no auxiliary outputs. Unsupported inclusion models fail even
at a zero horizon instead of appearing to succeed without a step.
The existing state-space, certified integer clock and hard node-budget
checks remain in force. The scheduler and public evolver two-tick
inclusion tests now exercise a reachable path through this gate.

This was a functional integration bug in v6.48, not a reason to weaken
the tests. Functional CI remains required on the repair commit, and
both continuous and discrete hybrid suites remain active. The
quantified-guard extension is deferred until that CI is verified.

## v6.50 — State-only guards and resets from bounded-input inclusion modes

The native discrete System model forbids free input symbols in guards,
progression invariants and transition resets. Without auxiliary outputs,
those predicates depend solely on the **pre-tick state**. The validated
shared source-guard classifier now handles such inclusion modes with the
same closed-side restriction used for deterministic source modes.
An unknown guard produces independent UPPER alternatives; URGENT uses
the pre-state complement for progression, while PERMISSIVE applies the
explicit progression invariant. Multiple events are not prioritized.

A bounded-input inclusion's progression adds independent INPUT coordinates
for the current tick as before. A transition from an inclusion instead
executes the state-only reset to a deterministic equation mode, consumes
one integer tick and **does not** select or append a current-tick input.
Previous input coordinates remain correlated with the source state.
All target inclusion modes and auxiliary-bearing inclusion sources
remain rejected, including at a zero horizon; this is enforced by a
single shared capability gate in predicate conversion, image application,
executor and scheduler. Unsupported inclusion dynamics never return a
partial UPPER orbit. LOWER remains unsupported.

**Two-family regression obligation:** discrete predicate tests verify
TRUE/FALSE/UNKNOWN, urgent complement, both constrained alternatives,
immutable source, input axes only on progression, exact tick and event
history. Step-image tests verify state-only inclusion reset; scheduler
and public-evolver tests verify UPPER reach/final and provenance.
The existing continuous hybrid guard classifier, continuous integration
and event/impact suites remain enabled and must be green together with
the discrete suites in Ubuntu, macOS, Windows, Coverage and Sanitizers.
Format remains non-blocking. No CI result can be transferred from
v6.49 to the v6.50 commit.

Next: verify functional CI; correct regressions without weakening
semantic expectations. Only then expand inclusion auxiliary or target
capabilities and eventually introduce approximate discrete evolvers.
Python binding and tutorial parity remain the final phase.

## v6.52 — Regression expectation corrected after v6.50 inclusion support

Debug CI for v6.50 (`0cdb905`) reports a failing case
`test_hybrid_discrete_time_predicate_evaluator: reject_inclusion_and_wrong_space`.
The test expected `std::invalid_argument` when classifying `x >= 0`
in an inclusion `x[n+1] = x[n] + u[n]`, with independent bounded input,
no auxiliary outputs, and pre-state `x in [0,1]`. That source capability
was intentionally implemented by v6.50; no exception is therefore correct.
The failure is an **outdated test expectation**, not evidence that the
supported classifier should be disabled.

The revised test requires three validated predicate outcomes (TRUE for
`x >= -1`, FALSE for `x >= 2`, and UNKNOWN at the `x >= 0`
boundary). It **retains negative coverage** with a separate inclusion
containing an input-dependent auxiliary assignment and with a mismatched
equation state space; both must throw. No numerical implementation is
modified and no accepted model is silently rejected. The rest of the
continuous and discrete regression suites remain in the CI gate.

The v6.50 result remains marked *failed*; the correction commit's CI
must complete on Ubuntu, macOS, Windows, Coverage and Sanitizers before
promoting this checkpoint to validated. Format remains non-blocking.
