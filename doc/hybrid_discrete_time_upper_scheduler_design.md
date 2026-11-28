# Deterministic hybrid discrete-time UPPER scheduler and orbit projection

Checkpoint: Evolution `multifunction`, design v6.44 (24 November 2026).

## Certified integer clock invariant

`HybridDiscreteTimeUpperScheduler::evolve` requires the certified
`HybridEnclosure::time_range()` of the initial set to be **contained in**
the singleton represented by `initial_index`. An intersection is not enough:
a non-singleton or uncertain clock cannot be silently interpreted as a
known integer sample. The same check is repeated for each dequeued node and
each one-tick successor. Failed checks raise `std::invalid_argument` even
when the requested horizon equals the initial index. The exact integer is
converted to `Real`, without truncating it to a machine integer or floating
point number.

## Orbit projection

`project_orbit` produces `HybridDiscreteTimeUpperOrbitProjection` with:
- `orbit.initial()`: the original enclosure;
- `orbit.reach()`: each reached non-root successor, including blocked and
  horizon leaves, without joining independent alternatives;
- `orbit.intermediate()`: each reached non-root enclosure strictly before
  the horizon, including blocked leaves;
- `orbit.final()`: only terminals marked `HORIZON`;
- `blocked`: only terminals marked `BLOCKED` (also when the root is blocked).

The node tree remains the authoritative source for parent indices, incoming
events and certified precondition evidence. `Orbit` is only a set-oriented
projection and does not encode provenance. A blocked leaf does not become a
horizon result; neither `final` nor `blocked` grants LOWER reachability.

No priority, zero-time jump chain, implicit guard restriction, set fusion or
sampling of difference inclusions is introduced. Unsupported inclusions still
raise rather than return a partial UPPER family. The `maximum_nodes` budget
remains a hard failure. The public `HybridDiscreteTimeEvolver<ValidatedTag>`
interface, full LOWER semantics and numerical refinements are deferred.

## Regression intent

The scheduler test checks mismatch at a zero horizon, a correct nonzero-time
restart, a mismatched shifted enclosure, an uncertain non-singleton clock,
an ordinary linear orbit, branching at the horizon, and a blocked root.
All previous tests of node budgets, event provenance, jumps and inclusion
rejection remain in place. The repository CI is the functional acceptance
gate; Format is intentionally not a gate until the redesign is complete.
