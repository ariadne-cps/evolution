# Trajectory-constant parameters in validated hybrid modes

Parameterized ODAE and ODAI modes use UPPER outer-flow evolution. This retains
their physical parameter binding and evaluates fields, guards, invariants,
resets and equation auxiliary outputs with the declared `[state, optional time,
parameters]` layout. Autonomous equations with parameters use the same clock
conditioning as ordinary equations; the clock is not a physical state.
No new first-entry, impact-direction or LOWER certificate is claimed.

The first parameterized mode binds its declared tuple once. Subsequent modes
that declare parameters must use the same ordered names and a physical domain
containing the bound domain. A parameter-free mode keeps the global binding
while its own functions consume only state and optional time. Partial/different
parameter tuples and domain refinement inferred from guards remain unsupported;
declare the full tuple, including unused parameters, when a mode uses it.

Target binding checks/initialization occur before applying a reset and its
target auxiliary map. Resets change physical state, not the constant parameter
map or global clock. Parameter assignments in resets retain System's existing
rejection. Transition metadata carries the target parameter declaration;
all jump publication paths use one preparation/reset operation.

The parameter argument count comes from an explicit named binding, never from
an arity guess. Guard constraints compose the complete argument map onto the
existing enclosure domain. Hybrid grouping must distinguish incompatible
physical tuples before requesting recombination. Compatible sources use the
physical-map transport already implemented for ordinary continuation.

Existing inclusion solver selection, certified provider-domain checks and
preconditioning limits remain. Parameterized modes require UPPER and reject
preconditioning. Nonparameterized autonomous equation modes retain their
analytic crossing path, including when carrying an unused physical tuple.

Seven dedicated regression cases verify:

- Autonomous/timed ODAE flow and repeated continuation on flat/tree backends,
  using the analytic solution of x'=p or x'=p+t and equation auxiliary identities.
- Autonomous/timed ODAI outer flow with structured and box methods on both
  backends, retaining the analytic extrema for p in [2,3], u in [-1/8,1/8].
- URGENT/PERMISSIVE parameter guards, invariants and immediate resets.
- Time/parameter URGENT and IMPACT guards, correlated resets and target
  equation outputs, zero dwell time and IMPACT equality on the boundary.
- First parameter binding at a target and preservation through an unparameterized
  mode that still uses the analytic crossing path.
- Incompatible names/domains, assignment to a parameter, LOWER and
  preconditioning rejection, and explicit inclusion method selection.
- Grouping of compatible checkpoints and exclusion of different physical tuples.

ODAI auxiliary declarations are not provided by the current System constructor;
this change does not invent them. Existing equation auxiliary declarations use
all four explicit layouts.

Local verification passes all 50 CTest suites (36 Evolution including Python,
13 native System and one native Integration), the complete test/binding/benchmark
build, a fresh install and both existing standalone tutorials. Before commit,
apply-format and check-format pass on 101 tracked C++ files.

Paving is owned by the separate `decoupling` branch. Python tutorial parity and
missing binding expansion remain the final Evolution activity.
