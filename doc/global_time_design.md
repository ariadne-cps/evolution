# External global time

Global time is an argument of the model, not a physical state variable.
`TimeVariable()` denotes the canonical global clock. A model can depend on
time without declaring `dot(t)=1`; guards, invariants, resets and auxiliary
outputs can use the same clock even when the differential field is autonomous.

## System declarations and hybrid argument spaces

```cpp
TimeVariable t;
RealVariable x("x"), u("u");
OrdinaryDifferentialAlgebraicEquation equation({dot(x) = t - x});
OrdinaryDifferentialAlgebraicInclusion inclusion(
    {dot(x) = t + u}, RealVariablesBox({-1 <= u <= 1}));
```

Both declarations have one physical state. The ordinary structured constructors
recognize the canonical clock when it occurs in an expression and has not
already been declared as a state, parameter, input or auxiliary variable. An
explicitly supplied time variable remains supported. Explicit variable roles
take precedence: a declared differential state named `t` keeps its physical
state role and is not silently converted into the external clock.

`HybridAutomatonInterface::time_variable(location)` declares the external clock.
For atomic automata the clock is also recognized in guards, invariants, resets,
auxiliary assignments and legacy differential assignments. Thus a mode's field
can be autonomous while its events are not. The physical state space stays
unchanged. `continuous_argument_space(location)` orders arguments as
`[state, optional time, constant model parameters]`. Free measurable inputs are
not guard/reset arguments. Equation fields and event functions use this space;
an inclusion still has no single-valued `dynamic_function`.

Resets may read global time but cannot reset it: their outputs are exactly the
target physical states. `check_mode` rejects a reset assigning to the external
clock. A separate model clock, such as a resettably elapsed state, remains a
physical state; enclosure dwell time is also distinct from global time.

## Numerical evolution and uncertain arrival

The existing enclosure already carries a global clock map `tau(s)`. A timed
step integrates the field on `tau(s)+h`, with local elapsed time `h`. The
private `AnchoredIntegrationStep` pairs the solver result with its argument
map, including the original source clock. It conditions the numerical step
on the clock origin and composes that argument back through `tau(s)`.
No synthetic system declaration, new physical state, persistent parameter or
time derivative is introduced. The conditioning origin is fixed only while
computing one numerical step; it is the updated source clock on the next step.
Model parameters retain their separate trajectory-constant semantics.

Ordinary and hybrid evolvers share this conditioning, including ordinary
evolution from a nonzero or uncertain initial clock. The ordinary workload
counts the requested elapsed duration; the hybrid horizon remains an absolute
clock bound. This construction supports parameter-dependent arrival times after resets and
retains state/time correlation for structured methods. Numerical bounds cover
the full argument maps, including points excluded by source constraints, so
composition never extrapolates beyond a certified domain. Explicitly selected
box integration keeps its existing outer Cartesian approximation contract.

Timed hybrid equation modes and all inclusion modes use the private
`hybrid_outer_flow` component for UPPER guard filtering and publication.
Equation and inclusion producers keep their own integration methods. Guards
and invariants compose with `[x(s),tau(s)]`; resets compose with the same
pre-reset map. Nonautonomous IMPACT candidates retain the closed guard
boundary, including zero-derivative contacts, without a directional exclusion.
No first-entry or whole-source coverage certificate is inferred.

The original autonomous ODAE crossing and timing strategies remain available.
A future certified nonautonomous crossing analysis must include the time term
`partial_t(g) + grad_x(g).f`; its ODAI counterpart needs suitable numerical
evidence. That analysis is not fabricated by the outer UPPER path.

## Auxiliary outputs and approximate evolution

`EnclosureFunctionArguments` explicitly distinguishes `STATE` and `STATE_TIME`.
Auxiliary maps, resets and constraints consume the declared layout, rather
than guessing it from an argument count. Time-aware auxiliary maps survive
subdivision, reconditioning and compatible recombination; they are evaluated
with the current clock after a reset. A state-only augmented paving cannot
represent nonempty time-dependent auxiliary outputs and rejects that operation;
joint state/time enclosure drawing remains available.

The continuous and hybrid equation simulators use RK4's native timed overload.
Every stage evaluates the appropriate external time; the method never appends
a differential clock to the state. Guards, invariants, resets and auxiliary
outputs use the sample's global clock. Single and multiple trajectories share
that convention. The final numerical step is clipped to the requested horizon;
explicit termination handles a non-dyadic horizon without a zero-step loop.

These trajectories remain approximate. Hybrid simulation retains its existing
strict-positive sampled guard rule and can miss within-step contacts. This
does not become validated event localization merely because time is explicit.

## Capability boundaries and validation

| Consumer | Autonomous | Explicit global time | Remaining limits |
| --- | --- | --- | --- |
| Ordinary equation evolver | yes | yes, including time- and fixed-parameter-dependent auxiliary outputs | Timed/parameterised preconditioning remains unsupported |
| Ordinary inclusion evolver | UPPER | UPPER with supported numerical field | Existing representation and solver capability limits |
| Hybrid equation evolver | Analytic path without parameters; UPPER outer-flow with parameters | UPPER outer-flow, including time/parameter events and outputs | Timed/parameterized LOWER and preconditioning rejected |
| Hybrid inclusion evolver | UPPER, including fixed parameters | UPPER, with structured or box method and fixed parameters | LOWER and existing method/domain capability limits |
| Ordinary/hybrid approximate evolvers | ODAE and structured ODAI, with fixed parameters | yes, including stage times, events and equation outputs | Explicit input policy required for ODAI; interval-Taylor sampling and general parameter sampling remain unsupported |

Typed hybrid composition retains its separate existing capability restriction.
The new clock support does not bypass it. The implemented
[paradigm API](evolver_paradigms_design.md) selects approximate trajectories or
validated enclosures explicitly. [Fixed-parameter auxiliary layouts](constant_parameter_auxiliaries_design.md)
also preserve physical parameter intervals through ordinary validated evolution.
[Hybrid parameter bindings](hybrid_constant_parameters_design.md) retain those
values across validated mode transitions, including state/time/parameter reset
arguments and target equation outputs.
Full Python/C++ tutorial parity remains the final Evolution activity.

Analytic regressions verify `x'=t`, source-clock correlations, auxiliary
identities through structural operations, `x'=t+u` with an uncertain arrival,
timed invariants/guards/resets on both autonomous and nonautonomous fields,
structured/box methods and flat/tree backends. Separate RK4 tests integrate
`t^2` from a nonzero initial time and check stage-time evaluation. Simulator
tests cover single/multiple trajectories, timed outputs, reset clock continuity
and a non-dyadic final horizon. No additional differential clock is present in
the tested physical state spaces.


Validation: macOS Debug with warnings as errors passes 46/46 integrated CTest
targets (32 Evolution including Python, 13 System C++ and the native Integration
RK4 suite). Formatting passes for 92 Evolution, 39 System and 42 Integration
C++ files. Installed C++ and existing Python tutorials are checked separately.
CI results belong to each published commit and are not inferred from local tests.
