# Explicit constant parameters in approximate ODAE evolution

Both ordinary and hybrid continuous approximate evolvers accept a constructor
with a RealValuation of named, explicitly selected model parameter values.
The original constructor is the empty selection. A parameterized mode requires
all its declared values; no implicit midpoint or random choice is introduced.

Values must be finite and their binary64 interval enclosure must be provably
contained in the declared domain. Values are copied into the evolver, checked against declared physical domains
and substituted into fields and equation auxiliary functions before numerical
evaluation. Hybrid guards, invariants and source resets use the same values.
Parameter-free modes retain the selection; later modes select their parameters
by name, allowing reordered or partial tuples when each selected value fits
the declaration. Validation occurs when entering a hybrid mode, before recording
its initial output or applying a transition into it. Ordinary selections require
exactly the declared names; hybrid selections can include later-mode names.

RK4 retains the physical state dimension and evaluates external time at every
stage. Parameters are constants in the prepared function, not differential
states. A selected value stays fixed for every trajectory and every invocation
on the same evolver, including clones and batches. Distinct selections can be
simulated with distinct evolvers. General random/grid parameter policies and
ODAI measurable-input sampling remain separate activities.

The returned orbit still labels state and equation auxiliary variables. The
constant_parameters getter exposes the selected valuation for provenance.
This is approximate numerical evolution, with the existing sampled event
semantics; it grants no validated UPPER/LOWER or first-contact guarantee.

Six regression cases cover autonomous/timed flow, named ordering, multiple trajectories,
cloning, parameter-dependent equation outputs/guards/invariants/resets, passage
through parameter-free modes, and missing/out-of-domain selections. Paving
changes remain on decoupling. Python binding expansion/tutorial parity remain
the final Evolution activity.

For example, an equation declaring p in [2,3] can be simulated with:

```cpp
ContinuousTimeEvolver<ApproximateTag> evolver(
    equation, RealValuation({{p, Real(2.5_x)}}));
```

The hybrid constructor accepts the same valuation form. Empty selections
remain valid for unparameterized systems. Parameter selection does not change
the existing polymorphic orbit signature or system type.

Local validation passes all 54 CTest suites: 37 Evolution including Python,
13 native System and four native Integration suites. The complete test/binding/
benchmark build, fresh install and existing standalone C++/Python tutorials
succeed. apply-format and check-format pass on 104 tracked C++ files before
publication. No submodule or paving changes are included.
