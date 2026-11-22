# System families: simulation and validated evolution

This is the implementation roadmap adopted in design v6.20 on 22 November
2026. Planned names in this document do not declare implemented APIs.

| Native family | Validated evolution | Numerical simulation | Current status |
| --- | --- | --- | --- |
| Discrete time (`IteratedMap`) | `DiscreteTimeEvolver` | `DiscreteTimeSimulator` | Evolver exists; simulator planned |
| Continuous time (ODAE / ODAI) | `ContinuousTimeEvolver` | `ContinuousTimeSimulator` | Evolver accepts both with explicit time; simulator accepts autonomous/nonautonomous ODAE |
| Hybrid discrete time | `HybridDiscreteTimeEvolver` | `HybridDiscreteTimeSimulator` | Native System model and both algorithms planned |
| Hybrid continuous time | `HybridContinuousTimeEvolver` | `HybridContinuousTimeSimulator` | ODAE/ODAI timed UPPER evolution; equation simulator supports external time; other extensions planned |

`IteratedMap` keeps its existing name and contract. Hybrid discrete time needs
its own native model in `ariadne-cps/system` before either algorithm is added.
The tick semantics must define guard evaluation before/after the update,
priorities and simultaneous events, resets, multiple jumps at the same tick,
termination and the domain of integer indices. Continuous event semantics must
not be silently reused on a time grid.

## Shared responsibilities and separate numerical kernels

Share system ownership, initial-coordinate handling, labels and auxiliary
outputs, configuration, termination, trajectory collection, progress and
diagnostics, branch provenance and result production/publication where the
contracts agree. Reuse flat/tree scheduling and concurrency infrastructure for
both hybrid families when appropriate; do not duplicate entire hybrid engines.
Subdivision and reconditioning require a valid numerical representation.

Keep the step kernels explicit: discrete iteration applies a map/relation at
an integer index, whereas continuous integration produces a real-time reach
tube and endpoint. First-entry maps, Lie derivatives, within-step maxima and
creep belong to continuous localization, not discrete ticks. Generalization
must expose these capabilities rather than implying equivalent certificates.

Validated evolvers produce sets with their supported approximation semantics.
Simulators produce approximate sampled trajectories. Sampling a trajectory
does not provide a validated enclosure, and API similarity cannot grant LOWER,
implicit DAE/DAI, multipatch or preconditioning support.

## Explicit continuous-to-discrete conversion

A discrete algorithm accepts a continuous system only through a converter
explicitly selected by its caller. No implicit ODAE/ODAI constructor or inverse
discrete-to-continuous conversion is planned.

A converter declares the step duration, time origin/grid, integration or
discretization scheme, parameter treatment and input policy. Approximate
simulation and validated step propagation have separate contracts. An ODAI
conversion must retain nondeterminism through a suitable typed discrete
relation/multimap, or reject the operation if `IteratedMap` cannot represent it.

Hybrid conversion additionally declares how within-step events, precedence,
resets, zero-time behavior and tangencies are handled. Endpoint guard samples
alone may miss crossings or contacts and cannot justify validated conservative
propagation. An explicitly approximate simulation conversion may use a weaker
contract, with its limitations visible to the caller.

## Implementation order and acceptance

1. Consolidate current hybrid continuous evolution, including ODAI UPPER
   IMPACT contacts with zero derivative, competing blockers, the closed horizon
   and parameter-dependent arrival. Keep unsupported LOWER cases explicit.
2. Complete continuous/hybrid time, fixed-parameter, input and auxiliary
   capabilities only for representations the numerical methods support.
3. Define and test the native hybrid discrete-time system in System.
4. Extract compatible hybrid infrastructure and implement the discrete hybrid
   evolver with declared UPPER/LOWER capabilities.
5. Add `DiscreteTimeSimulator`; extend `ContinuousTimeSimulator` with ODAI,
   explicit time, trajectory-constant parameters and reproducible input policies.
6. Implement `HybridDiscreteTimeSimulator` using the same native tick model;
   extend `HybridContinuousTimeSimulator` with the corresponding continuous
   capabilities and contact handling.
7. Implement the optional converters with separate simulation/validated
   contracts, testing within-step events and rejection of lost guarantees.
8. Verify the four-family/two-operation matrix: compatible and incompatible
   systems, integer/real horizons, reproducible seeds, inputs and parameters,
   simultaneous events, branch/reset behavior and numerical guarantees.
9. As the final Evolution activity, add missing Python bindings and give the
   Python tutorial functional parity with C++ for every available capability.
   Update installation, exports and documentation without compatibility aliases.

IMPACT and external global time advance the first two items; constant model
parameters in the hybrid path and further simulator capabilities remain open. All planned APIs, converters
and the native hybrid discrete-time model remain unimplemented. The detailed
[continuous hybrid handler design](hybrid_transition_handlers_design.md)
records the implemented event capabilities and their mathematical boundaries.


## Proposed paradigm axis (not yet implemented)

The proposed API migration makes the numerical paradigm part of each evolver
type, using `ApproximateTag` and `ValidatedTag`, for example
`ContinuousTimeEvolver<ApproximateTag>` and
`ContinuousTimeEvolver<ValidatedTag>`, and corresponding hybrid/discrete types.
The paradigm must select actual input/result types and available operations,
not just change a name. Approximate trajectories and validated enclosures
remain distinct; UPPER/LOWER is a separate semantic axis with explicit numerical
capabilities. Shared interfaces and compatible infrastructure should not force
the same step algorithm or certificates. The present canonical classes have
not been renamed in the global-time implementation; this proposal needs an
atomic API/consumer migration, with Python tutorial parity still the final task.
