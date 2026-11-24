# Evolution contracts by numerical paradigm

The numerical paradigm is an explicit template argument, without a default:

| Family | Approximate trajectories | Validated enclosures |
| --- | --- | --- |
| Continuous | `ContinuousTimeEvolver<ApproximateTag>` | `ContinuousTimeEvolver<ValidatedTag>` |
| Hybrid continuous | `HybridContinuousTimeEvolver<ApproximateTag>` | `HybridContinuousTimeEvolver<ValidatedTag>` |
| Discrete | Planned; incomplete primary template | `DiscreteTimeEvolver<ValidatedTag>` |
| Hybrid discrete | Native system and both algorithms planned | Native system and both algorithms planned |

The canonical family headers include the implemented specializations. Each
specialization owns its numerical implementation, in a corresponding
`*_approximate_evolver` or `*_validated_evolver` header/source. The former
simulator names and headers are removed, without compatibility aliases or
wrappers. Discrete iteration keeps its own numerical kernel.

## Distinct interfaces and results

`EvolverInterface<P, System, Initial, Termination>` has two explicit partial
specializations. Each exposes `Paradigm`, `SystemType`, `InitialType`,
`TerminationType`, `OrbitType`, and `ResultType`, plus polymorphic `system()`
and owning `clone()` access.

- The approximate interface accepts `orbit(initial, termination)`. Continuous
  initial data are lists of labelled approximate points; hybrid initial data
  are hybrid approximate points. Their orbits contain sampled curves. The
  native box/set and multiple-trajectory overloads remain available.
- The validated interface accepts `orbit(initial, termination, semantics)`.
  Initial data and results use labelled or hybrid enclosures. `EnclosureType`
  and `EnclosureListType` belong to this interface only. Existing box/set
  overloads, integration selection, flat/tree backends and certificates remain.

UPPER/LOWER is independent of the paradigm. An approximate trajectory cannot
request LOWER or UPPER, and a validated call must supply its set semantics.
The numerical capability gates still apply: for example ODAI and timed hybrid
flows currently support UPPER only. This migration adds no inner certificates.

Configurations and hybrid interfaces, bases and factories also carry the
paradigm. `EvolutionTraits<P, System>` associates only implemented combinations:
ODAE and ODAI have both continuous paradigms, with explicit sampling policies
required for approximate ODAI, and
`DifferenceAlgebraicEquation` has validated discrete evolution. Hybrid model traits select the
corresponding hybrid interface; its runtime model capability checks remain.
Structured approximate ODAI use an explicit input policy; interval-Taylor
multifunctions remain rejected. No approximate discrete algorithm is hidden
behind a fallback. Unsupported template combinations remain incomplete.
Tag inheritance does not enable dispatch to another paradigm specialization.

## Python and migration boundary

The existing bindings use the established Python template registry:
`ContinuousTimeEvolver[ApproximateTag]`,
`ContinuousTimeEvolver[ValidatedTag]`,
`DiscreteTimeEvolver[ValidatedTag]`, and
`HybridContinuousTimeEvolver[ApproximateTag]`.
The registered concrete classes have explicit names such as
`ContinuousTimeApproximateEvolver` and `ContinuousTimeValidatedEvolver`.
Existing configuration/result bindings are migrated with their owning types.
The hybrid validated bindings and full Python/C++ tutorial parity remain the
final Evolution activity, as previously agreed. An absent Python registry
entry does not claim an absent C++ implementation.

## Verification and remaining work

Compile-time contracts check paradigm identity, distinct initial/result types,
the two orbit signatures, unsupported systems and incomplete combinations.
Polymorphic runtime tests exercise continuous approximate/validated clocks and
hybrid approximate flow through owning clones. Existing numerical suites cover
the migrated specializations, including global time, inclusions and events.
Python tests check specialization selection, unsupported registry keys and
rejection of invalid orbit signatures. Installation and the existing C++/Python
tutorials must still work with the new names.

Next, complete missing numerical capabilities and the native hybrid discrete
model in the [system-family roadmap](system_families_design.md). Factor common
infrastructure only where input/result and certificate contracts agree. Keep
sampling policies, further parameter capabilities, contact handling, converters and final
Python tutorial parity as explicit unfinished activities.

Explicit fixed parameter values are now supported by ordinary and hybrid ODAE
approximate constructors, retaining the same paradigm interface and physical
state dimension. See [approximate parameter selection](approximate_constant_parameters_design.md).
The corresponding Python constructor bindings remain deferred to the final
binding/tutorial phase.

ContinuousDynamics is now shared metadata without a paradigm tag, extracted
into its own public header/source. Both continuous SystemType aliases use this
owning ODAE/ODAI view; numerical kernels and policies stay in the evolvers. The
original declarations are retained. Existing Python ODAE construction is ported
to the explicit equation overload; new sampling bindings remain deferred.
