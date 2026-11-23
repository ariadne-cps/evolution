# Cross-repository hybrid automaton interface migration (v6.29)

System provides:

- `HybridAutomatonInterface`: time-independent structural queries only.
- `ContinuousTimeHybridAutomatonInterface`: continuous time, ODAE/ODAI,
  guards/invariants/resets and continuous event semantics.
- `ContinuousTimeHybridAutomaton` and
  `CompositeContinuousTimeHybridAutomaton`: owned concrete models.

Evolution consumes the continuous specialization in
`HybridContinuousTimeEvolverInterface<ApproximateTag>` and
`HybridContinuousTimeEvolverInterface<ValidatedTag>`, the corresponding
factory, concrete evolution kernels and trait mappings. The module uses the
new System gitlink, and existing bindings/tests/tutorials are renamed.
The numeric algorithms are intentionally unchanged. No compatibility aliases
are introduced.

## Numerical boundaries

No discrete-time automaton has been implemented. It must use a separate
integer-time transition relation, not `HybridTime` or continuous first-contact
certificates. Reuse the structural API and generic scheduling and result
publication only where guarantees are identical. ODAI UPPER/LOWER capability
boundaries and sampled approximate trajectories remain unchanged.

## Validation

The coordinated commit requires System and Evolution CI verification across
platforms and a consuming build. This document does not claim the migration is
fully tested until the CI closes. The paving refactor remains in the separate
`decoupling` branch.
