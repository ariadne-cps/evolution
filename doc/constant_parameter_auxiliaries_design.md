# Constant parameters in ordinary validated auxiliary outputs

Ordinary validated equation evolution already carries trajectory-constant
parameter coordinates with `EnclosureVariableKind::PARAMETER`. Add auxiliary
layouts `[state, parameters]` and `[state, global time, parameters]` with an
explicit parameter count. Coordinate selection uses the last declared
`PARAMETER` axes, in their existing order, when no explicit binding is present.
With a named binding, it evaluates the [physical parameter map](constant_parameter_maps_design.md)
instead of interpreting normalized coordinates as physical values. Numerical
integration and auxiliary evaluation share that selection.

Attach the system auxiliary mapping after adjoining the model parameter
domain, for box, constrained-set and enclosure initial data. The clock remains
external and the physical state dimension is unchanged. A missing parameter
coordinate or inconsistent layout/count is an error. Reconditioning protects
the dependencies of named physical parameter maps; splitting/restriction
transports those maps to the corresponding source subdomains.

This step supports UPPER ordinary ODAE evolution with compatible integrators.
Existing LOWER, preconditioned, hybrid-parameter and approximate-sampling
capability gates remain. Named bindings support continuation under declarations
with the same ordered parameter names and a containing physical domain.

Multi-source recombination transports compatible named physical parameter
maps jointly with state and time. Unbound PARAMETER coordinates remain rejected
because their physical identities cannot be inferred after normalization.
Single-source identity results remain valid. Joint state/time/auxiliary
enclosure evaluation and drawing remain available. Paving work is deferred to
the separate `decoupling` branch.

The implemented regressions compare correlated auxiliary maps against analytic
solutions for autonomous and timed equations, including non-unit parameter
domains, multiple parameters, both initial-set forms and supplied enclosures.
They verify preservation after reconditioning and splitting, parameters used
only in auxiliary outputs, zero-time evaluation, and reject
inconsistent/missing bindings and unsupported projections/recombination, and
retain the existing numerical capability gates.
