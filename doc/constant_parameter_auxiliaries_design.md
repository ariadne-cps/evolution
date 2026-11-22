# Constant parameters in ordinary validated auxiliary outputs

Ordinary validated equation evolution already carries trajectory-constant
parameter coordinates with `EnclosureVariableKind::PARAMETER`. Add auxiliary
layouts `[state, parameters]` and `[state, global time, parameters]` with an
explicit parameter count. Coordinate selection uses the last declared
`PARAMETER` axes, in their existing order, and never interprets INPUT, ERROR
or TEMPORAL coordinates as constant parameters. Numerical integration and
auxiliary evaluation share that selection.

Attach the system auxiliary mapping after adjoining the model parameter
domain, for box, constrained-set and enclosure initial data. The clock remains
external and the physical state dimension is unchanged. A missing parameter
coordinate or inconsistent layout/count is an error. Reconditioning preserves
PARAMETER axes and their physical intervals; splitting/restriction keeps the
same explicit auxiliary binding.

This step supports UPPER ordinary ODAE evolution with compatible integrators.
Existing LOWER, preconditioned, hybrid-parameter and approximate-sampling
capability gates remain. It does not define persistent labelled parameter
bindings for resuming an orbit under a different system declaration.

A state-only augmented paving cannot represent an auxiliary value depending
on an external parameter or clock and rejects that conversion. The current
multi-source recombiner normalizes source coordinates and cannot preserve
this physical-parameter binding; reject parameter-dependent auxiliary
recombination explicitly, while a single-source identity result remains
valid. A future general parameter-map representation must preserve physical
parameter values through changes of basis before this restriction is lifted.
Joint state/time/auxiliary enclosure evaluation and drawing remain available.

The implemented regressions compare correlated auxiliary maps against analytic
solutions for autonomous and timed equations, including non-unit parameter
domains, multiple parameters, both initial-set forms and supplied enclosures.
They verify preservation after reconditioning and splitting, parameters used
only in auxiliary outputs, zero-time evaluation, and reject
inconsistent/missing bindings and unsupported projections/recombination, and
retain the existing numerical capability gates.
