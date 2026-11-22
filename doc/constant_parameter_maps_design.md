# Physical parameter maps and continuation

A trajectory-constant model parameter is a named physical value `p=P(s)`;
it is not the current coordinate `s_i`. Store an optional binding in Enclosure
containing ordered parameter names, the declared physical domain, and a
validated map on the enclosure parameter domain. Initial declaration binding
adjoins physical coordinates once. Continuation retains an existing binding
only when names/order agree and its declared domain is contained in the new
declaration domain. Incompatible declarations fail instead of resampling.

The mapping is immutable along a trajectory: flow and parameter extension
embed it, restriction/splitting restrict it, and reconditioning preserves every
coordinate on which it depends. Its uniform rounding errors must become
explicit persistent coordinates before integration or reduction; otherwise a
new error choice could hide a change in the physical parameter. Domain
constraints restrict the physical map to its declared bounds. Parameter-only
dependencies count towards the protected-coordinate budget in evolution.
Taylor parameter maps use a non-discarding sweeper, including coefficients
below the state solver's threshold. Alternative bases retain all coordinates;
exposing nonzero parameter errors currently requires a Taylor representation.

Multi-source recombination requires compatible named bindings. Normalize each
source map along with state, global time and dwell. Bound differences of the
joint `(state,time,dwell,physical parameters)` map against the reference map,
and add residual coordinates for all components. The existing constructive
coverage witness extends componentwise to physical parameters. Output domain
metadata is the hull of declared domains, and explicit constraints keep the
physical values within that hull. Recombination remains an outer operation:
it can introduce extra state/parameter combinations and grants no LOWER claim.

Raw PARAMETER axes without a binding cannot identify the physical tuple after
normalization and are rejected for multi-source recombination. Single-source
identity results remain valid. A Cartesian state extension by an interval or
box can transport one binding. A product of two Enclosure operands needs an
explicit reconciliation policy when either carries a binding.

Numerical field/provider domain checks and existing integrator capability
gates still apply. Compatible declaration metadata alone cannot validate
evaluation outside a certified interval/Taylor model domain. Hybrid parameter
evolution, approximate input/parameter policies, converters, and final Python
tutorial parity remain separate roadmap activities.

The implemented regressions verify analytic continuation after normalization on non-unit
physical domains, zero-horizon continuation without duplicate parameters,
nonzero global time, correlated initial state/parameters, reconditioning,
splitting/restriction, explicit parameter rounding errors, serial/parallel
recombination, explicit feasible witnesses from both source domains, multiple
parameter components/order, ODAI continuation and incompatible-declaration
rejection. Continuation under a wider declaration preserves the existing
parameter choices instead of widening them.

Paving changes are owned by the separate `decoupling` branch and are outside
this intervention. Python tutorial parity remains the final Evolution activity.
