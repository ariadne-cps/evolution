# Coverage verification notes

The coverage target measures the public Evolution headers and implementation
files using the normal test suite. Tests assert numeric values, bounds,
parameter dimensions, event histories, stopping conditions and public contract
failures. The watertank workload is a separate benchmark.

## Enclosure representation invariants

The former private `Enclosure::_check` repeatedly tested map dimensions after
operations that already establish them. Public constructors validate incoming
argument dimensions; map application also checks result compatibility with the
auxiliary mapping. Parameterised steps check both map dimensions and domains
before mutation. Constraint entry points check arity. Function factories create
maps on the supplied domain; composition and restriction preserve the parameter
arity. Parameter extension embeds the time, dwell and constraint maps while
appending the state coordinate and its variable kind together. Products and
reconditioning rebuild these corresponding members together.

These facts make the old post-operation dimension-failure branches unreachable
for supported function factories and successful public operations. They have
been removed. The public legacy `HybridEnclosure::_check` symbol remains as a
compatibility hook; the contained labelled enclosure enforces its contracts at
creation and update boundaries.

The reduced domain starts equal to the parameter domain. Restriction intersects
it with a subdomain, extension takes the same Cartesian product for both domains,
and reconditioning constructs corresponding new domains before installing the
new representation. The pinned hull contractor modifies a coordinate only via
`max(old_lower,candidate_lower)` and `min(old_upper,candidate_upper)`
(`function/procedure.cpp` and `procedure.tpl.hpp`). `FloatDP` implements these
operations using `std::max`/`std::min` with the old endpoint first. Thus a valid
endpoint remains valid even if an intermediate candidate is NaN, and contraction
preserves inclusion. The redundant NaN restoration and subset-failure checks
have been removed on this basis.

The local Taylor reconditioning helper receives sorted, distinct parameter
indices selected from `[0,argument_size)`. It is called only after establishing
that at least one parameter must be discarded. Its error index is either a state
component index below the number of error variables or that number itself (the
uniform-error sentinel). Assertions for contrary internal states were redundant.

Generic `Function` domains are `EuclideanDomain(argument_size)`; their dimension
checks therefore follow from the preceding arity preconditions. Function-patch
domains are bounded boxes and retain explicit domain compatibility checks.

## Template and macro interpretation

LLVM can emit constructor entry variants and multiple template instantiations
whose individual counters differ. The official function/line/branch summary is
the completion criterion; structured per-symbol data and macro expansions are
used to explain residuals. Region coverage is recorded but is not a target.
