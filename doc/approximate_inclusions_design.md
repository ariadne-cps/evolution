# Structured ODAI approximate evolution

Ordinary and hybrid approximate evolvers require an explicitly selected input
sampling policy for ODAI. Only structured effective fields are supported;
interval-Taylor multifunctions remain rejected because sampling their interval
coefficients does not select an admissible trajectory in general.

The policy receives the provider's normalized input box and a context with
trajectory index, accepted continuous-step index, global time and duration.
One input vector is chosen per RK4 step and held through all four stages.
Guards/invariants/resets use state, external time and constant parameters;
inputs are not added to physical state and do not leak into event arguments.
A reset consumes no continuous-step index. A new input is selected for the next
continuous step, including after changing input dimension or location.

MidpointInputSampling is an explicit constant-center choice.
UniformInputSampling(seed) uses a counter-based 64-bit mixer keyed by seed,
trajectory, step and coordinate, independent of workload scheduling. Policies
are cloned into the evolver and must be immutable/thread-safe. Output arity,
finiteness and membership in the numerical box are checked before integration.
Parameterized ODAI retain the same explicit RealValuation across steps/modes.

ContinuousDynamics becomes shared, paradigm-independent metadata for ODAE/ODAI.
Its ownership and views are unchanged; solvers and policies belong to evolvers.
Approximate ODAI construction requires a policy; the two-argument orbit signature
still returns sampled curves and cannot request UPPER/LOWER. These samples do
not provide a validated reach or within-step contact certificate.

Python's existing ODAE constructor binding is ported to the shared metadata
change. New sampling bindings and full tutorial parity remain the final phase.
Paving is outside this intervention.

## Construction and metadata

```cpp
ContinuousTimeEvolver<ApproximateTag> evolver(
    inclusion, UniformInputSampling(42u), RealValuation({{p, Real(2.5_x)}}));
HybridContinuousTimeEvolver<ApproximateTag> hybrid(
    automaton, MidpointInputSampling(), RealValuation({{p, Real(2.5_x)}}));
```

The original ODAE constructors remain. An approximate ordinary evolver's
SystemType is now ContinuousDynamics, the same owning metadata view used by
the validated evolver. The view has no paradigm tag: it owns declarations,
while the evolvers select numerical behavior and supported operations.

Trajectory indices follow the ordered initial-point list. Reordering that list
changes uniform streams. Step indices count every accepted continuous step,
including ODAE steps, and are retained through zero-time resets. Changing the
step configuration changes the piecewise input signal. Per-step samples are
finite binary64 values; convex interpolation avoids overflow for wide boxes.

Seven regression cases verify policy boxes and centers, nonlinear/time-dependent
flow with fixed parameters and alternating inputs, normalized provider
coordinates, exact numerical reproducibility across serial/parallel execution,
cloning and seeds, hybrid single/multiple paths and transitions that change input
dimension. For x'=p+t+u, y'=u^2 with p=5/2 and u alternating -1,+1 over four
steps of length 1/16, they obtain x(1/4)=21/32 and y(1/4)=1/4. A timed reset
from ODAE enters an autonomous ODAI with a time-only guard, then a two-input
nonlinear ODAI; the final state is 739/128. The tests also reject wrong input
arity, out-of-box/NaN samples, missing policies/parameters and interval-Taylor
representations, including at zero horizon where appropriate.

## Validation and source headers

Local verification passes all 55 CTest suites (38 Evolution including Python,
13 native System and four native Integration), the complete test/binding/
benchmark build, a fresh install and both existing standalone tutorials.
Before commit, apply-format and check-format pass on 109 tracked C++ files.

33 existing file prologues are aligned with the common Evolution format:
complete GPL-3.0-or-later text, separated copyright/banner lines and @file/@brief
metadata. Existing authors and years are retained; two stale renamed-test
banners are corrected. The header audit checks all 109 tracked C++ files,
including public/private headers and tests. New files use the same format.
Diagnostic output continues to use LOGGING. Submodules and paving are unchanged.

Built-in policies select round-to-nearest for interpolation and restore the
caller's rounding mode on return or exception. A regression compares identical
seed/context samples in upward and downward rounding over 32 steps and checks
that each caller mode is preserved. This prevents ambient arithmetic modes from
changing the input signal in the last bits.
