# Ariadne Evolution tutorial

The C++ tutorial is an installed-library consumer: **its only header is
`#include "ariadne-evolution.hpp"`**. All public system declarations,
numerical backends and supporting types needed by its examples must be
available through that umbrella header. The standalone CMake project links
only the installed `AriadneEvolution::ariadne-evolution` target.

## C++ worked examples

Examples are intentionally selected for distinct **implemented** C++
contracts and syntax, not as a Cartesian product of system families,
paradigms and set semantics.

| Worked example | Syntax it illustrates |
| --- | --- |
| Autonomous ODAE, sampled and validated | `ContinuousTimeEvolver<ApproximateTag>` with two-argument `orbit`; `ContinuousTimeEvolver<ValidatedTag>` with `Semantics::UPPER` and `Semantics::LOWER` |
| Nonautonomous ODAE, sampled | Explicit `TimeVariable` in a continuous field |
| Ordinary structured ODAI, sampled and validated | Explicit `MidpointInputSampling` vs validated full-input `Semantics::UPPER` |
| Ordinary difference equation, validated | `DiscreteTimeEvolver<ValidatedTag>` with both UPPER and LOWER |
| Ordinary difference inclusion, sampled and validated | `DiscreteTimeEvolver<ApproximateTag>` with explicit input policy vs `DiscreteTimeInclusionEvolver` with UPPER |
| Hybrid discrete map and reset alternatives | `HybridDiscreteTimeEvolver<ApproximateTag>` branch-complete sampled paths vs `HybridDiscreteTimeEvolver<ValidatedTag>` with an integer termination criterion and UPPER |
| Hybrid discrete inclusion with auxiliary event | One `ValidatedTag` permissive transition where the same `a=x+u` determines guard, progression invariant and reset; whole-box input correlations are retained |
| Continuous hybrid modes and flow | `ContinuousTimeHybridAutomaton`, approximate topology, validated short equation and inclusion evolutions |
| First-order PDE | One bounded two-dimensional hyperbolic problem using `FirstOrderPDE` and `first_order_pde` |

The same variables, initial sets and dynamics are reused where possible.
The PDE appears once; the other combinations are represented by the minimal
additional calls needed to show genuinely different signatures. Numerical
samples are not set enclosures, and `Semantics::LOWER` is not an INNER
approximation.

The Python tutorial demonstrates the currently exposed installed
`pyariadne_evolution` workflow: a validated difference equation and a
continuous hybrid transition model. It is kept separate from the C++ header
contract.

## Build and install

From the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPYTHON_ENFORCE_TESTS=ON
cmake --build build --parallel --target everything
cmake --install build --prefix install
```

## C++ tutorial

Build independently against the installed library:

```bash
cmake -S tutorial -B tutorial/build -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build tutorial/build --parallel
./tutorial/build/tutorial_evolution
```

On Windows the executable is `tutorial/build/tutorial_evolution.exe`.
Depending on the platform, configure the dynamic library search path to
include `install/lib` (or `install/bin` on Windows).

## Python tutorial

Point `PYTHONPATH` to the installed Python extension:
`install/lib/pythonX.Y/site-packages` on Unix, or
`install/Lib/site-packages` on Windows. Configure the library search path,
then run:

```bash
python tutorial/tutorial_evolution.py
```

Neither script requires pytest or graphics. Release CI builds and runs both
tutorials against installed artifacts on Ubuntu, macOS and Windows.
