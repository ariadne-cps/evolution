# Ariadne Evolution tutorial

The tutorials are non-graphical. The C++ tutorial exercises **five**
evolution capabilities with bounded scenarios:

1. **ContinuousTimeSimulator** — simulate a stable one-dimensional continuous
   vector field for 0.1 time units.
2. **ContinuousTimeEvolver** — compute a validated upper-semantics enclosure
   for the same field over 0.05 time units, then evolve an inclusion with
   measurable uncertain inputs using one explicitly selected reach method
   and one endpoint method.
3. **DiscreteTimeEvolver** — take one validated step of a contracting map.
4. **HybridContinuousTimeSimulator and HybridContinuousTimeEvolver** — specify a two-mode
   discrete-only automaton, create a simulator, then compute a short
   validated orbit for a separate continuous single-mode hybrid system,
   followed by an ODAI mode using explicit reach and endpoint approximations.
5. **FirstOrderPDE** — construct a validated, owning PDE model using
   `add_direction()` and run the actual two-dimensional first-order PDE
   solver and inspect its computed solution history. This is the most
   expensive section; the initial profile uses unit-frequency sinusoids.

The Python program builds a contracting `IteratedMap`, computes one validated
upper-semantics iteration using `DiscreteTimeEvolver`, and constructs the same
discrete-only `HybridAutomaton` used in C++. Both examples use the installed
`pyariadne_evolution` extension. `DiscreteTimeEvolverConfiguration` is exposed so that
its enclosure-radius and subdivision options can be configured.

As the final Evolution migration task, the Python tutorial will execute the
same workflows as C++, including `ContinuousTimeSimulator`,
`ContinuousTimeEvolver`, `HybridContinuousTimeEvolver`, and the `FirstOrderPDE`
solver. Any missing model, integrator or solver bindings will be added then,
after the underlying APIs have settled.

## Build and install

From the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPYTHON_ENFORCE_TESTS=ON
cmake --build build --parallel --target everything
cmake --install build --prefix install
```

## C++ tutorial

Build **separately** against the installed library:

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

Neither script requires pytest or graphics. The Release CI builds and runs
**both tutorials against the installed artifacts** on Ubuntu, macOS, and
Windows, to catch packaging and binding regressions.
