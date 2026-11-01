# Ariadne Evolution tutorial

Both examples are intentionally short and non-graphical. The C++ program
illustrates three kinds of system:

1. **Continuous** — simulate a stable one-dimensional vector field over
   0.1 time units, producing one trajectory.
2. **Discrete** — perform one upper-semantics iteration of a contracting map,
   producing a final enclosure.
3. **Hybrid** — construct a two-mode, discrete-only automaton and create a
   `HybridSimulator` with a `HybridTime` budget of one time unit and two
   discrete jumps. No expensive hybrid orbit is computed.

The Python program constructs the **same hybrid automaton**, instantiates
`HybridSimulator`, and checks `HybridTime` and evolution semantics. Its
continuous and discrete classes are reported but not instantiated: the
Python bindings for the System-side `VectorField` and `IteratedMap` models
do not yet expose their constructors. This is a binding limitation, not a
claim that those systems can be simulated in Python today.

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
