import pyariadne_evolution as a

def test_evolution_bindings():
    names = ("Semantics", "Enclosure", "VectorFieldSimulator",
             "IteratedMapEvolver", "IteratedMapEvolverConfiguration", "VectorFieldEvolver", "HybridTime",
             "HybridEnclosure", "HybridSimulator")
    for name in names:
        assert hasattr(a, name), f"Missing evolution binding: {name}"
    assert a.Semantics.UPPER != a.Semantics.LOWER


def test_iterated_map_configuration():
    from pyariadne_evolution import IteratedMap, IteratedMapEvolver, RealVariable, next
    x = RealVariable("x")
    evolver = IteratedMapEvolver(IteratedMap([next(x) << x / 2]))
    config = evolver.configuration()
    config.set_maximum_enclosure_radius(1.0)
    config.set_enable_subdivisions(False)
    assert str(config)
