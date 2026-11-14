import pyariadne_evolution as a

def test_evolution_bindings():
    names = ("Semantics", "Enclosure", "VectorFieldSimulator",
             "IteratedMapEvolver", "IteratedMapEvolverConfiguration", "VectorFieldEvolver", "HybridTime",
             "HybridEnclosure", "HybridSimulator", "recombine",
             "EnclosureRecombinerInterface", "EnclosureRecombiner",
             "ReferenceEnclosureRecombiner", "EnclosureRecombinationContext",
             "EnclosureRecombinationResult", "EnclosureRecombinationOutput")
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


def test_continuous_recombination():
    x = a.RealVariable("x")
    evolver = a.IteratedMapEvolver(a.IteratedMap([a.next(x) << x / 2]))
    evolver.configuration().set_maximum_enclosure_radius(1.0)
    initial = a.RealExpressionBoundedConstraintSet(
        [(a.dec_(0.9) <= x) & (x <= a.dec_(1.1))], [])
    orbit = evolver.orbit(initial, a.Integer(1), a.Semantics.UPPER)
    source = orbit.final()[0]
    result = a.recombine([source, source])
    assert isinstance(result, a.Enclosure)
    assert result.state_set() is not None
    assert isinstance(a.recombine([source]), a.Enclosure)
    concrete = a.ReferenceEnclosureRecombiner()
    handle = a.EnclosureRecombiner(concrete)
    handle = a.EnclosureRecombiner(handle)
    assert isinstance(concrete.clone(), a.ReferenceEnclosureRecombiner)
    del concrete
    batch = handle.recombine([source, source])
    assert len(batch.outputs) == 1
    assert batch.outputs[0].source_indices == [0, 1]
    assert isinstance(batch.outputs[0].enclosure, a.Enclosure)
    assert batch.discarded_sources == []
    assert batch.radius_limit_satisfied
    assert str(handle) == "ReferenceEnclosureRecombiner()"
    assert a.EnclosureRecombinationContext().maximum_enclosure_radius is None
    small = a.EnclosureRecombinationContext(0.01)
    assert small.maximum_enclosure_radius == 0.01
    constrained = handle.recombine([source, source], small)
    assert not constrained.radius_limit_satisfied
    assert constrained.outputs[0].source_indices == [0, 1]
    assert handle.recombine([source], a.EnclosureRecombinationContext(1.0)).radius_limit_satisfied
    assert handle.recombine([]).outputs == []
    assert handle.recombine([]).discarded_sources == []
    try:
        a.recombine([])
    except RuntimeError:
        pass
    else:
        raise AssertionError("An empty source list must be rejected")


if __name__ == "__main__":
    test_evolution_bindings()
    test_iterated_map_configuration()
    test_continuous_recombination()
