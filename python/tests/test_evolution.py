import pyariadne_evolution as a

def test_evolution_bindings():
    names = ("Semantics", "Enclosure", "VectorFieldSimulator",
             "OrdinaryDifferentialAlgebraicEquation", "OrdinaryDifferentialAlgebraicInclusion",
             "DiscreteTimeEvolver", "DiscreteTimeEvolverConfiguration", "ContinuousTimeEvolver", "HybridTime",
             "HybridEnclosure", "HybridSimulator", "recombine",
             "EnclosureRecombinerInterface", "EnclosureRecombiner",
             "ReferenceEnclosureRecombiner", "EnclosureRecombinationContext",
             "EnclosureRecombinationResult", "EnclosureRecombinationOutput")
    for name in names:
        assert hasattr(a, name), f"Missing evolution binding: {name}"
    assert a.Semantics.UPPER != a.Semantics.LOWER
    assert not hasattr(a, "VectorField")
    assert not hasattr(a, "DifferentialInclusion")
    for old_name in ("VectorFieldEvolver", "IteratedMapEvolver", "GeneralHybridEvolver"):
        assert not hasattr(a, old_name)


def test_iterated_map_configuration():
    from pyariadne_evolution import IteratedMap, DiscreteTimeEvolver, RealVariable, next
    x = RealVariable("x")
    evolver = DiscreteTimeEvolver(IteratedMap([next(x) << x / 2]))
    config = evolver.configuration()
    config.set_maximum_enclosure_radius(1.0)
    config.set_enable_subdivisions(False)
    assert str(config)


def test_continuous_recombination():
    x = a.RealVariable("x")
    evolver = a.DiscreteTimeEvolver(a.IteratedMap([a.next(x) << x / 2]))
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
    assert str(handle) == "ParallelReferenceEnclosureRecombiner()"
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


def test_ordinary_equation_consumers():
    import unittest
    x, t, p, output = (a.RealVariable(name) for name in ("x", "t", "p", "output"))
    equation = a.OrdinaryDifferentialAlgebraicEquation(
        [a.dot(x) << 1], [a.let(output) << (2*x)])
    integrator = a.TaylorPicardIntegrator(1e-5)
    evolver = a.ContinuousTimeEvolver(equation, integrator)
    initial = a.RealExpressionBoundedConstraintSet([(0 <= x) & (x <= 0)], [])
    orbit = evolver.orbit(initial, a.Real(1)/16, a.Semantics.UPPER)
    assert orbit.final().size() > 0
    assert orbit.final()[0].state_set().dimension() == 1
    assert orbit.final()[0].state_auxiliary_set().dimension() == 2
    assert isinstance(a.VectorFieldSimulator(equation), a.VectorFieldSimulator)
    unsupported = (
        a.OrdinaryDifferentialAlgebraicEquation(
            [a.dot(x) << (t - x)], a.RealVariablesBox({}), time=t),
        a.OrdinaryDifferentialAlgebraicEquation(
            [a.dot(x) << (p - x)], a.RealVariablesBox({p: (0, 1)})),
    )
    for field in unsupported:
        with unittest.TestCase().assertRaises(ValueError):
            a.ContinuousTimeEvolver(field, integrator)
        with unittest.TestCase().assertRaises(ValueError):
            a.VectorFieldSimulator(field)


if __name__ == "__main__":
    test_evolution_bindings()
    test_iterated_map_configuration()
    test_continuous_recombination()
    test_ordinary_equation_consumers()
