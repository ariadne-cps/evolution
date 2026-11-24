import pyariadne_evolution as a

def test_evolution_bindings():
    names = ("Semantics", "Enclosure", "ContinuousTimeEvolver", "ApproximateTag", "ValidatedTag",
             "OrdinaryDifferentialAlgebraicEquation", "OrdinaryDifferentialAlgebraicInclusion",
             "DiscreteTimeEvolver", "DiscreteTimeValidatedEvolverConfiguration", "HybridTime",
             "HybridEnclosure", "HybridContinuousTimeEvolver", "recombine",
             "EnclosureRecombinerInterface", "EnclosureRecombiner",
             "ReferenceEnclosureRecombiner", "EnclosureRecombinationContext",
             "EnclosureRecombinationResult", "EnclosureRecombinationOutput")
    for name in names:
        assert hasattr(a, name), f"Missing evolution binding: {name}"
    assert a.Semantics.UPPER != a.Semantics.LOWER
    assert not hasattr(a, "VectorField")
    assert not hasattr(a, "DifferentialInclusion")
    for old_name in ("VectorFieldEvolver", "DifferenceAlgebraicEquationEvolver", "GeneralHybridEvolver"):
        assert not hasattr(a, old_name)


def test_paradigm_contracts():
    import unittest
    approximate = a.ContinuousTimeEvolver[a.ApproximateTag]
    validated = a.ContinuousTimeEvolver[a.ValidatedTag]
    assert approximate is not validated
    assert a.DiscreteTimeEvolver[a.ValidatedTag] is a.DiscreteTimeValidatedEvolver
    assert a.HybridContinuousTimeEvolver[a.ApproximateTag] is a.HybridContinuousTimeApproximateEvolver
    for old_name in ("ContinuousTimeSimulator", "HybridContinuousTimeSimulator"):
        assert not hasattr(a, old_name)
    case = unittest.TestCase()
    with case.assertRaises(KeyError):
        a.DiscreteTimeEvolver[a.ApproximateTag]
    x = a.RealVariable("x")
    equation = a.OrdinaryDifferentialAlgebraicEquation([a.dot(x) << x])
    initial = a.RealExpressionBoundedConstraintSet([(a.dec_(0.9) <= x) & (x <= a.dec_(1.1))], [])
    sampled = approximate(equation)
    with case.assertRaises(TypeError):
        sampled.orbit(initial, a.Real(1)/16, a.Semantics.UPPER)
    certified = validated(equation, a.TaylorPicardIntegrator(1e-5))
    with case.assertRaises(TypeError):
        certified.orbit(initial, a.Real(1)/16)


def test_difference_algebraic_equation_configuration():
    from pyariadne_evolution import DifferenceAlgebraicEquation, DiscreteTimeEvolver, ValidatedTag, RealVariable, next
    x = RealVariable("x")
    evolver = DiscreteTimeEvolver[ValidatedTag](DifferenceAlgebraicEquation([next(x) << x / 2]))
    config = evolver.configuration()
    config.set_maximum_enclosure_radius(1.0)
    config.set_enable_subdivisions(False)
    assert str(config)


def test_continuous_recombination():
    x = a.RealVariable("x")
    evolver = a.DiscreteTimeEvolver[a.ValidatedTag](a.DifferenceAlgebraicEquation([a.next(x) << x / 2]))
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
    evolver = a.ContinuousTimeEvolver[a.ValidatedTag](equation, integrator)
    initial = a.RealExpressionBoundedConstraintSet([(0 <= x) & (x <= 0)], [])
    orbit = evolver.orbit(initial, a.Real(1)/16, a.Semantics.UPPER)
    assert orbit.final().size() > 0
    assert orbit.final()[0].state_set().dimension() == 1
    assert orbit.final()[0].state_auxiliary_set().dimension() == 2
    assert isinstance(a.ContinuousTimeEvolver[a.ApproximateTag](equation), a.ContinuousTimeEvolver[a.ApproximateTag])
    timed = a.OrdinaryDifferentialAlgebraicEquation(
        [a.dot(x) << (t - x)], a.RealVariablesBox({}), time=t)
    timed_evolver = a.ContinuousTimeEvolver[a.ValidatedTag](timed, integrator)
    timed_evolver.configuration().set_maximum_step_size(0.03125)
    timed_orbit = timed_evolver.orbit(initial, a.Real(1)/16, a.Semantics.UPPER)
    assert timed_orbit.final().size() > 0
    assert timed_orbit.reach().size() >= 2
    parameterised = a.OrdinaryDifferentialAlgebraicEquation(
        [a.dot(x) << (p - x)], a.RealVariablesBox({p: (0, 1)}))
    parameterised_evolver = a.ContinuousTimeEvolver[a.ValidatedTag](parameterised, integrator)
    parameterised_evolver.configuration().set_maximum_step_size(0.03125)
    parameterised_orbit = parameterised_evolver.orbit(initial, a.Real(1)/16, a.Semantics.UPPER)
    assert parameterised_orbit.final().size() > 0
    assert parameterised_orbit.reach().size() >= 2
    with unittest.TestCase().assertRaises(ValueError):
        parameterised_evolver.orbit(initial, a.Real(1)/16, a.Semantics.LOWER)
    timed_simulator = a.ContinuousTimeEvolver[a.ApproximateTag](timed)
    timed_simulator.configuration().set_step_size(0.03125)
    timed_simulator.orbit(initial, a.Real(1)/16)
    for field in (parameterised,):
        with unittest.TestCase().assertRaises(ValueError):
            a.ContinuousTimeEvolver[a.ApproximateTag](field)


def test_common_inclusion_consumers():
    import unittest
    x, t, p, u = (a.RealVariable(name) for name in ("x", "t", "p", "u"))
    system = a.OrdinaryDifferentialAlgebraicInclusion(
        [a.dot(x) << (t + p + u)], a.RealVariablesBox({u: (-1, 1)}),
        a.RealVariablesBox({p: (0, 1)}), time=t)
    integrator = a.TaylorPicardIntegrator(1e-6)
    selected = a.ContinuousTimeEvolver[a.ValidatedTag](
        system, integrator, a.InputApproximation(a.ZeroApproximation()),
        a.InputApproximation(a.ConstantApproximation()))
    boxes = a.ContinuousTimeEvolver[a.ValidatedTag](system, a.FirstOrderMultifunctionIntegrator())
    initial = a.RealExpressionBoundedConstraintSet([(0 <= x) & (x <= 0)], [])
    for evolver in (selected, boxes):
        evolver.configuration().set_maximum_step_size(0.03125)
        result = evolver.orbit(initial, a.Real(1)/16, a.Semantics.UPPER)
        assert result.final().size() > 0
        assert result.final()[0].state_set().dimension() == 1
        assert result.reach().size() >= 2
        with unittest.TestCase().assertRaises(ValueError):
            evolver.orbit(initial, a.Real(1)/16, a.Semantics.LOWER)


if __name__ == "__main__":
    test_evolution_bindings()
    test_paradigm_contracts()
    test_difference_algebraic_equation_configuration()
    test_continuous_recombination()
    test_ordinary_equation_consumers()
    test_common_inclusion_consumers()
