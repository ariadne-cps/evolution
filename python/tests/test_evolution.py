import pyariadne as a

def test_evolution_bindings():
    names = ("Semantics", "Enclosure", "VectorFieldSimulator",
             "IteratedMapEvolver", "VectorFieldEvolver", "HybridTime",
             "HybridEnclosure", "HybridSimulator")
    for name in names:
        assert hasattr(a, name), f"Missing evolution binding: {name}"
    assert a.Semantics.UPPER != a.Semantics.LOWER
