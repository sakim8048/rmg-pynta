import pytest

from rmgpynta import selection


def _rec(name, A=None, n=None, Ea=None):
    return {"reaction_string": name, "kinetics_A": A, "kinetics_n": n, "kinetics_Ea_kJmol": Ea}


def test_sensitivity_strategy_raises_with_explanation():
    with pytest.raises(selection.SensitivityNotSupported, match="surface_reactor"):
        selection.select_important_reactions([], strategy="sensitivity")


def test_flux_threshold_ranks_by_rate_constant_and_respects_top_n():
    records = [
        _rec("slow", A=1e8, n=0.0, Ea=150.0),
        _rec("fast", A=1e13, n=0.0, Ea=20.0),
        _rec("medium", A=1e10, n=0.0, Ea=80.0),
        _rec("no_kinetics"),  # excluded, not ranked as "slowest"
    ]
    selected = selection.select_by_rate_constant(records, temperature_k=900.0, top_n=2)
    assert [r["reaction_string"] for r in selected] == ["fast", "medium"]


def test_unrankable_reports_missing_kinetics():
    records = [_rec("has_kinetics", A=1e10, n=0.0, Ea=50.0), _rec("missing")]
    missing = selection.unrankable(records)
    assert [r["reaction_string"] for r in missing] == ["missing"]


def test_min_rate_filter():
    records = [_rec("fast", A=1e13, n=0.0, Ea=5.0), _rec("slow", A=1e8, n=0.0, Ea=200.0)]
    selected = selection.select_by_rate_constant(records, temperature_k=900.0, top_n=None, min_rate=1.0)
    assert [r["reaction_string"] for r in selected] == ["fast"]


def test_custom_callable_strategy():
    def keep_none(records, **kwargs):
        return []

    assert selection.select_important_reactions([_rec("x", 1, 0, 1)], strategy=keep_none) == []


def test_unknown_strategy_raises():
    with pytest.raises(ValueError):
        selection.select_important_reactions([], strategy="not_a_real_strategy")
