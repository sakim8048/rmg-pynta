"""Step 7: select which edge reactions are worth a Pynta DFT calculation.

Real constraint confirmed against ``rmgpy/rmg/input.py`` on this machine:
``surface_reactor(...)`` (what ``surfaceReactor(...)`` in an RMG ``input.py``
maps to) explicitly does this on a ``sensitivity`` argument::

    if sensitivity:
        raise NotImplementedError("Can't currently do sensitivity with surface reactors.")

Every real RMG job for this project (see the ``rmg_ch4_pt111/input.py``
referenced throughout the RMG-Pynta Feedback Loop artifact) uses
``surfaceReactor(...)``, not ``simpleReactor(...)``/``liquidReactor(...)``
(the two reactor types that do support ``sensitivity=[...]``). So RMG's
built-in sensitivity analysis -- what the original design sketch assumed
would be available -- cannot actually be used for this system. The default
strategy here is a flux/rate-constant threshold over the edge instead.
"""
from __future__ import annotations

import math
from typing import Callable, Iterable

R_KJ_PER_MOL_K = 8.314462618e-3


class SensitivityNotSupported(NotImplementedError):
    pass


def _rate_constant(record: dict, temperature_k: float) -> float | None:
    """k = A * T^n * exp(-Ea / R T), from an edge record's kinetics fields.

    Returns None if the record doesn't carry A/n/Ea (e.g. no kinetics was
    ever fit for this edge reaction, which happens for some estimation
    paths) -- such records are excluded from ranking rather than assigned
    an arbitrary rate.
    """
    A = record.get("kinetics_A")
    n = record.get("kinetics_n")
    Ea = record.get("kinetics_Ea_kJmol")
    if A is None or n is None or Ea is None:
        return None
    try:
        return A * (temperature_k**n) * math.exp(-Ea / (R_KJ_PER_MOL_K * temperature_k))
    except (OverflowError, ValueError):
        return None


def select_by_rate_constant(
    edge_records: Iterable[dict],
    temperature_k: float = 900.0,
    top_n: int | None = 25,
    min_rate: float | None = None,
) -> list:
    """Approximate importance ranking: fastest edge reactions at the
    reactor's operating temperature.

    This is a deliberately simple proxy, not a substitute for a true flux
    or leak-rate analysis (RMG computes internal edge-species leak fluxes
    during model generation, but a confirmed, stable API for pulling those
    numbers back out per-reaction was not available while writing this
    module -- see the TODO in ``drivers/run_rmg_job.py``). Records with no
    fitted kinetics are excluded rather than ranked arbitrarily; call
    ``unrankable(edge_records)`` to inspect what got left out.
    """
    ranked = []
    for r in edge_records:
        k = _rate_constant(r, temperature_k)
        if k is not None:
            ranked.append((k, r))
    ranked.sort(key=lambda pair: pair[0], reverse=True)
    if min_rate is not None:
        ranked = [(k, r) for k, r in ranked if k >= min_rate]
    if top_n is not None:
        ranked = ranked[:top_n]
    return [r for _, r in ranked]


def unrankable(edge_records: Iterable[dict]) -> list:
    return [r for r in edge_records if _rate_constant(r, 900.0) is None]


def select_important_reactions(
    edge_records: Iterable[dict],
    strategy: str | Callable = "flux_threshold",
    **kwargs,
) -> list:
    if strategy == "flux_threshold":
        return select_by_rate_constant(edge_records, **kwargs)
    if strategy == "sensitivity":
        raise SensitivityNotSupported(
            "RMG-Py's surface_reactor(...) raises NotImplementedError for "
            "sensitivity=[...] (rmgpy/rmg/input.py, surface_reactor, the "
            "'Can't currently do sensitivity with surface reactors' check) "
            "-- sensitivity analysis only works with simpleReactor/"
            "liquidReactor, neither of which represents this Pt "
            "surface-chemistry system. Use strategy='flux_threshold' "
            "(the default), or supply your own callable."
        )
    if callable(strategy):
        return strategy(edge_records, **kwargs)
    raise ValueError(f"Unknown selection strategy: {strategy!r}")
