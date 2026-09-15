"""Step 7: select which edge reactions are worth a Pynta DFT calculation.

Default strategy is a flux/rate-constant threshold over the edge
(``select_by_rate_constant``). A real sensitivity x uncertainty strategy is
also available (``strategy="sensitivity_uncertainty"``,
``select_by_sensitivity_uncertainty`` below) -- despite ``rmgpy/rmg/input.py``
hard-rejecting ``sensitivity=[...]`` for ``surfaceReactor(...)``::

    if sensitivity:
        raise NotImplementedError("Can't currently do sensitivity with surface reactors.")

that restriction is specific to the ``input.py`` convenience wrapper, not to
RMG-Py's solver: ``rmgpy.tools.uncertainty.Uncertainty.sensitivity_analysis()``
constructs ``rmgpy.solver.SurfaceReactor`` directly and it supports
sensitivity fine (confirmed against a real round-0 core from this project --
see ``drivers/run_uncertainty_job.py``, which drives that class from rmg_env
and needed one small upstream RMG-Py patch to get there).
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Callable, Iterable

R_KJ_PER_MOL_K = 8.314462618e-3


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


def _normalize_reaction_string(reaction_string: str) -> str:
    """Match ``drivers.run_uncertainty_job``'s normalization, so a ranking
    JSON's chemkin-format strings ("X(5)+CH4(1)<=>CH4X(13)") line up with
    this project's own plain-label ``reaction_string`` records
    ("X + CH4 <=> CH4X"). Also tolerates already-normalized input (idempotent)
    so records can be matched directly if ever needed.
    """
    lhs, rhs = reaction_string.split("<=>")

    def side(s):
        tokens = re.sub(r"\(\d+\)", "", s).split("+")
        return " + ".join(token.strip() for token in tokens)

    return f"{side(lhs)} <=> {side(rhs)}"


def select_by_sensitivity_uncertainty(
    records: Iterable[dict],
    ranking_path: "str | Path",
    top_n: int | None = 15,
) -> list:
    """Rank by RMG-Py's real sensitivity x uncertainty contribution to a
    chosen observable species, using a ranking JSON that
    ``drivers.run_uncertainty_job`` (run in rmg_env, since it needs
    ``rmgpy.tools.uncertainty``) must already have produced for this round --
    that driver isn't invoked from here, since this module stays pure Python
    (see the package's own subprocess-boundary design in
    ``rmgpynta/__init__.py``); see ``orchestrator.run_uncertainty_analysis``.

    A record's ``reaction_string`` (plain species labels, e.g.
    "X + CH4 <=> CH4X") is matched against the ranking's normalized chemkin
    strings. RMG occasionally has more than one "duplicate" reaction (same
    net stoichiometry, different reacting sites) sharing one
    ``reaction_string`` -- confirmed ~1% of a real core -- so the first
    matching record wins; both represent essentially the same TS-search job
    for Pynta's purposes.
    """
    ranking = json.loads(Path(ranking_path).read_text())
    by_string: dict = {}
    for record in records:
        by_string.setdefault(record["reaction_string"], record)

    selected = []
    for row in ranking:
        record = by_string.get(_normalize_reaction_string(row["reaction_string"]))
        if record is not None:
            selected.append(record)
        if top_n is not None and len(selected) >= top_n:
            break
    return selected


def select_important_reactions(
    edge_records: Iterable[dict],
    strategy: str | Callable = "flux_threshold",
    **kwargs,
) -> list:
    if strategy == "flux_threshold":
        return select_by_rate_constant(edge_records, **kwargs)
    if strategy == "sensitivity_uncertainty":
        return select_by_sensitivity_uncertainty(edge_records, **kwargs)
    if callable(strategy):
        return strategy(edge_records, **kwargs)
    raise ValueError(f"Unknown selection strategy: {strategy!r}")
