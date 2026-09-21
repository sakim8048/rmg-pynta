"""Corrected reimplementation of pynta.postprocessing.write_rmg_libraries.

Runs inside pynta_env. Standalone: imports only ``pynta`` and stdlib.

Real, confirmed bug in ``pynta/postprocessing.py`` (``write_rmg_libraries``,
verified against the installed checkout at ``~/pynta``): the loop that
builds ``reaction_text`` for ``reaction_library/reactions.py`` is guarded by
``if reaction_text == "":``, and the body-append + index increment sit
*inside* that guard --

    for ts,kinetics in ts_dict.items():
        ...
        if reaction_text == "":
            reaction_text += kin.create_RMG_header(...)
            reaction_text += kin.rmg_kinetics_text.replace("{index}",str(index))
            index += 1
            reaction_text += "\\n"

-- so after the first transition state's kinetics is appended,
``reaction_text`` is no longer ``""`` and every subsequent TS in
``ts_dict`` is silently skipped. The library never has more than one
reaction in it, however many TS the run actually computed. This
reimplements the same logic with the header written once and every valid
TS's kinetics appended (matching how the thermo-library loop right above it
in the same function already behaves correctly). This fix should also be
proposed to zadorlab/pynta directly, separately from this project.
"""
from __future__ import annotations

import os
import re
import shutil

from pynta.postprocessing import postprocess

# Real, confirmed bug in pynta/postprocessing.py: SurfaceConfiguration /
# GasConfiguration fit their NASA polynomial with wh.to_nasa(Tmin=298.15,
# ...), which literally renders as "Tmin=(298.15,'K')" (twice -- once on
# the first NASAPolynomial segment, once on the outer NASA object) in
# rmg_species_text. But rmgpy.thermo.nasa.NASA.to_thermo_data() (called
# from ThermoDatabase.correct_binding_energy while loading a library's
# training rules) hardcodes H298/S298 lookups at exactly 298 K, not
# 298.15 K -- 298 < 298.15 falls outside the segment's own bounds, so
# NASA.select_polynomial raises "ValueError: No valid NASA polynomial at
# temperature 298 K." Confirmed against a real round-1 RMG rerun using a
# Pynta-registered thermo library. Pynta's own placeholder/zero-Cp species
# template (same file) already writes "Tmin = (298.0, 'K')" for exactly
# this reason, so 298.0 is the convention to match, not 298.15.
_NASA_TMIN_298_15_RE = re.compile(r"Tmin=\(298\.15,'K'\)")


def _fix_nasa_tmin(rmg_species_text: str) -> str:
    return _NASA_TMIN_298_15_RE.sub("Tmin=(298.0,'K')", rmg_species_text)


def _lowest_energy(spc):
    """spc may be a single config object, or (per pynta's own, somewhat
    inconsistent handling in the original write_rmg_libraries) a dict or
    list of candidate conformers keyed/indexed by search order. Handle both
    shapes defensively rather than assume pynta's original code's apparent
    dict-vs-list inconsistency was intentional.
    """
    if isinstance(spc, dict):
        valid = {k: v for k, v in spc.items() if getattr(v, "valid", True)}
        if not valid:
            return None
        return spc[min(valid, key=lambda k: spc[k].energy)]
    if isinstance(spc, list):
        valid = [v for v in spc if getattr(v, "valid", True)]
        if not valid:
            return None
        return min(valid, key=lambda v: v.energy)
    return spc


_CRYSTAL_STRUCTURE_PREFIXES = ("fcc", "bcc", "hcp", "sc", "diamond", "rocksalt", "hexagonal")


def _rmg_facet(surface_type: str) -> str:
    """Translate Pynta's ``surface_type`` (e.g. ``"fcc111"``, ``"bcc110"``,
    ``"hcp0001"`` -- crystal structure + Miller index) into the facet label
    RMG's own surface binding-energy database actually uses.

    Real, confirmed bug: RMG-database/input/surface/libraries/metal.py keys
    its entries as ``metal + Miller index`` only (``"Pt111"``, ``"Pt211"``,
    ...), never with the crystal-structure prefix. ``ThermoDatabase.
    get_thermo_data`` builds its scaling lookup as ``db_label = entry.metal
    + entry.facet`` (rmgpy/data/thermo.py) from whatever this module writes
    into a library entry's ``metal``/``facet`` comment fields -- passing
    Pynta's raw ``surface_type`` straight through as ``facet`` produced
    ``db_label = "Pt" + "fcc111" = "Ptfcc111"``, which RMG's surface.py
    raises ``DatabaseError: Metal 'Ptfcc' not found in database`` on.
    Confirmed against a real round-1 RMG rerun using a Pynta-registered
    thermo library (metal="Pt", surface_type="fcc111").
    """
    for prefix in _CRYSTAL_STRUCTURE_PREFIXES:
        if surface_type.startswith(prefix):
            return surface_type[len(prefix):]
    return surface_type


def _lowest_barrier(kinetics):
    if isinstance(kinetics, dict):
        valid = {k: v for k, v in kinetics.items() if getattr(v, "valid", True)}
        if not valid:
            return None
        return kinetics[min(valid, key=lambda k: kinetics[k].barrier_f)]
    if isinstance(kinetics, list):
        valid = [v for v in kinetics if getattr(v, "valid", True)]
        if not valid:
            return None
        return min(valid, key=lambda v: v.barrier_f)
    return kinetics


def write_rmg_libraries_fixed(path, metal, facet, repeats, sites, site_adjacency, slab_path):
    """Reruns postprocess() and writes a *complete* RMG kinetics library.

    ``sites``/``site_adjacency``/``slab_path`` are required (not
    regenerated here) -- pass through the same values the just-completed
    ``pynta.main.Pynta`` job instance already computed
    (``job.sites``, ``job.site_adjacency``, ``job.slab_path`` after
    ``job.execute(...)`` -- see ``drivers/run_pynta_job.py``), since
    ``postprocess()`` needs exactly this data and re-deriving it
    independently here was not verified against a confirmed ACAT API.

    Returns the number of reactions actually written, so the caller can
    sanity-check it against ``len(ts_dict)`` rather than assume success.
    """
    spc_dict, ts_dict, spc_dict_thermo = postprocess(
        path, metal, facet, sites, site_adjacency, repeats=repeats, slab_path=slab_path, check_finished=True
    )

    # --- thermo_library.py: unchanged from pynta's own (correct) logic ---
    index = 0
    thermo_text = ""
    for _name, spc in spc_dict_thermo.items():
        minspc = _lowest_energy(spc)
        if minspc is None:
            continue
        if thermo_text == "":
            thermo_text += minspc.create_RMG_header("thermo_library", "", "", metal, _rmg_facet(facet))
        thermo_text += _fix_nasa_tmin(minspc.rmg_species_text.replace("{index}", str(index)))
        index += 1
        thermo_text += "\n"
    with open(os.path.join(path, "thermo_library.py"), "w") as f:
        f.write(thermo_text)

    # --- dictionary.txt: unchanged from pynta's own logic ---
    spc_dictionary_txt = "vacantX\n1 X u0 p0 c0\n\n"
    for name, spc in spc_dict.items():
        minspc = _lowest_energy(spc)
        if minspc is None:
            continue
        spc_dictionary_txt += name + "\n"
        spc_dictionary_txt += minspc.mol.to_adjacency_list()
        spc_dictionary_txt += "\n"

    # --- reactions.py: THE FIX -- append every valid TS, not just the first ---
    index = 0
    reaction_text = ""
    n_written = 0
    for _ts, kinetics in ts_dict.items():
        kin = _lowest_barrier(kinetics)
        if kin is None:
            continue
        if reaction_text == "":
            reaction_text += kin.create_RMG_header("reaction_library", lib_short_desc="", lib_long_desc="")
        reaction_text += kin.rmg_kinetics_text.replace("{index}", str(index))
        index += 1
        reaction_text += "\n"
        n_written += 1

    rxn_lib_dir = os.path.join(path, "reaction_library")
    if os.path.exists(rxn_lib_dir):
        shutil.rmtree(rxn_lib_dir)
    os.makedirs(rxn_lib_dir)
    with open(os.path.join(rxn_lib_dir, "reactions.py"), "w") as f:
        f.write(reaction_text)
    with open(os.path.join(rxn_lib_dir, "dictionary.txt"), "w") as f:
        f.write(spc_dictionary_txt)

    return n_written
