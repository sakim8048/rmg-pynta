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
import shutil

from pynta.postprocessing import postprocess


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
            thermo_text += minspc.create_RMG_header("thermo_library", "", "", metal, facet)
        thermo_text += minspc.rmg_species_text.replace("{index}", str(index))
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
