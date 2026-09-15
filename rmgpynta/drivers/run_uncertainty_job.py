#!/usr/bin/env python
"""Alternative to selection.py's rate-constant ranking: rank reactions by their
sensitivity x uncertainty contribution to a chosen observable species, using
RMG-Py's own rmgpy.tools.uncertainty.Uncertainty class.

Runs inside rmg_env:
  <rmg_env_python> run_uncertainty_job.py <chemkin_gas> <chemkin_surface> <dictionary> <output_dir> [extra_libraries_json]

``extra_libraries_json`` (optional) is a JSON file with
``{"extra_thermo_libraries": [...], "extra_reaction_libraries": [...]}`` --
the same round-tagged library names ``orchestrator.run_rmg`` spliced into
that round's ``input.py`` (``libraries.update_input_py_libraries``'s
``extra_thermo_libs``/``extra_kinetics_libs``). Required from round 1
onward: without it, ``extract_sources_from_model()`` won't recognize any
species/reaction whose thermo or kinetics actually came from a Pynta-derived
library registered in an earlier round, and raises on the resulting
unexpected source key.

REQUIRES a local patch to RMG-Py itself (~/RMG-Py/RMG-Py/rmgpy/data/thermo.py,
ThermoDatabase.get_thermo_data): without it, extract_sources_from_model() reliably
crashes with "ValueError: Method only valid for radicals." the first time it hits
a species whose overall radical character comes from one resonance structure while
another resonance structure in the same species.molecule list is a non-radical
valence tautomer (e.g. a singlet carbene) -- confirmed against this project's own
round 0 core, species 'NC=[C]'. estimate_radical_thermo_via_hbi requires
molecule.is_radical(), but the loop that calls it over every species.molecule entry
only checked molecule.reactive, not molecule.get_radical_count() > 0. The fix adds
that check, mirroring the function's own existing reactive-only filtering.

Why this needs to bypass ``surfaceReactor(...)``'s own sensitivity=[...] option:
RMG-Py's rmgpy/rmg/input.py hard-rejects it for surface reactors
(NotImplementedError("Can't currently do sensitivity with surface reactors.")).
Uncertainty.sensitivity_analysis() sidesteps that entirely by constructing
rmgpy.solver.SurfaceReactor directly, which supports sensitivity fine -- the
restriction only lives in the input-file convenience wrapper, not the solver.

Why matching by RMG reaction index doesn't work: Uncertainty.load_model() parses
the chemkin FILE into brand-new Reaction objects via rmgpy.chemkin.load_chemkin_file,
so their .index reflects chemkin file order, not the original generation-time index
this project's own core_reactions.json/edge_reactions.json record as "rmg_index".
Confirmed by direct comparison (reaction "X + CH4 <=> CH4X" is rmg_index 7 in our
JSON but chemkin-file-order index 326 here). Matching instead goes by a normalized
reaction string (species labels only, index-in-parens and "+"/"<=>" spacing
stripped) -- reaction_schema.py's records already use exactly that plain-label
format ("X + CH4 <=> CH4X"), confirmed to round-trip cleanly against real records.
"""
from __future__ import annotations

import json
import re
import shutil
import sys
from pathlib import Path

from rmgpy.tools.plot import parse_csv_data
from rmgpy.tools.uncertainty import Uncertainty

KINETICS_FAMILIES = [
    "Surface_Adsorption_vdW", "Surface_Adsorption_Dissociative", "Surface_Adsorption_Single",
    "Surface_Adsorption_Double", "Surface_Adsorption_Bidentate",
    "Surface_Dissociation", "Surface_Abstraction", "Surface_Abstraction_vdW",
    "default",
]
THERMO_LIBRARIES = ["surfaceThermoPt111", "primaryThermoLibrary", "DFT_QCI_thermo"]

# Mirrors rmg_input_template.py's surfaceReactor systems 1 (CH4) and 2 (NH3) --
# same T/P/termination/coverage/surface-volume-ratio, so the sensitivity run
# reflects the same conditions the core was actually generated under.
OBSERVABLE_CONDITIONS = {
    "CH4": {"mole_fraction": 0.10},
    "NH3": {"mole_fraction": 0.10},
}
T = (900, "K")
P = (1.0, "bar")
TERMINATION_TIME = (1.0, "s")
SURFACE_VOLUME_RATIO = (1.0e5, "m^-1")


def _normalize_reaction_string(chemkin_reaction_string: str) -> str:
    """'X(5)+CH4(1)<=>CH4X(13)' -> 'X + CH4 <=> CH4X', matching
    rmgpynta.drivers.run_rmg_job._reaction_record's plain-label format.
    """
    lhs, rhs = chemkin_reaction_string.split("<=>")

    def side(s):
        tokens = re.sub(r"\(\d+\)", "", s).split("+")
        return " + ".join(tokens)

    return f"{side(lhs)} <=> {side(rhs)}"


def _rank_reaction_contributions(unc: Uncertainty, sens_species, reaction_system_index: int = 0) -> list:
    """Reimplements the ranking core of Uncertainty.local_analysis() (each
    reaction's (sensitivity x uncertainty)**2 contribution to output variance)
    without its plotting side effect, and returns normalized reaction strings
    instead of the chemkin-file-order index local_analysis's own return value
    carries (see module docstring for why that index isn't usable here).
    """
    csvfile_path = Path(unc.output_directory) / "solver" / f"sensitivity_{reaction_system_index + 1}_SPC_{sens_species.index}.csv"
    _time, data_list = parse_csv_data(str(csvfile_path))
    ranked = []
    for data in data_list:
        if not data.reaction:
            continue
        rxn_index = int(data.index) - 1
        uncertainty = unc.kinetic_input_uncertainties[rxn_index]
        contribution = (data.data[-1] * uncertainty) ** 2
        ranked.append({
            "reaction_string": _normalize_reaction_string(data.reaction),
            "contribution": contribution,
        })
    ranked.sort(key=lambda r: r["contribution"], reverse=True)
    return ranked


def main(
    chemkin_gas: str,
    chemkin_surface: str,
    dictionary_path: str,
    output_dir: str,
    extra_libraries_json: str | None = None,
) -> None:
    output_dir = str(Path(output_dir).resolve())
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    extra_thermo_libraries: list = []
    extra_reaction_libraries: list = []
    if extra_libraries_json:
        extra = json.loads(Path(extra_libraries_json).read_text())
        extra_thermo_libraries = extra.get("extra_thermo_libraries", [])
        extra_reaction_libraries = extra.get("extra_reaction_libraries", [])

    unc = Uncertainty(output_directory=output_dir)
    unc.load_model(chemkin_path=chemkin_gas, dictionary_path=dictionary_path, surface_path=chemkin_surface)
    unc.load_database(
        kinetics_families=KINETICS_FAMILIES,
        kinetics_depositories=["training"],
        thermo_libraries=THERMO_LIBRARIES + extra_thermo_libraries,
        reaction_libraries=extra_reaction_libraries,
    )
    unc.extract_sources_from_model()
    unc.assign_parameter_uncertainties()

    by_label = {}
    for s in unc.species_list:
        by_label.setdefault(s.label, s)
    ar = by_label["Ar"]
    x = by_label["X"]

    # Best score seen per reaction across every observable species it was run
    # against, plus which species drove that score -- a reaction only needs to
    # matter for *one* target observable to be worth Pynta's time.
    best_by_reaction: dict = {}
    for name, cond in OBSERVABLE_CONDITIONS.items():
        sens_species = by_label[name]
        mole_fraction = cond["mole_fraction"]
        unc.sensitivity_analysis(
            initial_mole_fractions={sens_species: mole_fraction, ar: 1.0 - mole_fraction},
            sensitive_species=[sens_species],
            T=T, P=P, termination_time=TERMINATION_TIME,
            initial_surface_coverages={x: 1.0},
            surface_volume_ratio=SURFACE_VOLUME_RATIO,
        )
        # rmgpy.util.make_output_subdirectory (called by sensitivity_analysis
        # itself, every time) shutil.rmtree()s output_dir/solver before
        # writing to it -- confirmed against a real run: it silently destroyed
        # CH4's own sensitivity/thermo CSV+PNG the moment the NH3 iteration
        # started, though the ranking numbers below were already safely read
        # into memory by then, so only the on-disk diagnostics were lost, not
        # the actual selection. Copy this species' solver/ output aside
        # before the next observable's sensitivity_analysis() call wipes it.
        solver_copy = Path(output_dir) / f"solver_{name}"
        if solver_copy.exists():
            shutil.rmtree(solver_copy)
        shutil.copytree(Path(output_dir) / "solver", solver_copy)

        for row in _rank_reaction_contributions(unc, sens_species):
            key = row["reaction_string"]
            if key not in best_by_reaction or row["contribution"] > best_by_reaction[key]["contribution"]:
                best_by_reaction[key] = {**row, "sensitive_to": name}

    ranking = sorted(best_by_reaction.values(), key=lambda r: r["contribution"], reverse=True)
    Path(output_dir, "sensitivity_uncertainty_ranking.json").write_text(json.dumps(ranking, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:6])
