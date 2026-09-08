#!/usr/bin/env python
"""Step 6 (+ extraction for steps 7/8): run RMG, dump core/edge reactions to JSON.

Runs inside rmg_env: ``<rmg_env_python> run_rmg_job.py <input.py> <output_dir>``

Writes, into ``<output_dir>``:
  - everything RMG itself already writes there (``chemkin/``, ``seed/``,
    ``RMG.log``, ...) via its own normal ``RMG.execute()`` behavior
  - ``core_reactions.json`` / ``edge_reactions.json`` -- this driver's own
    dump of ``rmg.reaction_model.core.reactions`` / ``.edge.reactions``,
    in the intermediate schema ``rmgpynta.reaction_schema`` expects.

Standalone by design: imports only ``rmgpy`` and stdlib -- never
``rmgpynta`` itself.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from rmgpy.rmg.main import RMG  # confirmed construction pattern: see rmgpy/__main__.py's
# `rmg = RMG(input_file=args.file, output_directory=args.output_directory); rmg.execute(**kwargs)`


def _species_record(species) -> dict:
    """species.molecule[0].to_adjacency_list() is confirmed against pynta's
    own postprocessing.py, which calls the equivalent `.mol.to_adjacency_list()`
    on its own (compatible) Molecule objects.

    IMPORTANT, UNRESOLVED: this adjacency list is numbered from 1 within
    this species alone, and carries NO `*1`/`*2`/... reacting-atom labels --
    those only exist on a KineticsFamily's *template* atoms during reaction
    generation, not on the instantiated product/reactant Species objects
    here. Recovering which atom is *1 vs *2 (and keeping that consistent
    between the reactant and product side of one reaction) needs the
    family's reaction "recipe"/action list (bond BREAK/FORM/CHANGE actions
    with their template atom labels), which was not independently
    confirmed against a working example while writing this driver. Until
    that's implemented and verified, every record this driver emits is
    marked "atom_labels_verified": False, and
    rmgpynta.reaction_schema.write_pynta_reactions_yaml refuses to build a
    reaction.yaml from an unverified record unless explicitly overridden --
    see rmgpynta.reaction_schema.AtomLabelsNotVerified.
    """
    return {
        "label": species.label,
        "adjacency_list": species.molecule[0].to_adjacency_list(),
    }


def _reaction_record(rxn) -> dict:
    kinetics = getattr(rxn, "kinetics", None)
    record = {
        "rmg_index": getattr(rxn, "index", None),
        "chemkin_index": None,  # only assigned when chemkin is saved; matched up in provenance.py instead
        "reaction_string": (
            " + ".join(s.label for s in rxn.reactants) + " <=> " + " + ".join(s.label for s in rxn.products)
        ),
        "family": getattr(rxn, "family", None),
        "degeneracy": getattr(rxn, "degeneracy", 1.0),
        "reactants": [_species_record(s) for s in rxn.reactants],
        "products": [_species_record(s) for s in rxn.products],
        "kinetics_comment": getattr(kinetics, "comment", "") if kinetics is not None else "",
        "atom_labels_verified": False,  # see _species_record's docstring
    }
    # Best-effort surface-Arrhenius parameters for rmgpynta.selection's
    # rate-constant ranking, if this reaction has fitted kinetics at all
    # (edge reactions frequently don't, until/unless RMG explores them
    # further). Attribute names (.A.value_si etc.) follow rmgpy's standard
    # Arrhenius/SurfaceArrhenius Quantity convention but were not
    # independently re-confirmed against this specific checkout -- wrapped
    # defensively so a mismatch degrades to "no kinetics" rather than
    # crashing the whole extraction.
    try:
        record["kinetics_A"] = kinetics.A.value_si
        record["kinetics_n"] = kinetics.n.value_si
        record["kinetics_Ea_kJmol"] = kinetics.Ea.value_si / 1000.0
    except AttributeError:
        record["kinetics_A"] = record["kinetics_n"] = record["kinetics_Ea_kJmol"] = None
    return record


def main(input_file: str, output_dir: str) -> None:
    output_dir = str(Path(output_dir).resolve())
    rmg = RMG(input_file=input_file, output_directory=output_dir)
    rmg.execute()  # also triggers RMG's own chemkin/seed output as a normal side effect

    core = [_reaction_record(r) for r in rmg.reaction_model.core.reactions]
    edge = [_reaction_record(r) for r in rmg.reaction_model.edge.reactions]

    Path(output_dir, "core_reactions.json").write_text(json.dumps(core, indent=2))
    Path(output_dir, "edge_reactions.json").write_text(json.dumps(edge, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
