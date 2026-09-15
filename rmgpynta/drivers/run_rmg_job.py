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
import logging
import sys
from pathlib import Path

from rmgpy.rmg.main import RMG  # confirmed construction pattern: see rmgpy/__main__.py's
# `rmg = RMG(input_file=args.file, output_directory=args.output_directory); rmg.execute(**kwargs)`


def _species_record(species) -> dict:
    """species.molecule[0].to_adjacency_list() is confirmed against pynta's
    own postprocessing.py, which calls the equivalent `.mol.to_adjacency_list()`
    on its own (compatible) Molecule objects.

    If ``_label_reacting_atoms`` below successfully relabeled this species'
    molecule, ``to_adjacency_list()`` includes the family's real `*1`/`*2`/...
    labels on the atoms whose bonds this reaction actually breaks or forms.
    Otherwise this is the unlabeled fallback, numbered from 1 within this
    species alone -- see ``_reaction_record``'s ``atom_labels_verified``.

    This per-species text is only used for the reaction JSON's own
    "reactants"/"products" list (labels, individual reference); the combined
    reactant_block/product_block Pynta actually reads comes from
    ``_merged_adjacency_block`` instead, since concatenating per-species text
    directly reintroduces the bug that fixed (see its docstring).
    """
    return {
        "label": species.label,
        "adjacency_list": species.molecule[0].to_adjacency_list(),
    }


def _merged_adjacency_block(species_list) -> str:
    """Combine one reaction side's species into Pynta's single-block adjacency
    format -- one shared ``multiplicity 1`` header, then every atom of every
    species numbered continuously with bonds resolved against that shared
    numbering -- confirmed against real files under
    ``~/pynta-production/*.yaml`` (e.g. ``hb-yesdiffusion.yaml``).

    Concatenating each species' own (independently-numbered)
    ``to_adjacency_list()`` text, as this driver used to do indirectly via
    ``rmgpynta.reaction_schema._adjacency_block``, produces two bugs on any
    multi-species side: each species' own multiplicity line lands stray in
    the middle of the block, and atom numbering restarts at 1 per species
    instead of continuing, leaving bond references like ``{3,S}`` ambiguous
    between species. ``rmgpy.molecule.molecule.Molecule.merge`` (confirmed in
    RMG-Py's own source) combines molecules into one ``Molecule`` with a
    single continuous atom list, renumbering every bond reference to match --
    exactly the real RMG-side atom-mapping this driver's schema docstring
    said the renumbering needed.

    The header's multiplicity has to be the real one, not a placeholder:
    Pynta re-derives it right after parsing (``reactants.multiplicity =
    reactants.get_radical_count() + 1`` in pynta/tasks.py), but
    ``Molecule.from_adjacency_list`` -- the same RMG-Py parser Pynta calls --
    validates the stated multiplicity against Hund's rule *during* parsing
    and raises ``InvalidAdjacencyListError`` before Pynta's override ever
    runs if it's wrong (confirmed: a merged block with an unpaired-electron
    atom, common for abstraction intermediates, is rejected outright when
    hardcoded to "multiplicity 1"). ``Molecule.merge`` doesn't compute one
    for the merged object on its own, so it's set here with the same formula
    Pynta itself uses.
    """
    merged = species_list[0].molecule[0]
    for species in species_list[1:]:
        merged = merged.merge(species.molecule[0])
    merged.multiplicity = merged.get_radical_count() + 1
    lines = merged.to_adjacency_list().splitlines()
    if lines and lines[0].strip().startswith("multiplicity"):
        lines = lines[1:]
    return f"multiplicity {merged.multiplicity}\n" + "\n".join(lines) + "\n"


def _label_reacting_atoms(rxn, kinetics_families: dict) -> bool:
    """Relabel ``rxn``'s reactant/product species with the family's real
    `*1`/`*2`/... reacting-atom labels, via
    ``KineticsFamily.add_atom_labels_for_reaction`` (rmgpy/data/kinetics/
    family.py) -- confirmed against RMG-Py's own source: it mutates a
    reaction's reactant/product Species objects in place, labeling only the
    atoms whose bonds the family's recipe actually breaks/forms/changes,
    consistently between the reactant and product side.

    Works on deep copies of ``rxn.reactants``/``rxn.products`` rather than
    the live objects: those Species objects are shared across every other
    reaction in ``rmg.reaction_model`` that also involves them, and
    ``add_atom_labels_for_reaction`` narrows a shared species' ``.molecule``
    list down to the one resonance structure it just labeled -- mutating the
    live objects would silently corrupt labeling for any other
    not-yet-processed reaction sharing that species. Returns False (leaving
    the unlabeled fallback) if this family/reaction combination can't be
    labeled -- e.g. no matching family, or no isomorphic template match.
    """
    family = kinetics_families.get(getattr(rxn, "family", None))
    if family is None:
        return False
    try:
        rxn.reactants = [s.copy(deep=True) for s in rxn.reactants]
        rxn.products = [s.copy(deep=True) for s in rxn.products]
        family.add_atom_labels_for_reaction(rxn)
        return True
    except Exception:
        logging.warning("Could not recover atom labels for reaction %s (family %s)", rxn, family.label)
        return False


def _reaction_record(rxn, kinetics_families: dict) -> dict:
    atom_labels_verified = _label_reacting_atoms(rxn, kinetics_families)
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
        "reactant_block": _merged_adjacency_block(rxn.reactants),
        "product_block": _merged_adjacency_block(rxn.products),
        "kinetics_comment": getattr(kinetics, "comment", "") if kinetics is not None else "",
        "atom_labels_verified": atom_labels_verified,
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

    kinetics_families = rmg.database.kinetics.families
    core = [_reaction_record(r, kinetics_families) for r in rmg.reaction_model.core.reactions]
    edge = [_reaction_record(r, kinetics_families) for r in rmg.reaction_model.edge.reactions]

    Path(output_dir, "core_reactions.json").write_text(json.dumps(core, indent=2))
    Path(output_dir, "edge_reactions.json").write_text(json.dumps(edge, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
