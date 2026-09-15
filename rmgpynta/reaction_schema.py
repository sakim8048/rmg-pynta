"""Translate RMG reaction records into Pynta's ``reaction.yaml`` schema (steps 2 & 8).

This module is pure Python (no ``rmgpy`` or ``pynta`` import) and operates on
plain-dict "reaction records" -- the intermediate JSON that
``rmgpynta.drivers.run_rmg_job`` dumps from inside ``rmg_env`` after an RMG
run, where it has live ``Reaction``/``Species`` objects to pull adjacency
lists and family names from directly.

Reaction record schema (one dict per reaction, from the intermediate JSON)::

    {
      "chemkin_index": int,             # "Chemkin #7" in the annotated .inp comment
      "rmg_index": int,                 # "RMG #9" in the annotated .inp comment
      "reaction_string": str,           # e.g. "X(5)+CH4(1)<=>CH4X(13)"
      "family": str,                    # RMG kinetics family, e.g. "Surface_Adsorption_Dissociative"
      "degeneracy": float,
      "reactants": [{"label": str, "adjacency_list": str}, ...],  # per-species, for the "reaction: X + Y => Z" string only
      "products":  [{"label": str, "adjacency_list": str}, ...],  # ditto
      "reactant_block": str,            # this side's full Pynta-format block -- see reaction_record_to_pynta_entry
      "product_block": str,             # ditto
      "kinetics_comment": str,          # raw multi-line "! ..." comment block, for provenance
    }

Pynta ``reaction.yaml`` entry schema, confirmed against real files under
``~/pynta-production/*.yaml`` (e.g. ``hb-yesdiffusion.yaml``)::

    - index: 0
      reaction: "N2 + 2X => 2 NX"
      reaction_family: "Dissociative Adsorption"
      reactant: |
        multiplicity 1
        1 *1 N u0 p1 c0 {3,T}
        2 *2 N u0 p1 c0 {4,T}
        3 *3 X u0 p0 c0 {1,T}
        4 *4 X u0 p0 c0 {2,T}
      product: |
        multiplicity 1
        1 *1 N u0 p1 c0 {2,T}
        2 *2 N u0 p1 c0 {1,T}
        3 *3 X u0 p0 c0
        4 *4 X u0 p0 c0

Note ``index`` is cosmetic -- ``pynta.main.Pynta.__init__`` re-stamps every
entry's ``index`` by enumeration order when it loads ``rxns_file``
(``pynta/main.py``, the ``for i,r in enumerate(rxns_list): r["index"] = i``
block), so anything written here is only for human readability while
debugging a generated yaml.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import yaml

# RMG kinetics-family -> Pynta reaction_family. The first three entries were
# confirmed against families actually observed producing a matched Pynta job
# in this conversation's prior artifacts (Pt(111)/Pt(557)/Pt(211) audits) and
# in ~/pynta-production/*.yaml.
#
# The rest were added after tracing every use of reaction_family/family_name
# through ~/pynta/pynta/*.py: it is read in exactly two places
# (transitionstate.py's generate_constraints_harmonic_parameters, where the
# parameter is accepted but never referenced again, and postprocessing.py,
# where only the literal string "Surface_Migration" changes behavior -- it
# switches on a diffusion-specific validation path). Everywhere else it is
# inert: stored in info.json and echoed into a human-readable comment. The
# real TS-search strategy comes from structurally comparing reactant/product
# bonds (get_broken_formed_bonds + the SubgraphIsomorphicDecisionTree), not
# from this string. So for every family below (none of which is literally
# "Surface_Migration"), passing the RMG family name straight through is safe
# and matches how real runs under ~/pynta-production/*.yaml already recorded
# them (e.g. "Surface_Abstraction_vdW" and "Surface_Adsorption_Dissociative"
# appear verbatim as reaction_family values in files there).
RMG_TO_PYNTA_FAMILY = {
    "Surface_Adsorption_Dissociative": "Dissociative Adsorption",
    "Surface_Dissociation": "Dissociation",
    "Surface_Abstraction": "Surface Abstraction",
    "Surface_Adsorption_vdW": "Surface_Adsorption_vdW",
    "Surface_Adsorption_Single": "Surface_Adsorption_Single",
    "Surface_Adsorption_Double": "Surface_Adsorption_Double",
    "Surface_Adsorption_Bidentate": "Surface_Adsorption_Bidentate",
    "Surface_Abstraction_vdW": "Surface_Abstraction_vdW",
}

UNMAPPED_RMG_FAMILIES_SEEN_IN_INPUT_PY = (
    # From ~/rmg-production/rmg_ch4_pt111/input.py's kineticsFamilies list.
    # Still deliberately excluded from rmg_input_template.py -- both crashed
    # RMG itself with "Could not process a species with no reactive
    # structures" (see run.log.attempt1/2), independent of the Pynta-mapping
    # question, so there's nothing to map yet.
    "Surface_Bidentate_Dissociation",
    "Surface_Dissociation_to_Bidentate",
    "Surface_Monodentate_to_Bidentate",
    "Surface_vdW_to_Bidentate",
    "Surface_EleyRideal_Addition_Multiple_Bond",
)


def is_surface_family(rmg_family: str) -> bool:
    """Whether an RMG kinetics family involves the catalyst surface at all.

    RMG's own family-set naming convention (see RMG-database's
    input/kinetics/families/recommended.py: the ``surface`` and
    ``surface_development`` sets) prefixes every heterogeneous family with
    ``Surface_``; the gas-phase ``default`` set has none. Used to keep
    homogeneous gas-phase reactions (``H_Abstraction``, ``R_Recombination``,
    etc.) out of Pynta rounds -- Pynta's TS-search pipeline is built around
    adsorbed species and slab sites, so a reaction with no surface site has
    nothing for it to build a TS structure from.
    """
    return rmg_family.startswith("Surface_")


class UnsupportedReactionFamily(ValueError):
    """Raised when a reaction record's RMG family has no confirmed Pynta mapping."""


class AtomLabelsNotVerified(ValueError):
    """Raised when a reaction record's reacting-atom (``*1``..``*4``) labels
    were not confirmed recoverable on the rmg_env side.

    Pynta's reactant/product adjacency-list text labels the specific atoms
    that participate in the reaction (see the module docstring's real
    example: ``1 *1 N u0 p1 c0 {3,T}``). A generic RMG ``Species``'s
    ``molecule.to_adjacency_list()`` does not carry these labels -- they
    live on the family-generated ``Reaction``'s template atoms. Recovering
    them correctly (which atom is ``*1`` vs ``*2``, consistently between
    the reactant and product block) was not independently verified while
    writing this module (see ``drivers/run_rmg_job.py``'s extraction code).
    Shipping an unlabeled or mislabeled block to Pynta risks a TS search
    that either fails outright or silently explores the wrong bond -- so
    this is a hard failure by default rather than a best-effort guess.
    Pass ``allow_unverified_labels=True`` to override once you've confirmed
    the labeling logic against a real, working example.
    """


def translate_family(rmg_family: str) -> str:
    try:
        return RMG_TO_PYNTA_FAMILY[rmg_family]
    except KeyError:
        raise UnsupportedReactionFamily(
            f"RMG kinetics family {rmg_family!r} has no confirmed Pynta "
            f"reaction_family mapping. Confirmed mappings: "
            f"{sorted(RMG_TO_PYNTA_FAMILY)}. Add one to "
            f"RMG_TO_PYNTA_FAMILY in rmgpynta/reaction_schema.py only after "
            f"verifying it against a real, working Pynta run -- do not guess."
        ) from None


def reaction_record_to_pynta_entry(
    record: dict, index: int, allow_unverified_labels: bool = False
) -> dict:
    """Convert one intermediate reaction record into a Pynta rxns_file entry."""
    if not record.get("atom_labels_verified", False) and not allow_unverified_labels:
        raise AtomLabelsNotVerified(
            f"reaction record {record.get('reaction_string')!r} (family "
            f"{record.get('family')!r}) has no verified *1..*4 atom labels; "
            f"see rmgpynta.reaction_schema.AtomLabelsNotVerified."
        )
    lhs = " + ".join(s["label"] for s in record["reactants"])
    rhs = " + ".join(s["label"] for s in record["products"])
    return {
        "index": index,
        "reaction": f"{lhs} => {rhs}",
        "reaction_family": translate_family(record["family"]),
        "reactant": record["reactant_block"],
        "product": record["product_block"],
    }


def write_pynta_reactions_yaml(
    records: Iterable[dict], out_path: Path, allow_unverified_labels: bool = False
) -> Path:
    """Translate a list of intermediate reaction records into ``reaction.yaml``.

    Any record whose family is unmapped raises ``UnsupportedReactionFamily``,
    and any record without verified atom labels raises
    ``AtomLabelsNotVerified`` (see that class's docstring) -- both fail
    loud rather than silently dropping or mis-labeling a reaction, so a
    round never ships Pynta a smaller or chemically-wrong job than
    RMG/selection intended.
    """
    entries = [
        reaction_record_to_pynta_entry(r, i, allow_unverified_labels=allow_unverified_labels)
        for i, r in enumerate(records)
    ]
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        yaml.dump(entries, f, default_flow_style=False, sort_keys=False, width=100)
    return out_path


def synthesize_diffusion_entries(species_records: Iterable[dict], template_yaml: Path) -> list:
    """Build Pynta ``Diffusion``-family entries for adsorbed core species.

    RMG's own ``kineticsFamilies`` (see ~/rmg-production/rmg_ch4_pt111/input.py)
    has no site-to-site diffusion family -- RMG generates bond-forming/breaking
    chemistry, not a species hopping between otherwise-identical empty sites.
    Pynta's own reaction.yaml files *do* carry ``Diffusion``-family entries
    (confirmed: ``reaction_family: Diffusion`` in real yaml files such as
    ``hb-yesdiffusion.yaml``), so those have to be synthesized directly from
    each round's adsorbed core species rather than translated out of an RMG
    Reaction object.

    NOT IMPLEMENTED YET: doing this correctly means reproducing the exact
    ``*1..*4`` atom-labeling convention a diffusion entry uses (which atom
    hops, which site label is "from" vs "to"), and no verified real example
    of that adjacency-list text was available while writing this module --
    only the human-readable ``"NX + X => NX + X"`` string form (see the
    Pt(557) Kinetics Audit artifact). Guessing the label pattern risks
    writing a reaction.yaml Pynta's molecule parser rejects, or worse,
    accepts but mis-assigns the reacting atom. Before implementing this:
    open one real ``Diffusion``-family entry from
    ``~/pynta-production/*/hb-yesdiffusion.yaml`` (or pass its path as
    ``template_yaml`` here) and hand-verify the labeling against a working
    Pynta TS-search run, then encode that as a template substitution below.
    """
    raise NotImplementedError(
        "synthesize_diffusion_entries: no verified real Diffusion-entry "
        "adjacency-list example was available during this design pass -- "
        "see this function's docstring before implementing."
    )
