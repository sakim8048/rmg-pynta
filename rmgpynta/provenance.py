"""Parse RMG's annotated chemkin comments into a per-reaction provenance log (design doc §03).

RMG already records, as ``! ``-prefixed comment lines directly above each
reaction in its annotated chemkin output (``chemkin/chem_annotated-*.inp``,
``chemkin/chem_edge_annotated-*.inp`` -- see ``rmgpy.chemkin.save_chemkin_files``),
exactly how that reaction's kinetics were estimated: a direct library/seed
match, a matched training reaction, a rate-rule tree-node average, or a
template-based rate-rule estimate. This module turns those comments into a
structured, queryable record instead of leaving them as text a human has to
reread after every round.

The regexes below are built directly from real annotated-chemkin output at
``~/rmg-production/rmg_ch4_pt111/C2N_reactions/chemkin/chem_annotated-gas.inp``
and ``chem_annotated-surface.inp`` -- every comment-line shape they match was
observed verbatim in that file, not guessed from RMG-Py's documentation.
Species/thermo comments follow a similar ``! ``-prefixed convention in the
THERMO section, but no real example was available while writing this module
-- ``classify_species_thermo_comment`` is a best-effort pass over the
phrasings RMG-Py is documented to use and should be checked against a real
thermo block before its output is trusted the way the reaction-side output
can be.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterator

_REACTION_INDEX_RE = re.compile(r"^!\s*Reaction index:\s*Chemkin\s*#(\d+);\s*RMG\s*#(\d+)\s*$")
_TEMPLATE_REACTION_RE = re.compile(r"^!\s*Template reaction:\s*(\S+)\s*$")
_FLUX_PAIRS_RE = re.compile(r"^!\s*Flux pairs:\s*(.+?)\s*$")
_FAMILY_RE = re.compile(r"^!\s*family:\s*(\S+)\s*$")
_LIBRARY_RE = re.compile(r"^!\s*Library reaction:\s*(.+?)\s*$")
_SEED_RE = re.compile(r"^!\s*Seed mechanism:\s*(.+?)\s*$")
_MATCHED_REACTION_RE = re.compile(
    r"^!\s*Matched reaction\s+(\d+)\s+(.+?)\s+in\s+(\S+?)/training\s*$"
)
_MATCHED_RATE_RULE_RE = re.compile(r"^!\s*This reaction matched rate rule\s*\[(.+?)\]\s*$")
_AVERAGE_OF_RE = re.compile(r"^!\s*Average of\s*\[(.+?)\]\s*$")
_ESTIMATED_TEMPLATE_RE = re.compile(
    r"^!\s*Estimated using template\s*\[(.+?)\]\s*for rate rule\s*\[(.+?)\]\s*$"
)
_ESTIMATED_NODE_RE = re.compile(r"^!\s*Estimated from node\s+(\S+)\s+in family\s+(\S+?)\.?\s*$")
_DEGENERACY_RE = re.compile(r"^!\s*Multiplied by reaction path degeneracy\s*([\d.]+)\.?\s*$")

# Best-effort only -- see module docstring. Not verified against a real
# thermo block in this repo.
_THERMO_LIBRARY_RE = re.compile(r"^!\s*Thermo library:\s*(.+?)\s*$")
_THERMO_GAV_RE = re.compile(r"^!\s*Thermo group additivity estimation", re.IGNORECASE)

_REACTION_LINE_RE = re.compile(r"^\s*\S+\s*(?:<=>|=>)\s*\S+.*$")


def _iter_comment_blocks(lines: list) -> Iterator[tuple]:
    """Yield (comment_lines, reaction_line, trailer_lines) for each reaction record.

    ``trailer_lines`` captures modifier lines chemkin writes directly after a
    reaction (e.g. ``STICK``, ``DUPLICATE``) which carry no provenance
    information themselves but are kept for completeness.
    """
    block: list = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i].rstrip("\n")
        stripped = line.strip()
        if stripped.startswith("!"):
            block.append(line)
            i += 1
            continue
        if not stripped:
            i += 1
            continue
        if _REACTION_LINE_RE.match(stripped) and ("<=>" in stripped or "=>" in stripped):
            reaction_line = stripped
            i += 1
            trailer = []
            while i < n and lines[i].strip() and not lines[i].strip().startswith("!") and (
                "<=>" not in lines[i] and "=>" not in lines[i]
            ):
                # A bare modifier line (STICK, DUPLICATE, low-T/high-T Arrhenius
                # continuation, etc.) has no '=>' and isn't itself a new
                # reaction -- swallow it as a trailer of the reaction above.
                stripped_next = lines[i].strip()
                if stripped_next.upper() in {"STICK", "DUPLICATE", "MWON", "FORD", "REV"} or re.match(
                    r"^[\d.eE+\-]+(\s+[\d.eE+\-]+)*$", stripped_next
                ):
                    trailer.append(stripped_next)
                    i += 1
                else:
                    break
            yield block, reaction_line, trailer
            block = []
            continue
        # Any other non-comment, non-reaction, non-trailer line (section
        # keywords like REACTIONS/END, unit declarations, etc.) resets the
        # pending comment block rather than attaching it to something later.
        block = []
        i += 1


def _classify_kinetics(comment_lines: list) -> dict:
    record = {
        "chemkin_index": None,
        "rmg_index": None,
        "template_family": None,
        "family": None,
        "flux_pairs": None,
        "degeneracy": None,
        "kinetics_source": "unknown",
        "library_name": None,
        "training_reaction": None,
        "rate_rule_path": None,
        "raw_comment": "\n".join(comment_lines),
    }
    matched_reaction = None
    matched_rate_rule = None
    average_of = None
    estimated_template = None
    estimated_node = None
    library = None
    seed = None

    for line in comment_lines:
        if m := _REACTION_INDEX_RE.match(line):
            record["chemkin_index"] = int(m.group(1))
            record["rmg_index"] = int(m.group(2))
        elif m := _TEMPLATE_REACTION_RE.match(line):
            record["template_family"] = m.group(1)
        elif m := _FLUX_PAIRS_RE.match(line):
            pairs = [p.strip() for p in m.group(1).split(";") if p.strip()]
            record["flux_pairs"] = pairs
        elif m := _FAMILY_RE.match(line):
            record["family"] = m.group(1)
        elif m := _DEGENERACY_RE.match(line):
            record["degeneracy"] = float(m.group(1))
        elif m := _LIBRARY_RE.match(line):
            library = m.group(1)
        elif m := _SEED_RE.match(line):
            seed = m.group(1)
        elif m := _MATCHED_REACTION_RE.match(line):
            matched_reaction = {"number": int(m.group(1)), "text": m.group(2), "family": m.group(3)}
        elif m := _MATCHED_RATE_RULE_RE.match(line):
            matched_rate_rule = m.group(1)
        elif m := _AVERAGE_OF_RE.match(line):
            average_of = m.group(1)
        elif m := _ESTIMATED_TEMPLATE_RE.match(line):
            estimated_template = {"template": m.group(1), "rate_rule": m.group(2)}
        elif m := _ESTIMATED_NODE_RE.match(line):
            estimated_node = {"node": m.group(1), "family": m.group(2)}

    if record["family"] is None:
        record["family"] = record["template_family"]

    # Precedence follows how directly each phrasing ties the kinetics to a
    # specific, checkable source: an explicit library/seed entry is the
    # strongest claim, a matched training reaction is next (one real
    # reaction was found and reused, possibly averaged with a rate-rule
    # node), then a template/rate-rule estimate, then a bare tree-node
    # average with no matched reaction at all.
    if library is not None:
        record["kinetics_source"] = "library"
        record["library_name"] = library
    elif seed is not None:
        record["kinetics_source"] = "seed_mechanism"
        record["library_name"] = seed
    elif matched_reaction is not None:
        record["kinetics_source"] = "training_reaction"
        record["training_reaction"] = matched_reaction
        record["rate_rule_path"] = matched_rate_rule
    elif average_of is not None:
        record["kinetics_source"] = "rate_rule_average"
        record["rate_rule_path"] = average_of
    elif estimated_template is not None:
        record["kinetics_source"] = "rate_rule"
        record["rate_rule_path"] = estimated_template["rate_rule"]
    elif estimated_node is not None:
        record["kinetics_source"] = "rate_rule"
        record["rate_rule_path"] = estimated_node["node"]
    # else: stays "unknown" -- e.g. a bare reaction line with no comment
    # block above it at all (shouldn't happen with verboseComments=True,
    # see the real input.py's options(verboseComments=True) in the design
    # artifact's referenced RMG job, but handled rather than assumed away).

    return record


def classify_species_thermo_comment(comment_text: str) -> dict:
    """Best-effort only -- see module docstring. Verify before relying on this."""
    record = {"thermo_source": "unknown", "library_name": None, "raw_comment": comment_text}
    for line in comment_text.splitlines():
        if m := _THERMO_LIBRARY_RE.match(line.strip()):
            record["thermo_source"] = "library"
            record["library_name"] = m.group(1)
            return record
        if _THERMO_GAV_RE.match(line.strip()):
            record["thermo_source"] = "group_additivity"
            return record
    return record


def parse_annotated_chemkin(path: Path) -> list:
    """Parse one ``chem_annotated*.inp`` / ``chem_edge_annotated*.inp`` file.

    Returns a list of provenance records, one per reaction found in the
    file's ``REACTIONS`` section (in file order). Species/THERMO parsing is
    intentionally out of scope here -- this targets the reaction-kinetics
    provenance that steps 5-7 of the design actually need to reason about
    the model each round.
    """
    text = Path(path).read_text()
    lines = text.splitlines(keepends=True)
    try:
        start = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("REACTIONS"))
    except StopIteration:
        start = 0
    try:
        end = next(i for i, l in enumerate(lines) if l.strip().upper() == "END" and i > start)
    except StopIteration:
        end = len(lines)

    records = []
    for comment_lines, reaction_line, trailer in _iter_comment_blocks(lines[start:end]):
        record = _classify_kinetics(comment_lines)
        record["reaction_string"] = reaction_line.split()[0] if reaction_line else None
        record["reaction_line"] = reaction_line
        record["trailer"] = trailer
        record["source_file"] = str(path)
        records.append(record)
    return records


def write_provenance_log(records: list, out_path: Path) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(records, f, indent=2)
    return out_path


def summarize(records: list) -> dict:
    """The one rollup number worth watching round over round (design §03):
    what fraction of reactions are still rate-rule/group-additivity
    estimated rather than backed by a library (i.e. by DFT, once Pynta's
    libraries are registered -- see rmgpynta.libraries).
    """
    if not records:
        return {"n_reactions": 0}
    counts: dict = {}
    for r in records:
        counts[r["kinetics_source"]] = counts.get(r["kinetics_source"], 0) + 1
    n = len(records)
    estimated = counts.get("rate_rule", 0) + counts.get("rate_rule_average", 0)
    return {
        "n_reactions": n,
        "by_source": counts,
        "fraction_rate_rule_estimated": estimated / n,
        "fraction_library_backed": counts.get("library", 0) / n,
    }
