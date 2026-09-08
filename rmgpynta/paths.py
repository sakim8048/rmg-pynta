"""Round-indexed directory layout shared by every module in this package.

Layout (see README.md for the full picture)::

    work_dir/
      round_000/
        rmg/          input.py, chemkin/, seed/, RMG.log, edge_reactions.json, core_reactions.json
        pynta/        reaction.yaml, slab.xyz, sites.xyz, Adsorbates/, TS0.../
          thermo_library.py         # written by pynta's own PostprocessingTask
          reaction_library/         # reactions.py, dictionary.txt -- ditto (bug-fixed, see libraries.py)
        provenance/   reaction_provenance_round000.json
      round_001/ ...
      libraries/                    # round-tagged copies registered into RMG-database's search path
        thermo/round_000_thermo_library.py, ...
        kinetics/round_000_reaction_library/, ...
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RoundPaths:
    """All paths for one iteration of the RMG<->Pynta loop."""

    round_idx: int
    root: Path

    @property
    def round_dir(self) -> Path:
        return self.root / f"round_{self.round_idx:03d}"

    @property
    def rmg_dir(self) -> Path:
        return self.round_dir / "rmg"

    @property
    def pynta_dir(self) -> Path:
        return self.round_dir / "pynta"

    @property
    def provenance_dir(self) -> Path:
        return self.round_dir / "provenance"

    @property
    def input_py(self) -> Path:
        return self.rmg_dir / "input.py"

    @property
    def reaction_yaml(self) -> Path:
        return self.pynta_dir / "reaction.yaml"

    @property
    def core_reactions_json(self) -> Path:
        """Written by drivers/run_rmg_job.py after RMG.execute() completes."""
        return self.rmg_dir / "core_reactions.json"

    @property
    def edge_reactions_json(self) -> Path:
        """Written by drivers/run_rmg_job.py after RMG.execute() completes."""
        return self.rmg_dir / "edge_reactions.json"

    @property
    def selected_edge_json(self) -> Path:
        """Written by rmgpynta.selection after filtering the edge."""
        return self.rmg_dir / "selected_edge.json"

    @property
    def provenance_log(self) -> Path:
        return self.provenance_dir / f"reaction_provenance_round{self.round_idx:03d}.json"

    @property
    def pynta_thermo_library(self) -> Path:
        """Where Pynta's PostprocessingTask writes it (pynta/postprocessing.py:write_rmg_libraries)."""
        return self.pynta_dir / "thermo_library.py"

    @property
    def pynta_reaction_library_dir(self) -> Path:
        return self.pynta_dir / "reaction_library"

    def library_name(self) -> str:
        """Name used when registering this round's libraries with RMG-database."""
        return f"round_{self.round_idx:03d}"

    def ensure(self) -> "RoundPaths":
        for d in (self.rmg_dir, self.pynta_dir, self.provenance_dir):
            d.mkdir(parents=True, exist_ok=True)
        return self


def for_round(work_dir: Path, round_idx: int) -> RoundPaths:
    return RoundPaths(round_idx=round_idx, root=Path(work_dir)).ensure()


def libraries_root(work_dir: Path) -> Path:
    """Shared, cumulative, round-tagged library copies (rmgpynta.libraries.register_libraries)."""
    root = Path(work_dir) / "libraries"
    (root / "thermo").mkdir(parents=True, exist_ok=True)
    (root / "kinetics").mkdir(parents=True, exist_ok=True)
    return root
