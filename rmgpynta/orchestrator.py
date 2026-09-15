"""The RMGPynta orchestrator class: steps 1 through 8 + data provenance.

See ``rmgpynta/__init__.py`` for why this class never imports ``pynta`` or
``rmgpy`` directly -- every step that needs either one runs as a subprocess
via :mod:`rmgpynta.subprocess_utils` against a standalone script in
:mod:`rmgpynta.drivers`.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from rmgpynta import libraries, provenance, reaction_schema, selection
from rmgpynta.paths import RoundPaths, for_round
from rmgpynta.subprocess_utils import run_in_env

_DRIVERS_DIR = Path(__file__).parent / "drivers"
_RUN_PYNTA_DRIVER = _DRIVERS_DIR / "run_pynta_job.py"
_RUN_RMG_DRIVER = _DRIVERS_DIR / "run_rmg_job.py"
_RUN_UNCERTAINTY_DRIVER = _DRIVERS_DIR / "run_uncertainty_job.py"


@dataclass
class RoundResult:
    round_idx: int
    n_core_reactions: int
    n_edge_reactions: int
    n_selected: int
    library_name: str
    provenance_summary: dict


@dataclass
class RMGPynta:
    """Orchestrates one RMG<->Pynta feedback loop.

    Parameters
    ----------
    work_dir:
        Root directory; one ``round_NNN/`` subdirectory per iteration (see
        :mod:`rmgpynta.paths`).
    pynta_env_python, rmg_env_python:
        Interpreter paths for the two conda environments, e.g.
        ``~/.conda/envs/pynta_env/bin/python`` /
        ``~/.conda/envs/rmg_env/bin/python``.
    pynta_kwargs:
        Keyword arguments forwarded essentially verbatim into
        ``pynta.main.Pynta(...)`` every round (``metal``, ``surface_type``,
        ``label`` is set automatically per round, ``launchpad_path``,
        ``fworker_path``, ``queue_adapter_path``, ``repeats``, ``software``,
        etc.) -- ``path`` and ``rxns_file`` are set automatically per round
        and should not be included here. See the Pynta constructor's real
        (long) signature before filling this in; there is no attempt here
        to re-expose all ~40 of its parameters individually.
    rmg_input_template:
        Path to a base RMG ``input.py`` (``database(...)``,
        ``catalystProperties(...)``, ``species(...)``, ``surfaceReactor(...)``,
        ...) that already runs on its own for round 0. Each later round
        clones it with that round's Pynta-derived libraries spliced into
        ``thermoLibraries``/``reactionLibraries`` (see
        :mod:`rmgpynta.libraries`).
    rmg_database_path:
        Local ``RMG-database`` checkout; round libraries get registered
        under ``<rmg_database_path>/input/{thermo,kinetics}/libraries/``.
    selection_strategy, selection_kwargs:
        See :mod:`rmgpynta.selection`. Default ``"flux_threshold"``.
        ``"sensitivity_uncertainty"`` ranks by RMG-Py's real sensitivity x
        uncertainty contribution instead (see ``selection.py``'s docstring
        and ``run_uncertainty_analysis`` below) -- ``run_loop`` runs that
        driver automatically before round 0's selection when this strategy
        is set; ``selection_kwargs`` should not include ``ranking_path`` for
        that case, since ``run_loop`` fills it in itself.
    max_rounds:
        Upper bound on how many full loops (steps 2-8) to run.
    min_new_reactions:
        Stop early if a round's selection step (7) returns fewer than this
        many reactions -- the design doc's suggested convergence signal.
    allow_unverified_atom_labels:
        Passed through to
        ``rmgpynta.reaction_schema.write_pynta_reactions_yaml``. Defaults
        to False: see ``reaction_schema.AtomLabelsNotVerified`` for why
        this is a hard stop by default rather than a best-effort guess.
    """

    work_dir: Path
    pynta_env_python: Path
    rmg_env_python: Path
    pynta_kwargs: dict
    rmg_input_template: Path
    rmg_database_path: Path
    selection_strategy: "str | Callable" = "flux_threshold"
    selection_kwargs: dict = field(default_factory=dict)
    max_rounds: int = 5
    min_new_reactions: int = 1
    allow_unverified_atom_labels: bool = False

    def __post_init__(self):
        self.work_dir = Path(self.work_dir)

    # ---- steps 2 / 8 ----------------------------------------------------
    def build_reaction_yaml(self, records: list, rp: RoundPaths) -> Path:
        return reaction_schema.write_pynta_reactions_yaml(
            records, rp.reaction_yaml, allow_unverified_labels=self.allow_unverified_atom_labels
        )

    # ---- step 3 (+ corrected step 4) ------------------------------------
    def run_pynta(self, rp: RoundPaths) -> None:
        config = dict(self.pynta_kwargs)
        config["path"] = str(rp.pynta_dir)
        config["rxns_file"] = str(rp.reaction_yaml)
        config.setdefault("label", f"rmgpynta_round_{rp.round_idx:03d}")
        config_path = rp.pynta_dir / "_pynta_config.json"
        config_path.write_text(json.dumps(config, indent=2))
        run_in_env(self.pynta_env_python, _RUN_PYNTA_DRIVER, config_path)

    # ---- step 5 -----------------------------------------------------------
    def register_round_libraries(self, rp: RoundPaths) -> str:
        thermo, kinetics_dir = libraries.collect_pynta_libraries(rp)
        name = rp.library_name()
        libraries.register_libraries(self.rmg_database_path, thermo, kinetics_dir, name)
        libraries.archive_round_libraries(self.work_dir, rp, name)
        return name

    # ---- step 1 (once) / step 6 (every round) ------------------------------
    def run_rmg(self, rp: RoundPaths, extra_thermo_libs: list, extra_kinetics_libs: list) -> None:
        libraries.update_input_py_libraries(
            self.rmg_input_template, rp.input_py, extra_thermo_libs, extra_kinetics_libs
        )
        run_in_env(self.rmg_env_python, _RUN_RMG_DRIVER, rp.input_py, rp.rmg_dir)

    # ---- step 7 (sensitivity_uncertainty strategy only) ---------------------
    def run_uncertainty_analysis(
        self,
        rp: RoundPaths,
        edge: bool = False,
        extra_thermo_libs: list | None = None,
        extra_kinetics_libs: list | None = None,
    ) -> Path:
        """Runs drivers/run_uncertainty_job.py against this round's core (or
        edge, for round >= 1's select_round) chemkin output, producing
        sensitivity_uncertainty_ranking.json for
        selection.select_by_sensitivity_uncertainty to read. Needs
        rmgpy.tools.uncertainty, so it's a subprocess in rmg_env like run_rmg,
        never imported here directly (see rmgpynta/__init__.py).

        ``extra_thermo_libs``/``extra_kinetics_libs`` should be the same
        round-tagged library names ``run_rmg`` was called with for this
        round (``library_names`` in ``run_loop``) -- from round 1 onward the
        model's thermo/kinetics partly comes from Pynta-derived libraries
        registered in earlier rounds, and extract_sources_from_model() needs
        those loaded to attribute sources correctly. species_edge_dictionary.txt
        was confirmed to already include every core species too (not just
        edge-exclusive ones), so it works as a drop-in dictionary here.
        """
        prefix = "chem_edge_annotated" if edge else "chem_annotated"
        chemkin_dir = rp.rmg_dir / "chemkin"
        dictionary_name = "species_edge_dictionary.txt" if edge else "species_dictionary.txt"
        out_dir = rp.rmg_dir / ("uncertainty_edge" if edge else "uncertainty")

        extra_libraries_json = None
        if extra_thermo_libs or extra_kinetics_libs:
            extra_libraries_json = out_dir / "_extra_libraries.json"
            extra_libraries_json.parent.mkdir(parents=True, exist_ok=True)
            extra_libraries_json.write_text(json.dumps({
                "extra_thermo_libraries": extra_thermo_libs or [],
                "extra_reaction_libraries": extra_kinetics_libs or [],
            }))

        args = [
            chemkin_dir / f"{prefix}-gas.inp",
            chemkin_dir / f"{prefix}-surface.inp",
            chemkin_dir / dictionary_name,
            out_dir,
        ]
        if extra_libraries_json is not None:
            args.append(extra_libraries_json)
        run_in_env(self.rmg_env_python, _RUN_UNCERTAINTY_DRIVER, *args)
        return out_dir / "sensitivity_uncertainty_ranking.json"

    # ---- step 7 -------------------------------------------------------------
    def select_round(self, rp: RoundPaths, library_names: list) -> list:
        edge = json.loads(rp.edge_reactions_json.read_text())
        round_selection_kwargs = dict(self.selection_kwargs)
        if self.selection_strategy == "sensitivity_uncertainty":
            round_selection_kwargs["ranking_path"] = self.run_uncertainty_analysis(
                rp, edge=True, extra_thermo_libs=library_names, extra_kinetics_libs=library_names
            )
        selected = selection.select_important_reactions(
            edge, strategy=self.selection_strategy, **round_selection_kwargs
        )
        rp.selected_edge_json.write_text(json.dumps(selected, indent=2))
        return selected

    # ---- data provenance (design doc §03) ------------------------------------
    def log_provenance(self, rp: RoundPaths) -> dict:
        chemkin_dir = rp.rmg_dir / "chemkin"
        records = []
        if chemkin_dir.exists():
            for pattern in ("chem_annotated*.inp", "chem_edge_annotated*.inp"):
                for path in sorted(chemkin_dir.glob(pattern)):
                    records.extend(provenance.parse_annotated_chemkin(path))
        provenance.write_provenance_log(records, rp.provenance_log)
        return provenance.summarize(records)

    # ---- orchestration: step 1 -> (2..8)* -> stop condition --------------------
    def run_loop(self) -> list:
        results: list = []

        # Step 1, once: bootstrap the very first core from whatever
        # libraries the template's own database(...) block already lists.
        rp0 = for_round(self.work_dir, 0)
        self.run_rmg(rp0, extra_thermo_libs=[], extra_kinetics_libs=[])
        core = json.loads(rp0.core_reactions_json.read_text())

        # Pynta's TS-search pipeline is built around adsorbed species and
        # slab sites (see reaction_schema.is_surface_family's docstring) --
        # a homogeneous gas-phase reaction (e.g. H_Abstraction between two
        # gas species) has no surface site for it to build a TS structure
        # from, so those are dropped before build_reaction_yaml ever sees them.
        surface_core = [r for r in core if reaction_schema.is_surface_family(r["family"])]

        # Every later round ranks its edge reactions through select_round()
        # before handing them to Pynta (step 7) -- round 0's core never went
        # through that step, so round 1 used to receive the entire surface
        # core unranked (11,844 reactions in one real run). Apply the same
        # ranking here so round 1 starts from a comparably small, prioritized
        # set instead.
        round0_selection_kwargs = dict(self.selection_kwargs)
        if self.selection_strategy == "sensitivity_uncertainty":
            round0_selection_kwargs["ranking_path"] = self.run_uncertainty_analysis(rp0)
        reactions_for_pynta = selection.select_important_reactions(
            surface_core, strategy=self.selection_strategy, **round0_selection_kwargs
        )
        library_names: list = []

        for round_idx in range(1, self.max_rounds + 1):
            rp = for_round(self.work_dir, round_idx)

            self.build_reaction_yaml(reactions_for_pynta, rp)  # steps 2 / 8
            self.run_pynta(rp)  # step 3 (+ fixed step 4)
            library_names.append(self.register_round_libraries(rp))  # step 5
            self.run_rmg(rp, extra_thermo_libs=library_names, extra_kinetics_libs=library_names)  # step 6

            core = json.loads(rp.core_reactions_json.read_text())
            edge = json.loads(rp.edge_reactions_json.read_text())
            prov_summary = self.log_provenance(rp)  # §03
            selected = self.select_round(rp, library_names)  # step 7

            results.append(
                RoundResult(
                    round_idx=round_idx,
                    n_core_reactions=len(core),
                    n_edge_reactions=len(edge),
                    n_selected=len(selected),
                    library_name=library_names[-1],
                    provenance_summary=prov_summary,
                )
            )

            if len(selected) < self.min_new_reactions:
                break
            reactions_for_pynta = [r for r in selected if reaction_schema.is_surface_family(r["family"])]

        return results
