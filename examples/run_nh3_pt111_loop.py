#!/usr/bin/env python
"""Example wiring for the NH3/Pt(111) system this project was designed around.

This is illustrative, not a script meant to run unmodified -- fill in the
placeholder paths (marked below) for your own machine before running it.
Run with whatever Python has ``rmgpynta`` and ``PyYAML`` installed; it does
NOT need ``pynta`` or ``rmgpy`` importable itself (see rmgpynta/__init__.py).

CURRENT STATUS: only round 0 (the RMG-only bootstrap) can actually
complete today. Round 1 will raise ``AtomLabelsNotVerified`` at
``build_reaction_yaml`` -- drivers/run_rmg_job.py doesn't yet recover
verified *1..*4 reacting-atom labels from RMG Reaction objects, and
RMGPynta refuses to build a reaction.yaml without them by default (see
README.md's "Known gaps"). Likely to also hit UnsupportedReactionFamily
on some core reactions, since only 3 of the ~13 kinetics families in
rmg_ch4_pt111/input.py are mapped. Nothing here has been run end-to-end
yet -- this whole file has never been executed even once.
"""
from pathlib import Path

from rmgpynta import RMGPynta

HOME = Path.home()

orchestrator = RMGPynta(
    work_dir=HOME / "rmg-pynta-runs" / "nh3_pt111",
    # the real envs on this machine (~/miniconda3/envs/{pynta_env,rmg_env})
    pynta_env_python=HOME / "miniconda3" / "envs" / "pynta_env" / "bin" / "python",
    rmg_env_python=HOME / "miniconda3" / "envs" / "rmg_env" / "bin" / "python",
    # --- forwarded ~verbatim into pynta.main.Pynta(...) every round ---
    # See that class's real __init__ signature before changing these; only
    # `path`, `rxns_file`, and `label` are set automatically per round.
    pynta_kwargs=dict(
        metal="Pt",
        surface_type="fcc111",
        repeats=(3, 3, 4),
        software="MACECalculator",  # or "VASP" -- see the Pynta-VASP testing conversation
        # real files on this machine, but dated 2022-2023 -- confirm they
        # still point at a live MongoDB/queue allocation before using them;
        # ~/pynta-vasp/test/ has a more recently-touched my_launchpad.yaml +
        # FW_config.yaml if that testing's allocation is the one to use
        # instead (it has no matching fworker/qadapter files of its own).
        launchpad_path=str(HOME / "my_launchpad.yaml"),
        fworker_path=str(HOME / "my_fworker.yaml"),
        queue_adapter_path=str(HOME / "my_qadapter.yaml"),
    ),
    # --- a real, working RMG input.py for round 0, e.g. a copy of ---
    # ~/rmg-production/rmg_ch4_pt111/input.py
    rmg_input_template=HOME / "rmg-production" / "rmg_ch4_pt111" / "input.py",
    rmg_database_path=HOME / "RMG-Py" / "RMG-database",
    selection_strategy="flux_threshold",
    selection_kwargs=dict(temperature_k=900.0, top_n=15),
    max_rounds=5,
    min_new_reactions=1,
)

if __name__ == "__main__":
    results = orchestrator.run_loop()
    for r in results:
        print(
            f"round {r.round_idx}: core={r.n_core_reactions} edge={r.n_edge_reactions} "
            f"selected={r.n_selected} library={r.library_name!r} "
            f"rate_rule_fraction={r.provenance_summary.get('fraction_rate_rule_estimated')}"
        )
