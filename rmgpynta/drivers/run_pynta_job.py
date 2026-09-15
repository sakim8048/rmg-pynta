#!/usr/bin/env python
"""Step 3 (+ corrected step 4): submit a Pynta job, wait for it, fix its libraries.

Runs inside pynta_env: ``<pynta_env_python> run_pynta_job.py <config.json>``

``config.json`` is a plain dict of ``pynta.main.Pynta.__init__`` keyword
arguments (see that class's real signature, quoted in full in this
project's design notes -- ``path``, ``rxns_file``, ``surface_type``,
``metal``, ``label``, ``launchpad_path``, ``repeats``, ``software``, ...),
plus one optional extra key ``_execute_kwargs`` for
``pynta.main.Pynta.execute()``'s own keyword arguments
(``calculate_adsorbates``, ``calculate_transition_states``, ``launch``).

Standalone by design: imports only ``pynta``, stdlib, and
``postprocess_fix`` (same directory) -- never ``rmgpynta`` itself.

IMPORTANT -- shared-launchpad isolation: this project's fireworks get
tagged with ``spec._category = FIREWORKS_CATEGORY`` (see
``_tag_fireworks_with_category`` below) before anything launches them, and
``config["fworker_path"]`` must point to a FWorker file with
``category: 'rmgpynta_test'`` (not the shared ``~/my_fworker.yaml`` other
production runs use, whose ``query: '{}'`` matches literally any firework on
the launchpad). Confirmed the hard way: this project's ``run_pynta_job.py``,
run with ``queue=False`` via ``launch_multiprocess`` and the shared fworker,
pulled and ran fireworks belonging to the user's separate, unrelated DFT
production runs directly on this shared host instead of leaving them for
their own SLURM-based ``qlaunch`` process -- because both launchers were
racing the identical unscoped query against the same launchpad.
``FWorker.category`` (``fireworks/core/fworker.py``) merges ``spec._category``
into its query automatically once set, and ``Pynta.execute(...,
launch=False)`` adds a workflow without launching it (per that method's own
docstring), which is what makes tagging possible before any launcher --
ours or anyone else's category-scoped one -- can claim these fireworks.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from postprocess_fix import write_rmg_libraries_fixed  # noqa: E402

from pynta.main import Pynta  # noqa: E402

FIREWORKS_CATEGORY = "rmgpynta_test"


def _tag_fireworks_with_category(launchpad, label: str, category: str) -> None:
    """Tag every firework in the workflow named ``label`` with
    ``spec._category = category``, so a category-scoped FWorker (see
    ``rmg_pynta_test_fworker.yaml``) only ever picks up this project's own
    fireworks -- see this module's docstring for why that matters on a
    launchpad shared with unrelated production runs. Must run after
    ``Pynta.execute(launch=False)`` (workflow exists on the launchpad) and
    before ``Pynta.launch()`` (nothing has claimed a firework yet).
    """
    wf_doc = launchpad.workflows.find_one({"name": label})
    if not wf_doc:
        raise RuntimeError(
            f"No workflow named {label!r} found on the launchpad after "
            f"execute(launch=False) -- can't tag it for category isolation."
        )
    launchpad.update_spec(wf_doc["nodes"], {"_category": category})


def wait_for_fireworks_completion(launchpad_path, label, poll_interval=60, timeout=None):
    """Poll FireWorks for this job's workflow to finish.

    Pynta submits its workflow named exactly ``self.label`` (confirmed in
    ``pynta/main.py``'s ``execute()``: ``Workflow(self.fws, name=self.label)``),
    and with ``launch=True`` runs an infinite-mode rapidfire launcher (per
    that method's own docstring) that keeps spawning ready fireworks.

    CONFIRMED (previously flagged as unverified): ``LaunchPad.get_wf_summary_dict``
    takes a ``fw_id``, not a ``name`` -- calling it with ``name=label`` raises
    ``TypeError: get_wf_summary_dict() got an unexpected keyword argument 'name'``.
    Its ``'state'`` key does carry the expected ``COMPLETED``/``FIZZLED``
    vocabulary (rmgpy/core/launchpad.py: pulled straight from the workflow
    document's own ``state`` field) -- just needs an fw_id from the named
    workflow first, looked up once since ``self.label`` is stable across a
    workflow's lifetime.
    """
    from fireworks.core.launchpad import LaunchPad

    lp = LaunchPad.from_file(launchpad_path)
    wf_doc = lp.workflows.find_one({"name": label}, {"nodes": 1})
    if not wf_doc:
        raise RuntimeError(f"No workflow named {label!r} found on the launchpad to wait on.")
    fw_id = wf_doc["nodes"][0]

    start = time.time()
    while True:
        summary = lp.get_wf_summary_dict(fw_id)
        state = summary.get("state")
        if state == "COMPLETED":
            return
        if state == "FIZZLED":
            raise RuntimeError(
                f"FireWorks workflow {label!r} FIZZLED -- inspect with "
                f"`lpad get_wflows -n {label}` in pynta_env."
            )
        if timeout is not None and (time.time() - start) > timeout:
            raise TimeoutError(
                f"Timed out after {timeout}s waiting for workflow {label!r} "
                f"(last state: {state!r})"
            )
        time.sleep(poll_interval)


def main(config_path: str) -> None:
    config = json.loads(Path(config_path).read_text())
    execute_kwargs = config.pop("_execute_kwargs", {})
    execute_kwargs.setdefault("calculate_adsorbates", True)
    execute_kwargs.setdefault("calculate_transition_states", True)
    # Intercepted, never passed to execute() itself: launch always happens
    # (if at all) after _tag_fireworks_with_category, never as part of
    # execute() -- see this module's docstring.
    wants_launch = execute_kwargs.pop("launch", True)
    wait_kwargs = config.pop("_wait_kwargs", {})

    job = Pynta(**config)
    job.execute(launch=False, **execute_kwargs)

    if config.get("launchpad_path"):
        _tag_fireworks_with_category(job.launchpad, config["label"], FIREWORKS_CATEGORY)

    if wants_launch:
        job.launch()

    if wants_launch and config.get("launchpad_path"):
        wait_for_fireworks_completion(
            config["launchpad_path"], config["label"], **wait_kwargs
        )
    # If launch=False, the workflow was only added (and tagged) on the
    # launchpad, not run -- this driver returns without waiting, matching
    # Pynta's own execute(launch=False) semantics (queue management left to
    # the caller, e.g. a separate `qlaunch rapidfire` on an HPC scheduler --
    # which must itself be pointed at a category='rmgpynta_test' FWorker to
    # actually pick these up, same as rmg_pynta_test_fworker.yaml).

    n_written = write_rmg_libraries_fixed(
        path=config["path"],
        metal=config["metal"],
        facet=config["surface_type"],
        repeats=config.get("repeats", (3, 3, 4)),
        sites=job.sites,
        site_adjacency=job.site_adjacency,
        slab_path=job.slab_path,
    )
    Path(config["path"], "_rmgpynta_postprocess_result.json").write_text(
        json.dumps({"n_reactions_written": n_written})
    )


if __name__ == "__main__":
    main(sys.argv[1])
