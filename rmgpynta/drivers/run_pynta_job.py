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
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from postprocess_fix import write_rmg_libraries_fixed  # noqa: E402

from pynta.main import Pynta  # noqa: E402


def wait_for_fireworks_completion(launchpad_path, label, poll_interval=60, timeout=None):
    """Poll FireWorks for this job's workflow to finish.

    Pynta submits its workflow named exactly ``self.label`` (confirmed in
    ``pynta/main.py``'s ``execute()``: ``Workflow(self.fws, name=self.label)``),
    and with ``launch=True`` runs an infinite-mode rapidfire launcher (per
    that method's own docstring) that keeps spawning ready fireworks.

    NOTE: ``LaunchPad.get_wf_summary_dict``'s exact return shape was not
    independently confirmed against the installed ``fireworks`` version on
    this machine while writing this module -- verify the ``'state'`` key
    and the ``COMPLETED``/``FIZZLED`` spelling before relying on this in
    production, e.g. inside pynta_env:
    ``python -c "from fireworks import LaunchPad; help(LaunchPad.get_wf_summary_dict)"``
    """
    from fireworks.core.launchpad import LaunchPad

    lp = LaunchPad.from_file(launchpad_path)
    start = time.time()
    while True:
        summary = lp.get_wf_summary_dict(name=label)
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
    execute_kwargs.setdefault("launch", True)
    wait_kwargs = config.pop("_wait_kwargs", {})

    job = Pynta(**config)
    job.execute(**execute_kwargs)

    if execute_kwargs.get("launch") and config.get("launchpad_path"):
        wait_for_fireworks_completion(
            config["launchpad_path"], config["label"], **wait_kwargs
        )
    # If launch=False, the workflow was only added to the launchpad, not
    # run -- this driver returns without waiting, matching Pynta's own
    # execute(launch=False) semantics (queue management left to the caller,
    # e.g. a separate `qlaunch rapidfire` on an HPC scheduler).

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
