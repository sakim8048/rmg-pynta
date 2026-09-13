# rmgpynta

Orchestrator for a closed feedback loop between [RMG](https://github.com/ReactionMechanismGenerator/RMG-Py)
(Reaction Mechanism Generator) and [Pynta](https://github.com/zadorlab/pynta)
(automated DFT transition-state search): RMG proposes the reactions worth
knowing, Pynta computes them at DFT accuracy, and the result feeds back into
RMG's own thermo and kinetics libraries, round after round.

Design rationale, the step-by-step method, and the data-provenance design
live in the **RMG–Pynta Feedback Loop** artifact/Notion page this repo
implements; this README covers the code itself.

## Why this exists

Pynta already writes RMG-format libraries automatically after every run
(`pynta.tasks.PostprocessingTask` → `pynta.postprocessing.write_rmg_libraries`,
confirmed present in the run directories under `~/pynta-production/`).
Nothing reads those libraries back into an RMG `input.py` — the libraries
get generated and sit there. This package closes that loop.

## Package layout

```
rmgpynta/
  orchestrator.py     RMGPynta class -- run_loop() ties everything together
  reaction_schema.py  steps 2/8: RMG reaction records -> pynta reaction.yaml
  libraries.py        steps 4/5: collect + register Pynta's RMG libraries
  selection.py         step 7: pick which edge reactions to compute next
  provenance.py        parse RMG's annotated chemkin comments -> a per-round log
  paths.py              round-indexed work_dir layout
  subprocess_utils.py   run a driver script inside pynta_env / rmg_env
  drivers/
    run_pynta_job.py      runs INSIDE pynta_env -- imports `pynta`
    run_rmg_job.py         runs INSIDE rmg_env -- imports `rmgpy`
    postprocess_fix.py     corrected write_rmg_libraries (see below)
```

`rmgpynta` itself never imports `pynta` or `rmgpy` — they live in separate
conda environments on this machine (`pynta_env`, `rmg_env`) that don't
import each other either. Every step that needs one of them shells out to
that environment's own interpreter running one of the `drivers/` scripts,
which are standalone (stdlib + `pynta` or `rmgpy` only).

%## Install

%```bash
%pip install --cert ~/snl-ca.pem -e ".[test]"  # --cert needed on this machine: corporate TLS-interception proxy
%pytest
%```

%Tested against Python 3.13 (test run) and syntax-checked against 3.8 (the
%system `/usr/bin/python3`, and this repo's stated `requires-python`).

%This installs `rmgpynta` for whichever Python runs the orchestrator loop
%itself (needs only `PyYAML`). Point `pynta_env_python` / `rmg_env_python`
%at your existing conda envs -- `rmgpynta` does not need to be installed
%inside either of them.

## A real bug this project works around

`pynta.postprocessing.write_rmg_libraries` has a confirmed bug: the loop
that builds `reaction_library/reactions.py` is guarded by
`if reaction_text == "":` with the body-append *inside* that guard, so only
the *first* transition state's kinetics ever gets written, no matter how
many TS the run actually computed. `rmgpynta.drivers.postprocess_fix`
reimplements the same logic with that guard fixed, and every Pynta round run
through this package uses the corrected version instead of trusting Pynta's
own auto-written file. This fix should also be proposed to `zadorlab/pynta`
directly, independent of this project.

## A real RMG-Py constraint this project designs around

RMG-Py's `surface_reactor(...)` (what `surfaceReactor(...)` in an
`input.py` maps to) raises `NotImplementedError` if given `sensitivity=[...]`
(confirmed in `rmgpy/rmg/input.py`). Since this system is inherently a
surface-chemistry reactor, RMG's built-in sensitivity analysis cannot be
used for step 7 -- `rmgpynta.selection`'s default (`"flux_threshold"`) ranks
edge reactions by an approximate rate constant at the reactor temperature
instead, and `strategy="sensitivity"` raises immediately with an
explanation rather than silently doing the wrong thing.

## Known gaps -- read before running against a real system

These are left as explicit, documented failures rather than guessed-at
implementations, because guessing wrong here means either wasting real HPC
compute on a malformed Pynta job or silently corrupting a mechanism:

- **`reaction_schema.synthesize_diffusion_entries` is not implemented.**
  RMG's own kinetics families don't include site-to-site diffusion, but
  Pynta's real `reaction.yaml` files do carry `Diffusion`-family entries
  (e.g. `hb-yesdiffusion.yaml`). No verified real example of a diffusion
  entry's exact `*1`/`*2` atom-labeling was available while writing this
  module. Before implementing: open one real `Diffusion`-family entry under
  `~/pynta-production/*/hb-yesdiffusion.yaml` and hand-verify the labeling
  against a working Pynta TS-search run.
- **Reacting-atom (`*1`..`*4`) labels are not recovered from RMG `Reaction`
  objects.** `drivers/run_rmg_job.py` dumps every reaction record with
  `"atom_labels_verified": false`, and
  `reaction_schema.write_pynta_reactions_yaml` refuses to build a
  `reaction.yaml` from an unverified record unless you pass
  `allow_unverified_labels=True` (see `AtomLabelsNotVerified`'s docstring
  for what's actually needed to fix this properly -- it requires the
  generating `KineticsFamily`'s reaction "recipe"/action list, not just the
  instantiated product/reactant `Species`).
- **`RMG_TO_PYNTA_FAMILY` only maps 3 of the ~13 kinetics families** listed
  in a real `input.py` (`Surface_Dissociation`, `Surface_Adsorption_Dissociative`,
  `Surface_Abstraction`). The rest raise `UnsupportedReactionFamily` rather
  than guess a mapping -- see `reaction_schema.py`'s
  `UNMAPPED_RMG_FAMILIES_SEEN_IN_INPUT_PY`.
- **FireWorks completion polling (`run_pynta_job.wait_for_fireworks_completion`)**
  uses `LaunchPad.get_wf_summary_dict(...)['state']`, which was not
  independently confirmed against the installed `fireworks` version on this
  machine. Verify before relying on it for a long-running HPC job.
- **RMG library priority order** (does an earlier-listed
  `thermoLibraries`/`reactionLibraries` entry win over a later one?) is not
  confirmed. `libraries.update_input_py_libraries` inserts new, more
  DFT-informed libraries at the *front* of each list on the assumption that
  matters, but this should be verified against RMG-Py's actual database
  loading order.
- **Edge-reaction "importance" is approximated by rate constant at the
  reactor temperature**, not RMG's true internal leak-flux/edge-flux
  estimate (a confirmed, stable API for the latter was not available while
  writing this). Reasonable as a first cut; revisit if round-to-round
  selection looks off.
- **`pynta_env` and `rmg_env` are kept separate on purpose.** A merged
  environment was investigated and hits a real package collision
  (`pysidt` vs `pysidt-rmg`, two different distributions installing to the
  same import path) -- see [`docs/environment-merge-notes.md`](docs/environment-merge-notes.md)
  before attempting it again.

## Round directory layout

```
work_dir/
  round_000/                    # step 1, bootstrap: RMG only, no pynta/ yet
    rmg/    input.py, chemkin/, seed/, RMG.log, core_reactions.json, edge_reactions.json
  round_001/                    # first full loop of steps 2-8
    pynta/  reaction.yaml (built from round_000's core), thermo_library.py, reaction_library/
    rmg/    input.py (now lists round_001's registered libraries), core_reactions.json, edge_reactions.json, selected_edge.json
    provenance/  reaction_provenance_round001.json
  round_002/  reaction.yaml built from round_001's *selected edge*, ...
  libraries/  round-tagged archival copies of every round's thermo/kinetics libraries
```
