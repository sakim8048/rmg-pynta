# Bug Fix Log

Records bugs found and fixed across the RMG-Pynta project, one entry per day
a fix lands. Covers both repos:

- `rmg-pynta` (this repo, branch `RMG-Pynta`, created 2026-09-08)
- `pynta` (`/home/shikim/pynta`, branch `update_environment_preprocessing`),
  starting from 2026-09-08 (when `rmg-pynta` was created) onward

Newest first.

## 2026-09-21

### rmg-pynta

Commit [`83cb7e0`](https://github.com/sakim8048/rmg-pynta/commit/83cb7e0) — found via a real round-1 RMG rerun using a Pynta-registered thermo library.

- **`postprocess_fix.py`, `write_rmg_libraries_fixed`** — wrote Pynta's raw NASA fit text and facet string straight through into the generated thermo library, breaking downstream RMG loading two ways. First, Pynta fits with `Tmin=298.15 K`, but `rmgpy.thermo.nasa.NASA.to_thermo_data()` (called from `ThermoDatabase.correct_binding_energy` while loading the library) looks up H298/S298 at exactly 298 K — outside the segment's own bounds — raising `ValueError: No valid NASA polynomial at temperature 298 K`. Fixed by rewriting `Tmin=(298.15,'K')` to `Tmin=(298.0,'K')` in the generated species text, matching Pynta's own placeholder/zero-Cp species convention. Second, Pynta's `surface_type` (e.g. `"fcc111"`) includes the crystal-structure prefix, but RMG's surface binding-energy database (`RMG-database/input/surface/libraries/metal.py`) keys facets by Miller index only (e.g. `"Pt111"`); `ThermoDatabase.get_thermo_data`'s scaling lookup builds `db_label = entry.metal + entry.facet`, so the raw value produced `"Ptfcc111"` and `DatabaseError: Metal 'Ptfcc' not found in database`. Fixed by stripping the known crystal-structure prefixes (`fcc`, `bcc`, `hcp`, `sc`, `diamond`, `rocksalt`, `hexagonal`) before writing the library header's facet field.

## 2026-09-17

### rmg-pynta

Commit [`0f0b242`](https://github.com/sakim8048/rmg-pynta/commit/0f0b242) — found after `round_001`'s post-launch step sat blocked for two days straight in the `pynta-rmg-test` project.

- **`run_pynta_job.py`** — `Pynta.launch()` (`pynta/main.py`) runs FireWorks' `rapidfirequeue(..., nlaunches="infinite")`, which polls and resubmits forever and never returns on its own. The driver called it inline, ahead of `wait_for_fireworks_completion()`/`write_rmg_libraries_fixed()`, so once a round's actual FireWorks work reached state `COMPLETED` there was nothing left to launch but the process kept blocking anyway — confirmed on `round_001`: its workflow (536/536 fireworks) had already reached `COMPLETED`, while the process just kept logging "N jobs in queue / sleeping 60 secs" forever, never reaching the library-update step. `round_001` itself had to be recovered by hand (kill the stuck process, rerun just the wait+postprocess tail against the already-completed workflow). Fixed going forward by running `launch()` in its own subprocess (`--launch-only` mode, same file) while the parent polls completion independently and kills the launcher once the workflow finishes, then runs the library update automatically.

## 2026-09-15

### pynta

Commit [`593d2c52`](https://github.com/zadorlab/pynta/commit/593d2c52). Found while debugging fizzled `MolecularTSEstimate` fireworks (17143–17155) in the `pynta-rmg-test` project.

- **`pynta/utils.py`, `get_occupied_sites`** — raised a bare `ValueError` whenever an adsorbate atom had no candidate site within reach (e.g. `sites` came back empty for a gas-phase-only TS branch). Now skips that atom instead, matching how an out-of-cutoff nearest site is already skipped just below — "no candidate site" isn't an error, just nothing to mark occupied.
- **`pynta/geometricanalysis.py`, `generate_adsorbate_molecule`** — when computing `neighbor_sites` with `max_dist` set (every caller in the codebase passes `max_dist=np.inf` to mean "no distance filtering"), the code filtered candidate sites by proximity to `target_sites` (currently-occupied sites). If nothing was occupied yet (e.g. `get_unique_TS_structs`'s gas-phase-reactant branch, which builds the 2D graph on a bare slab *before* placing the adsorbate), `target_sites` was empty, so the loop appended nothing — `neighbor_sites` silently came back `[]` even under `max_dist=np.inf`. This broke an invariant already documented and relied on elsewhere (`pynta/utils.py`'s `_interaction_terms`: "admol was built with max_dist=inf, where generate_adsorbate_molecule keeps ninds=range(len(sites))"). Fixed by falling back to the full site list when there's no occupied site to anchor the distance filter.
- **`pynta/calculator.py`, `map_harmonically_forced`** — `os.makedirs(os.path.join(path,str(j)))` had no `exist_ok=True`, so any rerun of a TS estimate that reached the harmonic-mapping stage a second time crashed with `FileExistsError` on the leftover numbered directory from the first attempt (only two files, `harm.json`/`xtb.xyz`, are ever written there, both opened in overwrite mode, so reuse is safe). Fixed by adding `exist_ok=True`.
- **`pynta/calculator.py`, `get_energy_forces`** (4 occurrences, all copies of the same nested calculator class) — `energy = 0.0` accumulates `energy += E` across `atom_bond_potentials`/`site_bond_potentials`, then returns `energy[0][0]`. When both lists were empty (nothing to enforce for a candidate — e.g. a still-effectively-gas-phase species), `energy` stayed a bare Python `float` and `energy[0][0]` raised `TypeError: 'float' object is not subscriptable`. Fixed by returning `float(np.asarray(energy).reshape(-1)[0])`, which works whether `energy` ended up scalar or array-shaped.
- **`pynta/calculator.py`, `run_harmonically_forced`/`run_harmonically_forced_no_pbc`** — the `try/except: return None,None,None` around `opt.run(...)` swallowed the real exception silently, making every downstream failure ("no harmonically-mapped TS guesses survived") undiagnosable without re-running by hand. Added `traceback.print_exc()` before the fallback return so the actual cause shows up in the job log.
- **`pynta/tasks.py`, `MolecularTSEstimate.run_task`** — `Emin = np.min(np.array(Es))` crashed with a cryptic `ValueError: zero-size array to reduction operation minimum which has no identity` whenever every TS guess failed to produce a harmonic mapping (`Es` ends up empty). Added an explicit check that raises `ValueError("No harmonically-mapped TS guesses survived filtering for reaction: {rxn_name}")` instead.
- **`pynta/transitionstate.py`, `get_unique_TS_structs`** (line 182) — the actual root cause behind firework 17155 (`X + CH4 => CH4X`) reaching zero harmonic-mapped guesses. In the gas-phase-reactant branch (`num_surf_sites[0] == 0`), `add_adsorbate_to_site(adslab, adsorbate=adslab, ...)` passed the bare slab copy as *both* the base structure and the "adsorbate" being added onto it — the same Python object for both arguments. `add_adsorbate_to_site` aliases (not copies) its `adsorbate` argument, translates/rotates it in place, then does `atoms += ads`; since `atoms` and `ads` were the same mutated object, this duplicated every atom at coincident positions (zero interatomic distance), which is exactly why the MACE potential reported `Energy=inf`/`fmax=nan` at optimization step 0 and Sella failed immediately. Fixed by passing the actual isolated gas-phase molecule (`adss[0]`, the correct per-species structure already used the same way two lines later in this function) instead of `adslab`. Confirmed fixed by direct repro: `MolecularTSEstimate.run_task` on firework 17155's spec now converges (`Energy=-232.4 eV`, `fmax=0.0071`) and returns `SUCCESS`; firework 17155 subsequently `COMPLETED` through FireWorks.

### rmg-pynta

Commit [`5916a5b`](https://github.com/sakim8048/rmg-pynta/commit/5916a5b) — surfaced while getting a live run from RMG's round-0 core through to a real Pynta TS search.

- **`reaction_schema.py`** — `RMG_TO_PYNTA_FAMILY` was missing 5 of RMG's surface kinetics families (`Surface_Adsorption_vdW`, `_Single`, `_Double`, `_Bidentate`, `Surface_Abstraction_vdW`), so reactions in those families would hit `UnsupportedReactionFamily` instead of reaching Pynta. Fixed by mapping all 5 in, after verifying only the literal string `"Surface_Migration"` actually changes behavior anywhere in `pynta/*.py` (so the mapping is safe).
- **`orchestrator.py`** — gas-phase-only reactions (no surface site) were being passed straight to Pynta, which has nothing to build a TS from without one; now filtered out via `is_surface_family` before dispatch. Separately, round 0's core reactions were never being ranked before being handed to Pynta — one real run sent all 11,844 reactions instead of a prioritized handful. Now applies `select_important_reactions` to round 0 the same as later rounds.
- **`run_rmg_job.py`** — reacting-atom labels (`*1`..`*4`) weren't being recovered from RMG `Reaction` objects, causing `AtomLabelsNotVerified`. Fixed via `KineticsFamily.add_atom_labels_for_reaction`, called on deep copies (the live call mutates shared `Species` objects otherwise). Also, each side's combined adjacency block was built by concatenating per-species text, leaving stray multiplicity lines and non-continuous atom numbering; switched to `Molecule.merge()`. The block's multiplicity was hardcoded to `1`, which RMG's own parser rejects for radical-bearing blocks under Hund's rule — now computed from the real structure.
- **`run_pynta_job.py`** — fireworks created by this project weren't tagged with `spec._category` before anything could launch them. Confirmed concretely: with `queue=False`, the shared unscoped `FWorker`'s local `launch_multiprocess` pulled and ran a *separate production DFT workflow's* fireworks directly on this host instead of leaving them for their own SLURM allocation. Fixed by reordering to `execute(launch=False)` → tag → `launch()`, so a category-scoped `FWorker` only ever claims this project's own fireworks. Also fixed `wait_for_fireworks_completion`, which called `LaunchPad.get_wf_summary_dict` with a `name=` kwarg it doesn't accept — needs an `fw_id` instead.
- **`selection.py`, `run_uncertainty_job.py`** (new) — added a sensitivity x uncertainty selection strategy (via `rmgpy.tools.uncertainty.Uncertainty`, ranked against CH4/NH3 observables) and wired it into `orchestrator.select_round()` for round 1+, passing each round's accumulated library names through so Pynta-derived libraries get attributed correctly.

## 2026-09-09

### pynta

Commit [`cc7ae7c9`](https://github.com/zadorlab/pynta/commit/cc7ae7c9)

- **`pynta/geometricanalysis.py`** — `pynta/postprocessing.py`'s `get_TS` has called `validate_diffusion_TS` for `Surface_Migration` TSs since an earlier merge, but the function itself — and its helpers `_diffusion_endpoint_2D` and `get_diffusion_connecting_sites` — only ever existed on the `covdep_mc` branch's `geometricanalysis.py`, never ported over. Every run through `get_TS`/`get_kinetics`/`postprocess` crashed with `NameError: name 'validate_diffusion_TS' is not defined`. Ported verbatim from `origin/covdep_mc:pynta/geometricanalysis.py`; no changes needed since `generate_adsorbate_2D`'s signature already matched.

## 2026-09-08

### pynta

Commit [`f58de648`](https://github.com/zadorlab/pynta/commit/f58de648)

- **`pynta/tasks.py`, `MolecularTSEstimate`** — counted a reactant surface site as "empty" only when it had zero bonds outright; fixed to count it as empty when every bonded neighbor is itself a surface site.
- **`pynta/transitionstate.py`, `get_unique_TS_templates_site_pairings`** — raised when a required reaction-bond label was missing from `tsmol` (happens when a multidentate template merges onto a monodentate/vdW-matched geometry and drops a label); now skips that TS candidate instead of raising.

### rmg-pynta

Commit [`90d1a5c`](https://github.com/sakim8048/rmg-pynta/commit/90d1a5c) — `examples/run_nh3_pt111_loop.py`

- `launchpad_path`/`fworker_path`/`queue_adapter_path` pointed at `~/fw_config/`, which doesn't exist on this machine. Repointed at the real top-level configs (`~/my_launchpad.yaml` etc.), with a note that those are dated 2022-2023 and `~/pynta-vasp/test/` has a more recently touched (but incomplete) alternative set worth checking first. Also corrected the module docstring, which implied more of the pipeline worked than actually did — only round 0 (RMG bootstrap) could complete at the time; round 1 hit `AtomLabelsNotVerified` by design (fixed 2026-09-15, above).
