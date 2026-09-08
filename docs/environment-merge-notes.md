# Could pynta_env and rmg_env be merged?

Investigated 2026-09-08, prompted by wanting to skip the JSON hand-off
between `drivers/run_pynta_job.py` and `drivers/run_rmg_job.py` by running
both in one process. **Conclusion: not right now** — one real, unresolved
package collision, plus a smaller channel-policy conflict. Findings below
so this doesn't need re-investigating from scratch later.

## Python version: not actually a blocker

Current deployed versions differ (`pynta_env` 3.9.25, `rmg_env` 3.11.15),
but that's a solver-speed choice, not a hard constraint:

- `~/pynta/environment.yml`: `python =3.9` (pinned exactly)
- `~/RMG-Py/RMG-Py/CLAUDE.md`: *"Python is pinned >=3.9,<3.12"*, and
  `~/RMG-Py/RMG-Py/environment.yml` pins `python =3.11` with an inline
  comment: *"pinned to a single version to keep solver fast; CI rewrites
  this to the matrix version"*.

So **Python 3.9 is a real, shared intersection** — RMG-Py's CI already
tests against it. A merged env pinned to 3.9 is plausible from the
interpreter side alone.

## The actual blocker: `pysidt` vs `pysidt-rmg`

Both projects depend on a package that imports as `pysidt` -- but they
pull it from two different, separately-versioned distributions that
happen to share the same import namespace:

```
$ conda list -p ~/miniconda3/envs/pynta_env | grep pysidt
pysidt        1.1.0   py39h4e4d2e0_11   mjohnson541

$ conda list -p ~/miniconda3/envs/rmg_env | grep pysidt
pysidt-rmg    1.2.0   py_83             rmg
```

Both land at the identical import path in their respective envs:

```
pynta_env:  .../envs/pynta_env/lib/python3.9/site-packages/pysidt/__init__.py
rmg_env:    .../envs/rmg_env/lib/python3.11/site-packages/pysidt/__init__.py
```

This isn't a version-range conflict solvable by relaxing a pin -- they're
two separately-maintained package *distributions* (different conda
channels: `mjohnson541` vs `rmg`) that both install to `site-packages/pysidt/`.
Installing both in one environment means one silently overwrites the
other's files on disk; whichever "loses" leaves its dependent tool
running against the wrong version with no error at install time.

Related, likely-connected package: `~/pynta/environment.yml` also lists
`rmgmolecule` (0.3.0, `mjohnson541` channel) directly -- RMG-Py has no
equivalent external dependency here since `rmgpy.molecule` is vendored in
the main `rmgpy` package instead. Not a collision by itself, but evidence
the two projects have independently forked/renamed what was probably
originally shared code, which is the same pattern behind the `pysidt`
split.

## Secondary issue: contradictory channel policy

```
pynta_env  (~/pynta/environment.yml):        channels: [defaults, mjohnson541, conda-forge]
rmg_env    (~/RMG-Py/RMG-Py/environment.yml): channels: [conda-forge, rmg, nodefaults]
```

RMG-Py explicitly excludes the `defaults` channel (`nodefaults`); pynta
explicitly includes it first. Even where individual package versions
happen to line up, mixing these channel lists will make the solver's job
harder and can pull in inconsistent builds across the two ecosystems.

## What's confirmed *not* a problem

Direct version check (both envs, 2026-09-08):

| package | pynta_env | rmg_env | pynta's own constraint | RMG-Py's own constraint |
|---|---|---|---|---|
| numpy | 1.26.4 | 1.26.4 | unconstrained | `>=1.24,<2` |
| scipy | 1.13.1 | 1.17.1 | `>=1.1.0` | `>=1.13` |
| rdkit | 2025.03.5 | 2025.03.6 | not pinned directly (transitive) | `>=2024` |

These overlap fine. `rmg_env` has no `ase` at all (RMG-Py doesn't need
it), which is a non-issue for a merge, not a conflict.

## If this gets revisited

The concrete, checkable first step isn't "attempt a full merge" -- it's a
narrow spike: **does pynta actually work against `pysidt-rmg 1.2.0`
instead of its own pinned `pysidt 1.1.0`?** If RMG's fork is a compatible
superset (plausible, given the version bump and shared origin), pynta's
`environment.yml` could point at `rmg::pysidt-rmg` instead, removing the
collision without needing RMG-Py to change anything. That's a pynta-side
experiment (swap the dependency, run pynta's own test suite) rather than
a new merged-environment build.

## Current decision

Keep `pynta_env` / `rmg_env` separate; `rmgpynta`'s subprocess-boundary
design (see main `README.md`) stays as-is. This note exists so a future
attempt starts from "swap pynta's pysidt pin and test" instead of
re-discovering the collision from scratch.
