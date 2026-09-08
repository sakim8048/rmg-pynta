"""MOCK / EXPLANATORY FILE -- not meant to be run for real (see
run_nh3_pt111_loop.py for that). This exists to answer one question
concretely: "where does Pynta actually get imported?"

Answer: only in ZONE 3 below, and even then indirectly -- through
subprocess calls into driver scripts, not a top-of-file `import pynta`.
Neither ZONE 1 nor ZONE 2 imports anything, and that's not a gap this
project should fill -- it's how both file formats are designed to work.

    ZONE 1: RMG's input.py   -- a sandboxed DSL, exec()'d by RMG itself.
                                 Cannot import anything (see below).
    ZONE 2: Pynta's reaction.yaml -- plain data. YAML files don't import.
    ZONE 3: THIS file        -- the only place a human runs `python ...`
                                 directly, and the only place `import`
                                 statements for orchestration exist.
"""

# ============================================================================
# ZONE 1 -- RMG's input.py (a real one lives at
#           ~/rmg-production/rmg_ch4_pt111/input.py)
# ============================================================================
#
# This is NOT executed as a normal Python module. RMG's own loader
# (rmgpy/rmg/input.py:read_input_file) does this instead:
#
#     global_context = {'__builtins__': None}          # no builtins at all
#     local_context = {'database': ..., 'species': ..., 'surfaceReactor': ...}
#     exec(f.read(), global_context, local_context)
#
# `database(...)`, `species(...)`, `surfaceReactor(...)` etc. are names RMG
# injects into that exec() namespace before running the file's text -- they
# are NOT imported by the file itself. With `__builtins__` stripped, the
# file structurally cannot `import` anything even if you added the line.
# It is a configuration DSL that happens to use Python's call syntax, not
# a program. Mock contents (this is what's really inside a real one):
#
#     database(
#         thermoLibraries=['surfaceThermoPt111', 'primaryThermoLibrary'],
#         reactionLibraries=[],                 # <- rmgpynta fills this in
#         kineticsFamilies=['Surface_Dissociation', 'Surface_Adsorption_Dissociative', ...],
#     )
#     catalystProperties(metal='Pt111')
#     species(label='NH3', reactive=True, structure=SMILES('N'))
#     species(label='X', reactive=True, structure=adjacencyList("1 X u0"))
#     surfaceReactor(temperature=(900, 'K'), initialPressure=(1.0, 'bar'), ...)
#
# rmgpynta never edits this file's Python syntax -- it only text-splices
# new library *names* into the `thermoLibraries=[...]` / `reactionLibraries=[...]`
# lists between rounds (see rmgpynta.libraries.update_input_py_libraries).
# Nothing here ever mentions Pynta.


# ============================================================================
# ZONE 2 -- Pynta's reaction.yaml (auto-generated, never hand-written for
#           this workflow -- doesn't exist until round 1 runs)
# ============================================================================
#
# Also not Python -- plain YAML, parsed with `yaml.safe_load`. No imports
# are possible in YAML at all. Mock contents (one entry, real schema):
#
#     - index: 0
#       reaction: "NX + X => NHX + HX"
#       reaction_family: "Dissociation"
#       reactant: |
#         multiplicity 1
#         1 N u0 p1 c0 {2,D} {3,S}
#         2 X u0 p0 c0 {1,D}
#         3 H u0 p0 c0 {1,S}
#       product: |
#         multiplicity 1
#         1 N u0 p1 c0 {2,T}
#         2 X u0 p0 c0 {1,T}
#
# rmgpynta.reaction_schema.write_pynta_reactions_yaml() is what generates
# this, from ZONE 1's RMG output. Nothing here mentions RMG either -- by
# the time this file exists, RMG's job (for this round) is already done.


# ============================================================================
# ZONE 3 -- the actual "universal input": THIS file. The only one you run
#           directly (`python universal_input_explained.py`), and the only
#           place any `import` statement for orchestration appears.
# ============================================================================
from pathlib import Path

from rmgpynta import RMGPynta          # <-- the "import" you were looking for

HOME = Path.home()

orchestrator = RMGPynta(
    work_dir=HOME / "rmg-pynta-runs" / "mock",

    # --- these two paths are how THIS file reaches "into" Pynta and RMG ---
    # NOT by importing pynta/rmgpy here -- rmgpynta never does that (see
    # rmgpynta/__init__.py's docstring). Instead, RMGPynta shells out to
    # each interpreter as a *subprocess*, running a small driver script
    # that lives in rmgpynta/drivers/ and does the actual `import pynta`
    # or `import rmgpy` -- inside that environment, not this one.
    pynta_env_python=HOME / "miniconda3" / "envs" / "pynta_env" / "bin" / "python",
    rmg_env_python=HOME / "miniconda3" / "envs" / "rmg_env" / "bin" / "python",

    # points at a real ZONE 1 file -- read as text/DSL by RMG, as above
    rmg_input_template=HOME / "rmg-production" / "rmg_ch4_pt111" / "input.py",
    rmg_database_path=HOME / "RMG-Py" / "RMG-database",

    # ZONE 2 files get generated automatically per round -- nothing to
    # point at here; RMGPynta creates reaction.yaml itself, every round
    pynta_kwargs=dict(metal="Pt", surface_type="fcc111", repeats=(3, 3, 4)),

    max_rounds=5,
)

# So: three files, three formats, one direction of truth about imports --
# ZONE 1 (RMG DSL) imports nothing, ZONE 2 (YAML data) imports nothing,
# ZONE 3 (this file) is where `from rmgpynta import RMGPynta` lives, and
# `import pynta` / `import rmgpy` happen one level further down still, in
# rmgpynta/drivers/*.py, each inside its own subprocess.
if __name__ == "__main__":
    orchestrator.run_loop()
