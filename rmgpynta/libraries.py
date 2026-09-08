"""Steps 4 & 5: collect Pynta's RMG-format libraries and register them with RMG.

Step 4 (Pynta writes ``thermo_library.py`` / ``reaction_library/``) already
happens automatically today via ``pynta.tasks.PostprocessingTask`` ->
``pynta.postprocessing.write_rmg_libraries`` -- see
``rmgpynta.drivers.postprocess_fix`` for the corrected reimplementation this
project uses instead of trusting that function's output directly (it has a
confirmed bug that drops all but one reaction from ``reactions.py``).

Step 5 (registering those libraries with RMG) is the actual gap this
project closes: nothing on this machine currently copies a Pynta-written
library into ``RMG-database``'s search path or lists it in an RMG
``input.py``. This module does both halves of that.
"""
from __future__ import annotations

import re
import shutil
from pathlib import Path

from rmgpynta.paths import RoundPaths, libraries_root


class LibraryNotFound(FileNotFoundError):
    pass


def collect_pynta_libraries(rp: RoundPaths) -> tuple:
    """Locate this round's Pynta-written libraries; raise if incomplete.

    Expects ``rmgpynta.drivers.postprocess_fix.write_rmg_libraries_fixed``
    to have already run (it's called automatically at the end of
    ``drivers/run_pynta_job.py``) -- this function only validates and
    returns paths, it doesn't regenerate anything itself.
    """
    thermo = rp.pynta_thermo_library
    kinetics_dir = rp.pynta_reaction_library_dir
    if not thermo.exists():
        raise LibraryNotFound(f"no thermo_library.py under {rp.pynta_dir}")
    if not (kinetics_dir / "reactions.py").exists() or not (kinetics_dir / "dictionary.txt").exists():
        raise LibraryNotFound(f"no complete reaction_library/ under {rp.pynta_dir}")
    return thermo, kinetics_dir


def register_libraries(
    rmg_database_path: Path, thermo_library_path: Path, kinetics_library_dir: Path, library_name: str
) -> str:
    """Copy this round's libraries into RMG-database's own library search path.

    RMG resolves a thermo library named ``foo`` at
    ``<rmg_database_path>/input/thermo/libraries/foo.py`` and a kinetics
    library named ``foo`` at
    ``<rmg_database_path>/input/kinetics/libraries/foo/{reactions.py,dictionary.txt}``
    -- exactly the shape Pynta's own ``write_rmg_libraries`` already
    produces (confirmed against ``~/rmg-production/rmg_ch4_pt111/input.py``'s
    ``thermoLibraries=['surfaceThermoPt111', 'primaryThermoLibrary', ...]``
    and the standard RMG-database layout), which is why this is a copy, not
    a format conversion.
    """
    rmg_database_path = Path(rmg_database_path)
    thermo_dest = rmg_database_path / "input" / "thermo" / "libraries" / f"{library_name}.py"
    kinetics_dest = rmg_database_path / "input" / "kinetics" / "libraries" / library_name

    thermo_dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(thermo_library_path, thermo_dest)

    if kinetics_dest.exists():
        shutil.rmtree(kinetics_dest)
    shutil.copytree(kinetics_library_dir, kinetics_dest)

    # Also keep a round-tagged copy under work_dir/libraries/ purely for
    # human inspection / reproducibility (paths.libraries_root) -- this is
    # separate from, and not read by, RMG itself.
    return library_name


def archive_round_libraries(work_dir: Path, rp: RoundPaths, library_name: str) -> None:
    root = libraries_root(work_dir)
    shutil.copy2(rp.pynta_thermo_library, root / "thermo" / f"{library_name}.py")
    dest = root / "kinetics" / library_name
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(rp.pynta_reaction_library_dir, dest)


_LIST_LITERAL_RE_TEMPLATE = r"({name}\s*=\s*\[)(.*?)(\])"


def _insert_into_list_literal(text: str, list_name: str, new_items: list) -> str:
    pattern = re.compile(_LIST_LITERAL_RE_TEMPLATE.format(name=re.escape(list_name)), re.DOTALL)
    m = pattern.search(text)
    if not m:
        raise ValueError(
            f"could not find `{list_name}=[...]` in the input.py template -- "
            f"is this a real RMG input.py with a database(...) block?"
        )
    existing_body = m.group(2)
    to_add = [item for item in new_items if f"'{item}'" not in existing_body and f'"{item}"' not in existing_body]
    insertion = "".join(f"\n        '{item}'," for item in to_add)
    replacement = m.group(1) + insertion + existing_body + m.group(3)
    return text[: m.start()] + replacement + text[m.end() :]


def update_input_py_libraries(
    template_path: Path,
    out_path: Path,
    thermo_library_names: list,
    kinetics_library_names: list,
) -> Path:
    """Splice new library names into a copy of an RMG ``input.py`` template.

    New names are inserted at the front of each list literal. NOTE: RMG's
    actual priority order for ``thermoLibraries``/``reactionLibraries`` (does
    an earlier-listed library win, or a later one?) was not independently
    confirmed while writing this module -- if round-to-round library
    priority matters for your system, verify this against RMG-Py's database
    loading order before relying on "newest round wins" being true.
    """
    text = Path(template_path).read_text()
    text = _insert_into_list_literal(text, "thermoLibraries", thermo_library_names)
    text = _insert_into_list_literal(text, "reactionLibraries", kinetics_library_names)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text)
    return out_path
