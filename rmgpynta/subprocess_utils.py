"""Run the standalone scripts in ``rmgpynta.drivers`` inside pynta_env/rmg_env.

Pynta and RMG-Py live in separate conda environments and don't import each
other (see the package docstring in ``rmgpynta/__init__.py``), so every
step that needs either one runs as a subprocess of that environment's own
Python interpreter, executing a self-contained driver script that imports
only stdlib plus ``pynta`` or ``rmgpy`` -- never ``rmgpynta`` itself, since
this package isn't assumed to be installed in either environment.
"""
from __future__ import annotations

import subprocess
from pathlib import Path


class DriverError(RuntimeError):
    def __init__(self, cmd, result):
        self.cmd = cmd
        self.result = result
        super().__init__(
            f"{' '.join(str(c) for c in cmd)} failed (exit {result.returncode})\n"
            f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
        )


def run_in_env(
    python_exe: Path,
    script_path: Path,
    *args,
    cwd: Path = None,
    timeout: float = None,
    check: bool = True,
) -> subprocess.CompletedProcess:
    cmd = [str(python_exe), str(script_path), *[str(a) for a in args]]
    result = subprocess.run(
        cmd, cwd=str(cwd) if cwd else None, timeout=timeout, capture_output=True, text=True
    )
    if check and result.returncode != 0:
        raise DriverError(cmd, result)
    return result
