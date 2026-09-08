"""rmgpynta: orchestrator for a closed feedback loop between RMG and Pynta.

RMG (Reaction Mechanism Generator) and Pynta (automated DFT transition-state
search) already talk in one direction on this machine: every completed Pynta
job's postprocessing step writes RMG-format thermo/kinetics libraries
(``thermo_library.py`` / ``reaction_library/``) into its own run directory
(see ``pynta.tasks.PostprocessingTask`` -> ``pynta.postprocessing.write_rmg_libraries``).
Nothing currently reads those libraries back into an RMG ``input.py``.

This package closes that loop: RMG proposes a core (and, each round, an
edge), a filtered subset of the edge is translated into a Pynta
``reaction.yaml``, Pynta computes it, the resulting libraries are registered
back with RMG, and RMG reruns with better-informed libraries.

Pynta and RMG-Py live in separate conda environments (``pynta_env``,
``rmg_env``) on this machine and do not import each other -- Pynta's own
kinetics/thermo classes come from the standalone ``molecule``/``pysidt``
packages, not ``rmgpy``. So this package never imports ``pynta`` or
``rmgpy`` directly; it drives each one as a subprocess via the small,
self-contained scripts in :mod:`rmgpynta.drivers`, which run *inside* the
relevant environment and are the only code in this repo that needs
``pynta`` or ``rmgpy`` importable.
"""

from rmgpynta.orchestrator import RMGPynta, RoundResult

__version__ = "0.1.0"
__all__ = ["RMGPynta", "RoundResult"]
