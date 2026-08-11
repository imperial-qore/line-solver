"""
Third-party solver wrappers (client side).

These solvers do not implement analysis in native Python; they serialize the
model to an external artifact and shell out to a separate tool, exchanging data
across a serialized (JSON/XML) boundary rather than an in-process object:

- ``solver_jmt``  -> Java Modelling Tools (JMT.jar) via .jsimg
- ``solver_lqns`` -> the ``lqns``/``lqsim`` binaries
- ``solver_qns``  -> the ``qnsolver`` binary (itself calls LQNS)
- ``solver_ldes`` -> the LINE Discrete Event Simulator (ldes.jar / native
  binary) via ``model.json`` in and ``result.json`` out (fully JSON-mediated)

Kept here to separate results that depend on an external installation from the
native analytical/simulation engines. Public import paths are unchanged: the
``Solver*`` classes are re-exported from ``line_solver.solvers`` and
``line_solver``.
"""

from .solver_jmt import SolverJMT, SolverJMTOptions
from .solver_lqns import SolverLQNS, LQNSOptions, LQNSResult
from .solver_qns import SolverQNS, QNSOptions, QNSResult
from .solver_ldes.ldes_options import LDESOptions, LDESResult

try:
    from .solver_ldes.solver_ldes import SolverLDES
except ImportError:
    SolverLDES = None

__all__ = [
    'SolverJMT', 'SolverJMTOptions',
    'SolverLQNS', 'LQNSOptions', 'LQNSResult',
    'SolverQNS', 'QNSOptions', 'QNSResult',
    'SolverLDES', 'LDESOptions', 'LDESResult',
]
