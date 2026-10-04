"""
LDES solver for LINE — subprocess/JSON backend.

Runs simulation via `java -jar ldes.jar solve model.json` and parses
the JSON result, using the same simulation engine as MATLAB and
python-wrapper.

Example usage:
    from line_solver.solvers.solver_ldes import SolverLDES, LDESOptions

    options = LDESOptions(seed=23000, samples=200000)
    solver = SolverLDES(model, options)
    result = solver.runAnalyzer()
    table = solver.getAvgTable()
"""

from .ldes_options import LDESOptions, LDESResult
from .solver_ldes import SolverLDES

__all__ = [
    'LDESOptions',
    'LDESResult',
    'SolverLDES',
]
