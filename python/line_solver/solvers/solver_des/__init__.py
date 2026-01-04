"""
SimPy-based Discrete Event Simulation (DES) solver for LINE.

This package provides a native Python DES solver using SimPy, achieving
full parity with the Java/Kotlin DES implementation.

Features:
- 20+ scheduling disciplines (FCFS, LCFS, PS, DPS, GPS, SRPT, etc.)
- Phase-type distributions (PH, APH, MAP, BMAP)
- Fork-join networks with quorum support
- MSER-5 warmup detection
- Overlapping batch means confidence intervals
- Multiclass, multiserver queueing networks

Example usage:
    from line_solver.solvers.solver_des import SolverDES, DESOptions

    options = DESOptions(seed=23000, samples=200000)
    solver = SolverDES(model, options)
    result = solver.runAnalyzer()
    table = solver.getAvgTable()
"""

from .des_options import DESOptions, DESResult

# Lazy imports for simpy-dependent modules to allow importing submodules
# (like scheduling) without requiring simpy to be installed
def __getattr__(name):
    """Lazy import for simpy-dependent components."""
    if name == 'SolverDES':
        from .solver_des import SolverDES
        return SolverDES
    elif name == 'SimPySimulator':
        from .simulator import SimPySimulator
        return SimPySimulator
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    'DESOptions',
    'DESResult',
    'SolverDES',
    'SimPySimulator',
]
