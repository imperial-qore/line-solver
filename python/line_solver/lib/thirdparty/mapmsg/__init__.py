"""
MAPMsG - MAP/M/s+G Call Center Model Solver.

Python port of the MAPMsG MATLAB library by O. Gursoy, K. A. Mehr, N. Akar.

This package provides steady-state and first-passage-time analysis for
MAP/M/s+G queues (call center models with MAP arrivals, exponential service,
s servers, and generally distributed patience times).

Reference:
    O. Gursoy, K. A. Mehr, N. Akar, "The MAP/M/s + G Call Center Model with
    Generally Distributed Patience Times: Steady-state Solution and First
    Passage Time Distribution."
"""

from .cme_parameter_calculator import cme_parameter_calculator, MESystem
from .additive_decomposition import additive_decomposition
from .mrmfq_solver import mrmfq_solver
from .mapmsg_compiler import (
    MAPMsGResult,
    solve_steady_state,
    solve_first_passage_virtual,
    solve_first_passage_actual,
    STEADY_STATE,
    FIRST_PASSAGE_VIRTUAL,
    FIRST_PASSAGE_ACTUAL,
)

__all__ = [
    "cme_parameter_calculator",
    "MESystem",
    "additive_decomposition",
    "mrmfq_solver",
    "MAPMsGResult",
    "solve_steady_state",
    "solve_first_passage_virtual",
    "solve_first_passage_actual",
    "STEADY_STATE",
    "FIRST_PASSAGE_VIRTUAL",
    "FIRST_PASSAGE_ACTUAL",
]
