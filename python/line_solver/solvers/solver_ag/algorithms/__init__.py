"""Algorithms of the agent-based (RCAT) solver.

Every algorithm here decomposes the network into agents -- one per
(station, class) pair -- and solves the fixed point over the reversed rates of
the synchronizing actions. They differ only in how the reversed rate of an
action is read off the active agent's stationary vector, and in whether an open
agent is truncated or carried to infinity through its matrix-geometric tail.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


@dataclass
class AGResult:
    """Result structure for AG algorithms.

    Attributes:
        QN: (M, K) Queue lengths
        UN: (M, K) Server utilizations
        RN: (M, K) Response times
        TN: (1, K) Throughputs
        CN: (1, K) Cycle times (closed networks)
        XN: (1, K) System throughputs
        totiter: Sweeps taken by the reversed-rate fixed point
        method: Algorithm name
        runtime: Computation time (seconds)
    """
    QN: np.ndarray
    UN: np.ndarray
    RN: np.ndarray
    TN: np.ndarray
    CN: Optional[np.ndarray] = None
    XN: Optional[np.ndarray] = None
    totiter: int = 0
    method: str = ""
    runtime: float = 0.0


class AGAlgorithm(ABC):
    """Abstract base class for AG algorithms."""

    @abstractmethod
    def solve(self, sn, options) -> AGResult:
        """Solve the network.

        Args:
            sn: NetworkStruct
            options: Solver options

        Returns:
            AGResult with QN, UN, RN, TN metrics
        """
        pass

    @staticmethod
    @abstractmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """Check if algorithm can solve this network.

        Args:
            sn: NetworkStruct

        Returns:
            (can_solve, reason_if_not)
        """
        pass


from .ag_builder import RCATModel, RCATSolver, build_rcat_model  # noqa: E402
from .ag_inap import INAPAlgorithm, INAPPlusAlgorithm, INAPInfAlgorithm  # noqa: E402

__all__ = [
    'AGResult',
    'AGAlgorithm',
    'RCATModel',
    'RCATSolver',
    'build_rcat_model',
    'INAPAlgorithm',
    'INAPPlusAlgorithm',
    'INAPInfAlgorithm',
]
