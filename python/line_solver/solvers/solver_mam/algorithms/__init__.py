"""
MAM Algorithm implementations.

This module contains the core solving algorithms for SolverMAM:
- dec_source: Decomposition with MMAP arrivals (default)
- dec_mmap: Service-scaled departure processes
- dec_poisson: Poisson approximation
- mna_open: Matrix-normalizing approximation (open networks)
- mna_closed: Matrix-normalizing approximation (closed networks)
- ldqbd_solver: Level-dependent QBD (single-class closed)
- rcat_*: RCAT and AutoCAT methods
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional, Dict, Tuple
import numpy as np


@dataclass
class MAMResult:
    """Result structure for MAM algorithms.

    Attributes:
        QN: (M, K) Queue lengths
        UN: (M, K) Server utilizations
        RN: (M, K) Response times
        TN: (1, K) Throughputs
        CN: (1, K) Cycle times (closed networks)
        XN: (1, K) System throughputs
        totiter: Total iterations
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


class MAMAlgorithm(ABC):
    """Abstract base class for MAM algorithms."""

    @abstractmethod
    def solve(self, sn, options) -> MAMResult:
        """Solve the network.

        Args:
            sn: NetworkStruct
            options: Solver options

        Returns:
            MAMResult with QN, UN, RN, TN metrics
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


from .dec_source import DecSourceAlgorithm, DecSourceOptions
from .dec_mmap import DecMMAPAlgorithm, DecMMAPOptions
from .dec_poisson import DecPoissonAlgorithm
from .mna_open import MNAOpenAlgorithm
from .mna_closed import MNAClosedAlgorithm
from .ag_builder import RCATModel, RCATSolver, build_rcat_model
from .ag_inap import INAPAlgorithm, INAPPlusAlgorithm

__all__ = [
    'MAMResult',
    'MAMAlgorithm',
    'DecSourceAlgorithm',
    'DecSourceOptions',
    'DecMMAPAlgorithm',
    'DecMMAPOptions',
    'DecPoissonAlgorithm',
    'MNAOpenAlgorithm',
    'MNAClosedAlgorithm',
    'RCATModel',
    'RCATSolver',
    'build_rcat_model',
    'INAPAlgorithm',
    'INAPPlusAlgorithm',
]
