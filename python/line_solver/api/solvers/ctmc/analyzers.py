"""
CTMC Solver analyzers.

Native Python implementation of CTMC solver analyzers that orchestrate
method selection and provide the main entry point for CTMC analysis.

Port from:

"""

import math
import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
import time

from .memory_guard import ctmc_memory_gate, DEFAULT_SAFETY_FRACTION
from ...sn import NetworkStruct, SchedStrategy, NodeType
from .handler import (
    solver_ctmc,
    solver_ctmc_basic,
    SolverCTMCOptions,
    SolverCTMCReturn,
)


@dataclass
class CTMCResult:
    """
    Result of CTMC solver analysis.

    Attributes:
        QN: Mean queue lengths (M x K)
        UN: Utilizations (M x K)
        RN: Response times (M x K)
        TN: Throughputs (M x K)
        CN: Cycle times (1 x K)
        XN: System throughputs (1 x K)
        pi: Steady-state distribution
        infgen: Infinitesimal generator matrix
        space: State space matrix
        runtime: Runtime in seconds
        method: Method used
    """
    QN: Optional[np.ndarray] = None
    UN: Optional[np.ndarray] = None
    RN: Optional[np.ndarray] = None
    TN: Optional[np.ndarray] = None
    CN: Optional[np.ndarray] = None
    XN: Optional[np.ndarray] = None
    pi: Optional[np.ndarray] = None
    infgen: Optional[np.ndarray] = None
    space: Optional[np.ndarray] = None
    runtime: float = 0.0
    method: str = ""


def solver_ctmc_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverCTMCOptions] = None
) -> CTMCResult:
    """
    CTMC Analyzer - main entry point for CTMC analysis.

    Analyzes queueing networks by constructing and solving the
    underlying continuous-time Markov chain.

    Supported methods:
        - 'default': State-space enumeration
        - 'basic': Explicit state enumeration

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        CTMCResult with all performance metrics

    Raises:
        ValueError: For unsupported configurations
    """
    start_time = time.time()

    if options is None:
        options = SolverCTMCOptions()

    method = options.method.lower()
    result = CTMCResult()

    # see _kb/06-solver-catalog.md (CTMC main section) for the memory-gate rationale
    M = sn.nstations
    K = sn.nclasses

    if sn.njobs is not None:
        N = sn.njobs.flatten()
        # worst-case log state-space size (stars-and-bars per class); open
        # classes are truncated at the (already-resolved) cutoff.
        log_nstates = 0.0
        for k in range(K):
            nk = (options.cutoff if np.isinf(N[k]) else int(N[k]))
            # log C(nk + M - 1, nk)
            log_nstates += (math.lgamma(1 + nk + M - 1)
                            - math.lgamma(1 + M - 1)
                            - math.lgamma(1 + nk))

        force = bool(getattr(options, 'force', False))
        safety = float(getattr(options, 'memory_safety_fraction',
                               DEFAULT_SAFETY_FRACTION))
        ok, msg = ctmc_memory_gate(log_nstates, force=force,
                                   verbose=bool(options.verbose),
                                   safety_fraction=safety)
        if not ok:
            raise MemoryError(msg)

    # Select and execute method
    if method in ['default', 'basic']:
        ret = solver_ctmc_basic(sn, options)
        result.method = 'basic'
    else:
        # Unknown method - use basic
        if options.verbose:
            print(f"Warning: Unknown CTMC method '{method}'. Using basic.")
        ret = solver_ctmc_basic(sn, options)
        result.method = 'basic'

    # Copy results
    if ret is not None:
        result.QN = ret.Q
        result.UN = ret.U
        result.RN = ret.R
        result.TN = ret.T
        result.CN = ret.C
        result.XN = ret.X
        result.pi = ret.pi
        result.infgen = ret.infgen
        result.space = ret.space

    # Clean up NaN values
    if result.QN is not None:
        result.QN = np.nan_to_num(result.QN, nan=0.0)
    if result.UN is not None:
        result.UN = np.nan_to_num(result.UN, nan=0.0)
    if result.RN is not None:
        result.RN = np.nan_to_num(result.RN, nan=0.0)
    if result.TN is not None:
        result.TN = np.nan_to_num(result.TN, nan=0.0)
    if result.CN is not None:
        result.CN = np.nan_to_num(result.CN, nan=0.0)
    if result.XN is not None:
        result.XN = np.nan_to_num(result.XN, nan=0.0)

    result.runtime = time.time() - start_time

    return result


__all__ = [
    'CTMCResult',
    'solver_ctmc_analyzer',
]
