"""
SSA Solver analyzers.

Native Python implementation of SSA solver analyzers that orchestrate
method selection and provide the main entry point for simulation analysis.

Port from:

"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
import time

from ...sn import NetworkStruct, SchedStrategy, NodeType
from .handler import (
    solver_ssa,
    solver_ssa_basic,
    SolverSSAOptions,
    SolverSSAReturn,
)


@dataclass
class SSAResult:
    """
    Result of SSA solver analysis.

    Attributes:
        QN: Mean queue lengths (M x K)
        UN: Utilizations (M x K)
        RN: Response times (M x K)
        TN: Throughputs (M x K)
        CN: Cycle times (1 x K)
        XN: System throughputs (1 x K)
        Q_ci: Queue length confidence intervals
        U_ci: Utilization confidence intervals
        R_ci: Response time confidence intervals
        T_ci: Throughput confidence intervals
        total_time: Total simulated time
        samples: Number of samples collected
        runtime: Runtime in seconds
        method: Method used
    """
    QN: Optional[np.ndarray] = None
    UN: Optional[np.ndarray] = None
    RN: Optional[np.ndarray] = None
    TN: Optional[np.ndarray] = None
    CN: Optional[np.ndarray] = None
    XN: Optional[np.ndarray] = None
    Q_ci: Optional[np.ndarray] = None
    U_ci: Optional[np.ndarray] = None
    R_ci: Optional[np.ndarray] = None
    T_ci: Optional[np.ndarray] = None
    total_time: float = 0.0
    samples: int = 0
    runtime: float = 0.0
    method: str = ""


def solver_ssa_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverSSAOptions] = None
) -> SSAResult:
    """
    SSA Analyzer - main entry point for simulation analysis.

    Analyzes queueing networks using discrete-event simulation
    with the Stochastic Simulation Algorithm (Gillespie method).

    Supported methods:
        - 'default': Serial Gillespie simulation
        - 'serial': Serial simulation (explicit)
        - 'parallel': Parallel simulation (not yet implemented)
        - 'nrm': Next Reaction Method

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SSAResult with all performance metrics

    Raises:
        ValueError: For unsupported configurations
    """
    start_time = time.time()

    if options is None:
        options = SolverSSAOptions()

    method = options.method.lower()
    result = SSAResult()

    # Select and execute method
    if method in ['default', 'serial', 'nrm']:
        ret = solver_ssa_basic(sn, options)
        result.method = 'serial'
    elif method in ['parallel', 'ssa.parallel']:
        # Not yet implemented, fall back to serial
        if options.verbose:
            print("Warning: Parallel SSA not implemented. Using serial.")
        ret = solver_ssa_basic(sn, options)
        result.method = 'serial'
    else:
        # Unknown method - use serial
        if options.verbose:
            print(f"Warning: Unknown SSA method '{method}'. Using serial.")
        ret = solver_ssa_basic(sn, options)
        result.method = 'serial'

    # Copy results
    if ret is not None:
        result.QN = ret.Q
        result.UN = ret.U
        result.RN = ret.R
        result.TN = ret.T
        result.CN = ret.C
        result.XN = ret.X
        result.Q_ci = ret.Q_ci
        result.U_ci = ret.U_ci
        result.R_ci = ret.R_ci
        result.T_ci = ret.T_ci
        result.total_time = ret.total_time
        result.samples = ret.samples

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
    'SSAResult',
    'solver_ssa_analyzer',
]
