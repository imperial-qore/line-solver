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
        - 'parallel': Worker-count-invariant replica averaging (R=8 seeded replicas)
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

    # see _kb/06-solver-catalog.md (SSA main section) -- SPN routes to NRM directly
    from ...sn import NodeType as _NodeType
    _is_spn = any(sn.nodetype[i] == _NodeType.TRANSITION for i in range(sn.nnodes))

    # Select and execute method
    if _is_spn:
        from .nrm import solver_ssa_nrm
        ret = solver_ssa_nrm(sn, options)
        result.method = 'nrm'
    elif method == 'nrm':
        # see _kb/06-solver-catalog.md (SSA: "NRM engine internals")
        from .nrm import (solver_ssa_nrm, _ALLOWED_SCHED, _fcr_nrm_ok,
                          _impatience_nrm_ok, _phase_nrm_ok,
                          _cache_nrm_ok)
        supported = all(sn.sched[ist] in _ALLOWED_SCHED for ist in range(sn.nstations))
        # see _kb/06-solver-catalog.md (SSA: "NRM engine now supports FCR directly")
        if supported and not _fcr_nrm_ok(sn):
            supported = False
        # QUEUE_LENGTH balking and memoryless reneging only
        if supported and not _impatience_nrm_ok(sn):
            supported = False
        # see _kb/06-solver-catalog.md (SSA main section) -- phase-type
        # expansion at INF/PS and non-preemptive buffered families only
        if supported and not _phase_nrm_ok(sn):
            supported = False
        # Cache nodes are simulated natively; only the retrieval (delayed-hit)
        # system still needs the serial engine.
        if supported and not _cache_nrm_ok(sn):
            supported = False
        if supported:
            ret = solver_ssa_nrm(sn, options)
            result.method = 'nrm'
        else:
            ret = solver_ssa_basic(sn, options)
            result.method = 'serial'
    elif method in ['default', 'serial']:
        ret = solver_ssa_basic(sn, options)
        result.method = 'serial'
    elif method in ['parallel', 'para', 'ssa.parallel']:
        # Worker-count-invariant replica averaging (mirrors MATLAB
        # solver_ssa_analyzer_parallel); see handler.solver_ssa_parallel.
        from .handler import solver_ssa_parallel
        ret = solver_ssa_parallel(sn, options)
        result.method = 'parallel'
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
