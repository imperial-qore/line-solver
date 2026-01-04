"""
MAM Solver analyzers.

Native Python implementation of MAM solver analyzers that orchestrate
method selection and provide the main entry point for matrix-analytic analysis.

Port from:


"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
import time

from ...sn import (
    NetworkStruct,
    SchedStrategy,
    NodeType,
    sn_is_open_model,
    sn_is_closed_model,
)
from .handler import solver_mam, solver_mam_basic, SolverMAMOptions, SolverMAMReturn


@dataclass
class MAMResult:
    """
    Result of MAM solver analysis.

    Attributes:
        QN: Mean queue lengths (M x K)
        UN: Utilizations (M x K)
        RN: Response times (M x K)
        TN: Throughputs (M x K)
        CN: Cycle times (1 x K)
        XN: System throughputs (1 x K)
        AN: Arrival rates (M x K)
        WN: Waiting times (M x K)
        iter: Number of iterations
        runtime: Runtime in seconds
        method: Method used
        lG: Log normalization constant (if applicable)
    """
    QN: Optional[np.ndarray] = None
    UN: Optional[np.ndarray] = None
    RN: Optional[np.ndarray] = None
    TN: Optional[np.ndarray] = None
    CN: Optional[np.ndarray] = None
    XN: Optional[np.ndarray] = None
    AN: Optional[np.ndarray] = None
    WN: Optional[np.ndarray] = None
    iter: int = 0
    runtime: float = 0.0
    method: str = ""
    lG: float = 0.0


def _has_fcfs_scheduling(sn: NetworkStruct) -> bool:
    """
    Check if network has FCFS scheduling.

    Args:
        sn: Network structure

    Returns:
        True if any station has FCFS scheduling
    """
    sched_dict = sn.sched if sn.sched else {}
    for i in range(sn.nstations):
        station_sched = sched_dict.get(i)
        if station_sched == SchedStrategy.FCFS:
            return True
    return False


def _has_ps_scheduling(sn: NetworkStruct) -> bool:
    """
    Check if network has PS (Processor Sharing) scheduling.

    Args:
        sn: Network structure

    Returns:
        True if any station has PS scheduling
    """
    sched_dict = sn.sched if sn.sched else {}
    for i in range(sn.nstations):
        station_sched = sched_dict.get(i)
        if station_sched == SchedStrategy.PS:
            return True
    return False


def _has_hol_scheduling(sn: NetworkStruct) -> bool:
    """
    Check if network has HOL (Head-of-Line) scheduling.

    Args:
        sn: Network structure

    Returns:
        True if any station has HOL scheduling
    """
    sched_dict = sn.sched if sn.sched else {}
    for i in range(sn.nstations):
        station_sched = sched_dict.get(i)
        if station_sched == SchedStrategy.HOL:
            return True
    return False


def solver_mam_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverMAMOptions] = None
) -> MAMResult:
    """
    MAM Analyzer - main entry point for matrix-analytic analysis.

    Selects appropriate MAM method based on network characteristics
    and solver options, then performs the analysis.

    Supported methods:
        - 'default': Automatic method selection (dec.source)
        - 'dec.source': Source decomposition
        - 'dec.poisson': Poisson approximation
        - 'dec.mmap': MMAP decomposition
        - 'mna': Matrix-analytic method
        - 'inap': Iterative numerical approximation
        - 'exact': Exact analysis via RCAT

    Args:
        sn: Network structure
        options: Solver options (method, tolerance, verbosity)

    Returns:
        MAMResult with all performance metrics

    Raises:
        RuntimeError: For unsupported configurations (e.g., mixed models)
    """
    start_time = time.time()

    if options is None:
        options = SolverMAMOptions()

    method = options.method.lower()

    # Check for mixed model (not supported)
    is_open = sn_is_open_model(sn)
    is_closed = sn_is_closed_model(sn)

    if is_open and is_closed:
        raise RuntimeError("SolverMAM does not support mixed models with both open and closed classes.")

    result = MAMResult()
    ret = None

    # Select and execute method
    if method == 'default' or method == 'dec.source':
        if options.verbose:
            print("Using dec.source method, calling solver_mam_basic")
        ret = solver_mam_basic(sn, options)
        result.method = 'dec.source'

    elif method == 'dec.poisson':
        if options.verbose:
            print("Using dec.poisson method with space_max=1")
        options.space_max = 1
        ret = solver_mam_basic(sn, options)
        result.method = 'dec.poisson'

    elif method == 'dec.mmap':
        if options.verbose:
            print("Using dec.mmap method")
        ret = solver_mam(sn, options)
        result.method = 'dec.mmap'

    elif method == 'mna':
        if options.verbose:
            if is_closed:
                print("Using MNA method for closed model")
            else:
                print("Using MNA method for open model")
        ret = solver_mam(sn, options)
        result.method = 'mna'

    elif method in ['inap', 'exact']:
        if options.verbose:
            print(f"Using RCAT method: {method}")
        ret = solver_mam(sn, options)
        result.method = method

    else:
        # Unsupported method - use default
        if options.verbose:
            print(f"Warning: Unsupported MAM method '{method}'. Using dec.source.")
        ret = solver_mam_basic(sn, options)
        result.method = 'dec.source'

    # Copy results from handler return
    if ret is not None:
        result.QN = ret.Q
        result.UN = ret.U
        result.RN = ret.R
        result.TN = ret.T
        result.CN = ret.C
        result.XN = ret.X
        result.AN = ret.A if ret.A is not None else (ret.T.copy() if ret.T is not None else None)
        result.WN = ret.W if ret.W is not None else (ret.R.copy() if ret.R is not None else None)
        result.iter = ret.it

    # Handle external arrivals (source stations)
    if sn.sched is not None and result.TN is not None:
        rates = np.asarray(sn.rates) if hasattr(sn, 'rates') and sn.rates is not None else None
        if rates is not None:
            sched_dict = sn.sched if isinstance(sn.sched, dict) else {}
            for i in range(sn.nstations):
                station_sched = sched_dict.get(i)
                if station_sched == SchedStrategy.EXT:
                    for k in range(sn.nclasses):
                        if i < rates.shape[0] and k < rates.shape[1]:
                            result.TN[i, k] = rates[i, k]

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
    'MAMResult',
    'solver_mam_analyzer',
]
