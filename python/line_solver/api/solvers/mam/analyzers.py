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
from .mmap_fj import solver_mam_basic_mmap


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
        - 'default': Automatic method selection (ldqbd on a single-class closed
          Delay+Queue, else dec.source)
        - 'dec.source': Source decomposition
        - 'dec.poisson': Poisson approximation
        - 'dec.mmap': MMAP decomposition
        - 'dec.source.mmap': MMAP fork-join decomposition with mmap_max
          synchronization at the joins
        - 'mna': Matrix-normalizing approximation (auto-selects mna_open or
          mna_closed; rejects mixed models)
        - 'ldqbd': Level-dependent QBD (single-class Delay/Queue networks)
        - 'inap': RCAT iterative numerical approximation
        - 'inapplus': RCAT INAP+ (weighted rates, no normalization)
        - 'inapinf': RCAT INAP with matrix-geometric isolated components
        - 'exact': Exact analysis via RCAT

    Args:
        sn: Network structure
        options: Solver options (method, tolerance, verbosity)

    Returns:
        MAMResult with all performance metrics

    Raises:
        RuntimeError: For unsupported configurations (e.g., mna method with mixed models)
    """
    start_time = time.time()

    if options is None:
        options = SolverMAMOptions()

    # Set ETAQA truncation default if not specified
    if not hasattr(options, 'etaqa_trunc') or options.etaqa_trunc is None:
        options.etaqa_trunc = 8

    # Discrete-time (slotted) models are recognized from the distributions and
    # routed to the Q-MAM discrete-time algorithms. The test runs before any
    # phase-type conversion, which would fit a continuous surrogate to a
    # Geometric and erase the lattice; see _kb/06-solver-catalog.md for the
    # LAS-DA convention.
    from ...sn.predicates import sn_is_discrete_time
    is_dt, slot_length, _dt_info = sn_is_discrete_time(sn, options)
    if is_dt:
        from .dt import solver_mam_dt
        dt_ret = solver_mam_dt(sn, options, slot_length)
        dt_result = MAMResult()
        dt_result.QN = dt_ret.QN
        dt_result.UN = dt_ret.UN
        dt_result.RN = dt_ret.RN
        dt_result.TN = dt_ret.TN
        dt_result.CN = dt_ret.CN
        dt_result.XN = dt_ret.XN
        dt_result.AN = dt_ret.TN.copy()
        dt_result.WN = dt_ret.RN.copy()
        dt_result.iter = dt_ret.totiter
        dt_result.method = dt_ret.method
        dt_result.runtime = time.time() - start_time
        return dt_result

    method = options.method.lower()

    # Check model type
    is_open = sn_is_open_model(sn)
    is_closed = sn_is_closed_model(sn)

    result = MAMResult()
    ret = None

    # Single-class closed Delay+Queue: the LD-QBD is exact (the level-dependent
    # arrival rate (N-n)*lambda captures the population constraint that
    # dec.source only approximates). Mirrors the MATLAB/JAR default dispatch.
    if method == 'default':
        from ....solvers.solver_mam.algorithms import ldqbd_is_closed_delay_queue
        if ldqbd_is_closed_delay_queue(sn):
            method = 'ldqbd'

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

    elif method == 'dec.source.mmap':
        if options.verbose:
            print("Using dec.source.mmap method, calling solver_mam_basic_mmap")
        ret = solver_mam_basic_mmap(sn, options)
        result.method = 'dec.source.mmap'

    elif method in ['mna', 'mna_open', 'mna_closed', 'inap', 'inapplus',
                    'inapinf', 'exact', 'ldqbd']:
        # MNA, the RCAT family (inap/inapplus/inapinf/exact) and the LD-QBD all
        # have dedicated algorithms; mirrors the MATLAB solver_mam_analyzer
        # dispatch to solver_mna_*/solver_mam_ag/solver_mam_ldqbd. Imported
        # lazily: the algorithm package imports back into
        # api.solvers.mam.handler, so a module-level import would cycle.
        if method == 'mna':
            if is_open and is_closed:
                raise RuntimeError(
                    "The mna method in SolverMAM does not support mixed models.")
            method_impl = 'mna_closed' if is_closed else 'mna_open'
        else:
            method_impl = method
        if options.verbose:
            print(f"Using {method_impl} method")
        from ....solvers.solver_mam.algorithms import (
            INAPAlgorithm, INAPPlusAlgorithm, INAPInfAlgorithm, LDQBDAlgorithm,
            MNAOpenAlgorithm, MNAClosedAlgorithm,
        )
        algo_class = {
            'mna_open': MNAOpenAlgorithm,
            'mna_closed': MNAClosedAlgorithm,
            'inap': INAPAlgorithm,
            'exact': INAPAlgorithm,
            'inapplus': INAPPlusAlgorithm,
            'inapinf': INAPInfAlgorithm,
            'ldqbd': LDQBDAlgorithm,
        }[method_impl]
        algo_result = algo_class().solve(sn, options)
        result.QN = algo_result.QN
        result.UN = algo_result.UN
        result.RN = algo_result.RN
        result.TN = algo_result.TN
        result.CN = algo_result.CN
        result.XN = algo_result.XN
        result.AN = algo_result.TN.copy() if algo_result.TN is not None else None
        result.WN = algo_result.RN.copy() if algo_result.RN is not None else None
        result.iter = algo_result.totiter
        result.method = method_impl

    else:
        raise RuntimeError(f"Unknown method: {method}")

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
