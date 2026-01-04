"""
NC Solver analyzers.

Native Python implementation of NC solver analyzers that orchestrate the
core handlers and provide interpolation for non-integer populations.

Port from:


"""

import numpy as np
from math import floor, ceil
from typing import Optional
from dataclasses import dataclass, field
import time

from ...sn import NetworkStruct
from .handler import solver_nc, solver_ncld, SolverNCReturn, SolverNCLDReturn, SolverOptions


# Fine tolerance for numerical comparisons
FINE_TOL = 1e-12


@dataclass
class NCResultProb:
    """Probability results from NC solver."""
    logNormConstAggr: Optional[float] = None
    marginal: Optional[np.ndarray] = None
    joint: Optional[float] = None
    itemProb: Optional[np.ndarray] = None


@dataclass
class NCResult:
    """
    Result of NC solver analysis.

    Attributes:
        QN: Mean queue lengths (M x K)
        UN: Utilizations (M x K)
        RN: Response times (M x K)
        TN: Throughputs (M x K)
        CN: Cycle times (1 x K)
        XN: System throughputs (1 x K)
        lG: Log normalizing constant
        STeff: Effective service times
        it: Number of iterations
        runtime: Runtime in seconds
        method: Method used
        solver: Solver name
        prob: Probability results
    """
    QN: Optional[np.ndarray] = None
    UN: Optional[np.ndarray] = None
    RN: Optional[np.ndarray] = None
    TN: Optional[np.ndarray] = None
    CN: Optional[np.ndarray] = None
    XN: Optional[np.ndarray] = None
    lG: float = 0.0
    STeff: Optional[np.ndarray] = None
    it: int = 0
    runtime: float = 0.0
    method: str = ""
    solver: str = "NC"
    prob: NCResultProb = field(default_factory=NCResultProb)


def solver_nc_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverOptions] = None
) -> NCResult:
    """
    Main NC solver analyzer.

    Performs NC analysis with interpolation for non-integer populations.
    If population values are not integers, interpolates between floor and
    ceiling values.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        NCResult with all performance metrics

    Raises:
        RuntimeError: For unsupported configurations
    """
    start_time = time.time()

    if options is None:
        options = SolverOptions()

    # Check for multiserver with exact method
    nservers = sn.nservers
    if nservers is not None:
        nservers_finite = nservers.copy()
        nservers_finite = nservers_finite[np.isfinite(nservers_finite)]
        if len(nservers_finite) > 0 and np.max(nservers_finite) > 1 and options.method == 'exact':
            raise RuntimeError(
                "NC solver cannot provide exact solutions for open or mixed queueing networks. "
                "Remove the 'exact' option."
            )

    # Create floor/ceiling copies for non-integer populations
    njobs = sn.njobs
    if njobs is None:
        njobs = np.zeros((1, sn.nclasses))

    njobs_floor = np.floor(njobs)
    njobs_ceil = np.ceil(njobs)
    eta = np.abs(njobs - njobs_floor)

    # Check for non-integer populations
    non_integer_job = np.any(eta > FINE_TOL)

    result = NCResult()

    if non_integer_job:
        # Interpolate between floor and ceiling
        sn_floor = sn.copy()
        sn_ceil = sn.copy()
        sn_floor.njobs = njobs_floor
        sn_ceil.njobs = njobs_ceil

        ret_floor = solver_nc(sn_floor, options)
        ret_ceil = solver_nc(sn_ceil, options)

        result.runtime = ret_floor.runtime + ret_ceil.runtime

        # Interpolate results
        if ret_floor.Q is not None and ret_ceil.Q is not None:
            result.QN = ret_floor.Q + eta * (ret_ceil.Q - ret_floor.Q)
        if ret_floor.U is not None and ret_ceil.U is not None:
            result.UN = ret_floor.U + eta * (ret_ceil.U - ret_floor.U)
        if ret_floor.R is not None and ret_ceil.R is not None:
            result.RN = ret_floor.R + eta * (ret_ceil.R - ret_floor.R)
        if ret_floor.T is not None and ret_ceil.T is not None:
            result.TN = ret_floor.T + eta * (ret_ceil.T - ret_floor.T)
        if ret_floor.X is not None and ret_ceil.X is not None:
            result.XN = ret_floor.X + eta * (ret_ceil.X - ret_floor.X)

        # Interpolate lG
        result.lG = ret_floor.lG + np.sum(eta) * (ret_floor.lG - ret_ceil.lG)

        result.it = ret_floor.it + ret_ceil.it
        result.method = ret_ceil.method
    else:
        # Integer populations - direct computation
        ret = solver_nc(sn, options)
        result.QN = ret.Q.copy() if ret.Q is not None else None
        result.UN = ret.U.copy() if ret.U is not None else None
        result.RN = ret.R.copy() if ret.R is not None else None
        result.TN = ret.T.copy() if ret.T is not None else None
        result.XN = ret.X.copy() if ret.X is not None else None
        result.lG = ret.lG
        result.it = ret.it
        result.method = ret.method

    # Calculate cycle times using Little's Law: C(k) = N(k) / X(k)
    if result.XN is not None and njobs is not None:
        K = sn.nclasses
        result.CN = np.zeros((1, K))
        for k in range(K):
            njobs_flat = njobs.flatten()
            xn_flat = result.XN.flatten()
            if k < len(njobs_flat) and k < len(xn_flat) and xn_flat[k] > 0:
                result.CN[0, k] = njobs_flat[k] / xn_flat[k]

    result.runtime = time.time() - start_time
    result.solver = "NC"

    return result


def solver_ncld_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverOptions] = None
) -> NCResult:
    """
    Load-dependent NC solver analyzer.

    Performs load-dependent NC analysis with interpolation for non-integer
    populations.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        NCResult with all performance metrics

    Raises:
        RuntimeError: For unsupported configurations
    """
    start_time = time.time()

    if options is None:
        options = SolverOptions()

    # Check for multiserver with exact method and open classes
    nservers = sn.nservers
    if nservers is not None:
        nservers_finite = nservers.copy()
        nservers_finite = nservers_finite[np.isfinite(nservers_finite)]
        has_infinite_jobs = sn.njobs is not None and np.any(np.isinf(sn.njobs))
        if (len(nservers_finite) > 0 and np.max(nservers_finite) > 1 and
                has_infinite_jobs and options.method == 'exact'):
            raise RuntimeError(
                "NC solver cannot provide exact solutions for open or mixed queueing networks. "
                "Remove the 'exact' option."
            )

    # Create floor/ceiling copies for non-integer populations
    njobs = sn.njobs
    if njobs is None:
        njobs = np.zeros((1, sn.nclasses))

    njobs_floor = np.floor(njobs)
    njobs_ceil = np.ceil(njobs)
    eta = np.abs(njobs - njobs_floor)

    # Check for non-integer populations
    non_integer_job = np.any(eta > FINE_TOL)

    result = NCResult()

    if non_integer_job:
        if options.method == 'exact':
            raise RuntimeError(
                "NC load-dependent solver cannot provide exact solutions for fractional populations."
            )

        # Interpolate between floor and ceiling
        sn_floor = sn.copy()
        sn_ceil = sn.copy()
        sn_floor.njobs = njobs_floor
        sn_ceil.njobs = njobs_ceil

        ret_floor = solver_ncld(sn_floor, options)
        ret_ceil = solver_ncld(sn_ceil, options)

        # Interpolate results
        if ret_floor.Q is not None and ret_ceil.Q is not None:
            result.QN = ret_floor.Q + eta * (ret_ceil.Q - ret_floor.Q)
        if ret_floor.U is not None and ret_ceil.U is not None:
            result.UN = ret_floor.U + eta * (ret_ceil.U - ret_floor.U)
        if ret_floor.R is not None and ret_ceil.R is not None:
            result.RN = ret_floor.R + eta * (ret_ceil.R - ret_floor.R)
        if ret_floor.T is not None and ret_ceil.T is not None:
            result.TN = ret_floor.T + eta * (ret_ceil.T - ret_floor.T)
        if ret_floor.X is not None and ret_ceil.X is not None:
            result.XN = ret_floor.X + eta * (ret_ceil.X - ret_floor.X)

        # Interpolate lG
        result.lG = ret_floor.lG + np.sum(eta) * (ret_floor.lG - ret_ceil.lG)

        result.it = ret_floor.it + ret_ceil.it
        result.method = ret_ceil.method
    else:
        # Integer populations - direct computation
        ret = solver_ncld(sn, options)
        result.QN = ret.Q
        result.UN = ret.U
        result.RN = ret.R
        result.TN = ret.T
        result.XN = ret.X
        result.lG = ret.lG
        result.it = ret.it
        result.method = ret.method

    # Calculate cycle times using Little's Law: C(k) = N(k) / X(k)
    if result.XN is not None and njobs is not None:
        K = sn.nclasses
        result.CN = np.zeros((1, K))
        for k in range(K):
            njobs_flat = njobs.flatten()
            xn_flat = result.XN.flatten()
            if k < len(njobs_flat) and k < len(xn_flat) and xn_flat[k] > 0:
                result.CN[0, k] = njobs_flat[k] / xn_flat[k]

    result.runtime = time.time() - start_time
    result.solver = "NCLD"

    return result


__all__ = [
    'NCResult',
    'NCResultProb',
    'solver_nc_analyzer',
    'solver_ncld_analyzer',
]
