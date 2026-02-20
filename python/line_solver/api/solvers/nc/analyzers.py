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


def solver_nc_lossn_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverOptions] = None
) -> NCResult:
    """
    NC solver analyzer for open loss networks with FCR.

    Analyzes open queueing networks with a single multiclass Delay node
    inside a Finite Capacity Region (FCR) with DROP policy using the
    Erlang fixed-point approximation.

    Args:
        sn: Network structure with FCR configuration
        options: Solver options

    Returns:
        NCResult with all performance metrics

    Reference:
        MATLAB: solver_nc_lossn_analyzer.m
    """
    from ...lossn import lossn_erlangfp

    start_time = time.time()

    if options is None:
        options = SolverOptions()

    K = sn.nclasses
    M = sn.nstations

    # 1. Extract arrival rates from Source
    nu = np.zeros(K)
    rates = sn.rates
    if rates is None:
        raise RuntimeError("Network structure has no rates defined")

    for r in range(K):
        refstat = int(sn.refstat[r]) if hasattr(sn, 'refstat') and sn.refstat is not None else 0
        nu[r] = rates[refstat, r]  # arrival rate at source

    # 2. Find delay station in FCR and extract constraints
    region_matrix = sn.region[0]  # First (and only) FCR

    # Find stations with FCR constraints (non-negative values)
    stations_in_fcr = []
    for i in range(M):
        # Check if any per-class constraint or global constraint is non-negative
        has_class_constraint = np.any(region_matrix[i, :K] >= 0)
        has_global_constraint = region_matrix[i, K] >= 0 if region_matrix.shape[1] > K else False
        if has_class_constraint or has_global_constraint:
            stations_in_fcr.append(i)

    if not stations_in_fcr:
        raise RuntimeError("No stations found in FCR")

    delay_idx = stations_in_fcr[0]

    # Extract global and per-class max jobs
    global_max = region_matrix[delay_idx, K] if region_matrix.shape[1] > K else -1
    class_max = region_matrix[delay_idx, :K].copy()

    # Handle unbounded constraints (replace -1 with large value)
    if global_max < 0:
        global_max = 1e6
    class_max[class_max < 0] = 1e6

    # 3. Build A matrix (J x K) where J = K+1 links
    # Link 0: global constraint (all classes contribute)
    # Links 1..K: per-class constraints (only class r contributes to link r+1)
    J = K + 1
    A = np.zeros((J, K))
    A[0, :] = 1.0  # global link: all classes contribute
    for r in range(K):
        A[r + 1, r] = 1.0  # per-class link: only class r contributes

    # 4. Build C vector (J,)
    C_vec = np.zeros(J)
    C_vec[0] = global_max
    C_vec[1:] = class_max

    # 5. Call lossn_erlangfp
    qlen, loss, E, niter = lossn_erlangfp(nu, A, C_vec)

    # 6. Convert to standard outputs
    Q = np.zeros((M, K))
    U = np.zeros((M, K))
    T = np.zeros((M, K))
    R = np.zeros((M, K))

    # At delay node: qlen is effective throughput (after loss)
    for r in range(K):
        T[delay_idx, r] = qlen[r]  # effective throughput = arrival rate * (1 - loss)
        mu_r = rates[delay_idx, r]  # service rate at delay
        if mu_r > 0:
            Q[delay_idx, r] = qlen[r] / mu_r  # Little's law: Q = X * S
            R[delay_idx, r] = 1.0 / mu_r  # response time = service time (infinite server)
            U[delay_idx, r] = T[delay_idx, r] / mu_r  # "utilization" for delay

    result = NCResult()
    result.QN = Q
    result.UN = U
    result.RN = R
    result.TN = T
    result.XN = qlen.reshape(1, -1)  # system throughput per class
    result.CN = np.zeros((1, K))  # cycle time not applicable for open networks
    result.lG = np.nan  # no normalizing constant for loss networks
    result.it = niter
    result.method = 'erlangfp'
    result.runtime = time.time() - start_time
    result.solver = "NC"

    return result


__all__ = [
    'NCResult',
    'NCResultProb',
    'solver_nc_analyzer',
    'solver_ncld_analyzer',
    'solver_nc_lossn_analyzer',
]
