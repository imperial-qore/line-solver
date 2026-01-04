"""
MVA Solver handler.

Native Python implementation of MVA solver handler that orchestrates
chain aggregation, core MVA computation, and result disaggregation.

Port from:

"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple, List
import time

from ...sn import (
    NetworkStruct,
    SchedStrategy,
    NodeType,
    sn_get_demands_chain,
    sn_deaggregate_chain_results,
    sn_has_product_form,
)
from ...pfqn.mva import pfqn_mva
from ...pfqn.mvald import pfqn_mvams


def _sched_matches(sched_value, *strategies) -> bool:
    """Check if a scheduling strategy matches any of the given strategies.

    Handles comparison across different enum class implementations by comparing
    enum names as strings.
    """
    if sched_value is None:
        return False

    # Get the name of the scheduling strategy
    if hasattr(sched_value, 'name'):
        sched_name = sched_value.name
    else:
        sched_name = str(sched_value)

    # Compare against each target strategy
    for strategy in strategies:
        if hasattr(strategy, 'name'):
            target_name = strategy.name
        else:
            target_name = str(strategy)
        if sched_name == target_name:
            return True

    return False


@dataclass
class SolverMVAOptions:
    """Options for MVA solver."""
    method: str = 'exact'
    tol: float = 1e-8
    verbose: bool = False


@dataclass
class SolverMVAReturn:
    """
    Result of MVA solver handler.

    Attributes:
        Q: Mean queue lengths (M x K)
        U: Utilizations (M x K)
        R: Response times (M x K)
        T: Throughputs (M x K)
        C: Cycle times (1 x K)
        X: System throughputs (1 x K)
        lG: Log normalizing constant
        runtime: Runtime in seconds
        method: Method used
        it: Number of iterations
    """
    Q: Optional[np.ndarray] = None
    U: Optional[np.ndarray] = None
    R: Optional[np.ndarray] = None
    T: Optional[np.ndarray] = None
    C: Optional[np.ndarray] = None
    X: Optional[np.ndarray] = None
    lG: float = 0.0
    runtime: float = 0.0
    method: str = "exact"
    it: int = 0


def solver_mva(
    sn: NetworkStruct,
    options: Optional[SolverMVAOptions] = None
) -> SolverMVAReturn:
    """
    MVA solver handler.

    Performs Mean Value Analysis by:
    1. Aggregating class-level parameters into chains
    2. Separating delay stations from queueing stations
    3. Computing chain-level performance using pfqn_mvams
    4. Disaggregating results back to class level

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverMVAReturn with all performance metrics

    Raises:
        RuntimeError: For unsupported configurations (non-product-form,
                     LCFS without LCFS-PR, etc.)
    """
    start_time = time.time()

    if options is None:
        options = SolverMVAOptions()

    M = sn.nstations
    K = sn.nclasses

    # Check for LCFS scheduling - not supported in this version
    sched_dict = sn.sched if sn.sched else {}
    for i in range(M):
        station_sched = sched_dict.get(i)
        if _sched_matches(station_sched, SchedStrategy.LCFS, SchedStrategy.LCFSPR):
            raise RuntimeError("LCFS queueing networks are not supported in this version.")

    # For non-LCFS models, check product-form requirement
    if not sn_has_product_form(sn):
        raise RuntimeError(
            "Unsupported exact MVA analysis, the model does not have a product form"
        )

    # Separate infinite server (delay) stations from queueing stations
    infSET = []  # Delay stations
    qSET = []    # Queueing stations

    for i in range(M):
        station_sched = sched_dict.get(i)

        if _sched_matches(station_sched, SchedStrategy.EXT):
            continue
        elif _sched_matches(station_sched, SchedStrategy.INF):
            infSET.append(i)
        elif _sched_matches(station_sched, SchedStrategy.PS, SchedStrategy.LCFSPR,
                           SchedStrategy.FCFS, SchedStrategy.SIRO):
            qSET.append(i)
        else:
            sched_name = str(station_sched) if station_sched else "unknown"
            raise RuntimeError(f"Unsupported exact MVA analysis for {sched_name} scheduling")

    # Get chain-level demands
    chain_result = sn_get_demands_chain(sn)
    Lchain = chain_result.Lchain
    STchain = chain_result.STchain
    Vchain = chain_result.Vchain
    alpha = chain_result.alpha
    Nchain = chain_result.Nchain.flatten()
    refstatchain = chain_result.refstatchain

    nservers = sn.nservers
    if nservers is None:
        nservers = np.ones(M)
    else:
        nservers = nservers.flatten()

    C_chains = sn.nchains

    # Initialize chain-level result matrices
    Uchain = np.zeros((M, C_chains))
    Tchain = np.zeros((M, C_chains))
    Wchain = np.zeros((M, C_chains))
    Qchain = np.zeros((M, C_chains))
    Xchain = np.zeros((1, C_chains))
    C_result = np.zeros((1, C_chains))
    lG = 0.0

    # Handle open classes
    ocl = [i for i in range(len(Nchain)) if np.isinf(Nchain[i])]
    lambda_chain = np.zeros((1, C_chains))

    for r in ocl:
        ref_stat = int(refstatchain[r, 0]) if refstatchain.ndim > 1 else int(refstatchain[r])
        if STchain[ref_stat, r] > 0:
            lambda_chain[0, r] = 1.0 / STchain[ref_stat, r]
        Qchain[ref_stat, r] = np.inf

    # Classes with non-zero population
    rset = [i for i in range(C_chains) if Nchain[i] != 0.0]

    # Build Lp and Zp for queueing and delay stations
    Lp = np.zeros((len(qSET), C_chains))
    for idx, q_idx in enumerate(qSET):
        for j in range(C_chains):
            Lp[idx, j] = STchain[q_idx, j] * Vchain[q_idx, j]

    Zp = np.zeros((len(infSET), C_chains))
    for idx, inf_idx in enumerate(infSET):
        for j in range(C_chains):
            Zp[idx, j] = STchain[inf_idx, j] * Vchain[inf_idx, j]

    # Get nservers for queueing stations
    nserversp = np.array([nservers[q] for q in qSET]).reshape(-1, 1) if qSET else np.ones((1, 1))

    # Compute Z (think times) as sum of delay demands
    Z_total = Zp.sum(axis=0) if Zp.size > 0 else np.zeros(C_chains)

    # Call MVA core algorithm
    if len(qSET) > 0:
        Xchain_out, Qpf, Upf, Cpf, lG = pfqn_mvams(
            lambda_chain.flatten(), Lp, Nchain, Z_total,
            mi=np.ones(len(qSET)),
            S=nserversp.flatten().astype(int)
        )

        # Map results back to full station indices
        for idx, q_idx in enumerate(qSET):
            for j in range(C_chains):
                Qchain[q_idx, j] = Qpf[idx, j] if Qpf.ndim > 1 else Qpf[idx]

        Xchain = Xchain_out.reshape(1, -1) if Xchain_out.ndim == 1 else Xchain_out
    else:
        # No queueing stations - compute throughput from delays only
        for r in range(C_chains):
            if np.isfinite(Nchain[r]) and Z_total[r] > 0:
                Xchain[0, r] = Nchain[r] / Z_total[r]

    # Compute queue lengths at delay stations
    for idx, inf_idx in enumerate(infSET):
        for j in range(C_chains):
            Qchain[inf_idx, j] = Xchain[0, j] * STchain[inf_idx, j] * Vchain[inf_idx, j]

    # Compute waiting times
    ccl = [i for i in range(len(Nchain)) if np.isfinite(Nchain[i])]
    for r in rset:
        for k in infSET:
            Wchain[k, r] = STchain[k, r]
        for k in qSET:
            if np.isinf(nservers[k]):
                Wchain[k, r] = STchain[k, r]
            else:
                if Vchain[k, r] == 0.0 or Xchain[0, r] == 0.0:
                    Wchain[k, r] = 0.0
                else:
                    Wchain[k, r] = Qchain[k, r] / (Xchain[0, r] * Vchain[k, r])

    # Compute response times and throughputs
    for r in rset:
        W_sum = np.sum(Vchain[:, r] * Wchain[:, r])
        if W_sum == 0.0:
            Xchain[0, r] = 0.0
        else:
            if np.isinf(Nchain[r]):
                C_result[0, r] = W_sum
            elif Nchain[r] == 0.0:
                Xchain[0, r] = 0.0
                C_result[0, r] = 0.0
            else:
                C_result[0, r] = W_sum
                Xchain[0, r] = Nchain[r] / C_result[0, r]

        for k in range(M):
            Qchain[k, r] = Xchain[0, r] * Vchain[k, r] * Wchain[k, r]
            Tchain[k, r] = Xchain[0, r] * Vchain[k, r]

    # Compute utilizations
    for k in range(M):
        for r in rset:
            if np.isinf(nservers[k]):
                Uchain[k, r] = Vchain[k, r] * STchain[k, r] * Xchain[0, r]
            else:
                Uchain[k, r] = Vchain[k, r] * STchain[k, r] * Xchain[0, r] / nservers[k]

    # Utilization capping for FCFS/PS stations
    for k in range(M):
        station_sched = sched_dict.get(k)
        for r in range(C_chains):
            if Vchain[k, r] * STchain[k, r] > options.tol:
                if _sched_matches(station_sched, SchedStrategy.FCFS, SchedStrategy.PS):
                    Urow_sum = np.sum(Uchain[k, :])
                    if Urow_sum > 1 + options.tol:
                        denom = np.sum(Vchain[k, :] * STchain[k, :] * Xchain[0, :])
                        if denom > 0:
                            Uchain[k, r] = (min(1.0, Urow_sum) * Vchain[k, r] *
                                           STchain[k, r] * Xchain[0, r] / denom)

    # Compute Rchain from Qchain and Tchain
    Rchain = np.zeros((M, C_chains))
    for i in range(M):
        for j in range(C_chains):
            if Tchain[i, j] > 0:
                Rchain[i, j] = Qchain[i, j] / Tchain[i, j]

    # Clean up non-finite values
    Xchain = np.where(~np.isfinite(Xchain), 0.0, Xchain)
    Uchain = np.where(~np.isfinite(Uchain), 0.0, Uchain)
    Qchain = np.where(~np.isfinite(Qchain), 0.0, Qchain)
    Rchain = np.where(~np.isfinite(Rchain), 0.0, Rchain)

    # Zero out results for zero population chains
    Nzero = [i for i in range(C_chains) if Nchain[i] == 0.0]
    for j in Nzero:
        Xchain[0, j] = 0.0
        Uchain[:, j] = 0.0
        Qchain[:, j] = 0.0
        Rchain[:, j] = 0.0
        Tchain[:, j] = 0.0
        Wchain[:, j] = 0.0

    # Disaggregate chain results to class level
    deagg_result = sn_deaggregate_chain_results(
        sn, Lchain, None, STchain, Vchain, alpha,
        None, None, Rchain, Tchain, None, Xchain
    )

    result = SolverMVAReturn(
        Q=deagg_result.Q,
        U=deagg_result.U,
        R=deagg_result.R,
        T=deagg_result.T,
        C=deagg_result.C,
        X=deagg_result.X,
        lG=lG,
        runtime=time.time() - start_time,
        method=options.method,
        it=0
    )

    return result


__all__ = [
    'solver_mva',
    'SolverMVAReturn',
    'SolverMVAOptions',
]
