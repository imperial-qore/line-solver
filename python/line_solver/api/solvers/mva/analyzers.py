"""
MVA Solver analyzers.

Native Python implementation of MVA solver analyzers that orchestrate
method selection and provide the main entry point for MVA analysis.

Port from:

"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional
import time

from ...sn import (
    NetworkStruct,
    SchedStrategy,
    sn_has_product_form,
    sn_has_fractional_populations,
)
from .handler import solver_mva, SolverMVAOptions, SolverMVAReturn


@dataclass
class MVAResult:
    """
    Result of MVA solver analysis.

    Attributes:
        QN: Mean queue lengths (M x K)
        UN: Utilizations (M x K)
        RN: Response times (M x K)
        TN: Throughputs (M x K)
        CN: Cycle times (1 x K)
        XN: System throughputs (1 x K)
        AN: Arrival rates (M x K)
        WN: Waiting times (M x K)
        logNormConstAggr: Log normalizing constant
        iter: Number of iterations
        runtime: Runtime in seconds
        method: Method used
    """
    QN: Optional[np.ndarray] = None
    UN: Optional[np.ndarray] = None
    RN: Optional[np.ndarray] = None
    TN: Optional[np.ndarray] = None
    CN: Optional[np.ndarray] = None
    XN: Optional[np.ndarray] = None
    AN: Optional[np.ndarray] = None
    WN: Optional[np.ndarray] = None
    logNormConstAggr: float = np.nan
    iter: int = 0
    runtime: float = 0.0
    method: str = ""


def _has_lcfs_lcfspr_config(sn: NetworkStruct) -> bool:
    """
    Check if network has special LCFS/LCFS-PR configuration.

    This configuration has a product-form solution but is not
    recognized by sn_has_product_form.

    Args:
        sn: Network structure

    Returns:
        True if network has both LCFS and LCFS-PR stations
    """
    has_lcfs = False
    has_lcfspr = False

    sched_dict = sn.sched if sn.sched else {}
    for i in range(sn.nstations):
        station_sched = sched_dict.get(i)

        if station_sched == SchedStrategy.LCFS:
            has_lcfs = True
        elif station_sched == SchedStrategy.LCFSPR:
            has_lcfspr = True

    return has_lcfs and has_lcfspr


def solver_mva_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverMVAOptions] = None
) -> MVAResult:
    """
    MVA Analyzer - main entry point for MVA analysis.

    Selects appropriate MVA method based on network characteristics
    and solver options, then performs the analysis.

    Supported methods:
        - 'exact', 'mva': Exact MVA using population recursion
        - 'default': Automatic method selection
        - 'amva', 'bs', 'qd', 'sqni', etc.: Approximate MVA methods
        - 'qna': Queueing Network Analyzer (for general networks)

    Args:
        sn: Network structure
        options: Solver options (method, tolerance, verbosity)

    Returns:
        MVAResult with all performance metrics

    Raises:
        RuntimeError: For unsupported methods or configurations
    """
    start_time = time.time()

    if options is None:
        options = SolverMVAOptions()

    method = options.method.lower()

    # Remove 'amva.' prefix if present
    if method.startswith('amva.'):
        method = method[5:]

    result = MVAResult()
    ret = None

    if method in ['exact', 'mva']:
        ret = solver_mva(sn, options)
        result.iter = 0

    elif method == 'qna':
        # QNA for general networks
        ret = solver_qna(sn, options)

    elif method == 'default':
        # Automatic method selection
        if _has_lcfs_lcfspr_config(sn):
            # LCFS/LCFS-PR configuration has product-form solution
            ret = solver_mva(sn, options)
            method = 'exact'
        elif (sn.nchains <= 4 and
              np.sum(sn.njobs) <= 20 and
              sn_has_product_form(sn) and
              not sn_has_fractional_populations(sn)):
            # Small network with product-form - use exact MVA
            ret = solver_mva(sn, options)
            method = 'exact'
        else:
            # Use approximate MVA
            ret = solver_amva(sn, options)
            method = ret.method if hasattr(ret, 'method') else 'amva'

    elif method in ['amva', 'bs', 'qd', 'qli', 'fli', 'lin',
                    'qdlin', 'sqni', 'gflin', 'egflin']:
        ret = solver_amva(sn, options)
        method = ret.method if hasattr(ret, 'method') else method

    else:
        # Unsupported method
        result.QN = None
        result.UN = None
        result.RN = None
        result.TN = None
        result.CN = None
        result.XN = None
        result.AN = None
        result.WN = None
        result.logNormConstAggr = np.nan
        result.iter = 0
        result.runtime = time.time() - start_time
        result.method = method
        if options.verbose:
            print(f"Warning: Unsupported SolverMVA method: {method}")
        return result

    # Copy results from handler return
    if ret is not None:
        result.QN = ret.Q
        result.UN = ret.U
        result.RN = ret.R
        result.TN = ret.T
        result.CN = ret.C
        result.XN = ret.X
        result.AN = ret.T.copy() if ret.T is not None else None  # Arrival rates = throughputs
        result.WN = ret.R.copy() if ret.R is not None else None  # Waiting times ≈ response times
        result.logNormConstAggr = ret.lG
        result.iter = ret.it

    result.runtime = time.time() - start_time
    result.method = method

    return result


def solver_amva(
    sn: NetworkStruct,
    options: Optional[SolverMVAOptions] = None
) -> SolverMVAReturn:
    """
    Approximate MVA solver.

    Uses approximate MVA methods for larger networks where exact
    MVA would be computationally expensive.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverMVAReturn with performance metrics
    """
    from ...pfqn.mva import pfqn_aql, pfqn_bs, pfqn_sqni
    from .amvald import solver_amvald, AmvaldOptions

    start_time = time.time()

    if options is None:
        options = SolverMVAOptions()

    method = options.method.lower()
    if method.startswith('amva.'):
        method = method[5:]

    M = sn.nstations
    K = sn.nclasses

    # Get chain-level demands
    from ...sn import sn_get_demands_chain, sn_deaggregate_chain_results
    chain_result = sn_get_demands_chain(sn)
    Lchain = chain_result.Lchain
    STchain = chain_result.STchain
    Vchain = chain_result.Vchain
    alpha = chain_result.alpha
    Nchain_float = chain_result.Nchain.flatten()

    # Check if this is a mixed model (has both open and closed chains)
    has_open = np.any(np.isinf(Nchain_float))
    has_closed = np.any(np.isfinite(Nchain_float) & (Nchain_float > 0))
    is_mixed = has_open and has_closed

    # Use solver_amvald for mixed models, pure open models, or advanced methods (lin, qdlin, egflin, etc.)
    # Pure open networks (has_open and not has_closed) must also use solver_amvald
    # to match MATLAB behavior which falls through to solver_amvald at line 268
    if has_open or method in ('lin', 'qdlin', 'fli', 'gflin', 'egflin', 'qd', 'qli'):
        # Use solver_amvald which properly handles mixed models with class switching
        amvald_options = AmvaldOptions(
            method=method if method in ('lin', 'qdlin', 'fli', 'gflin', 'egflin', 'qd', 'qli') else 'egflin',
            iter_tol=options.iter_tol if hasattr(options, 'iter_tol') and options.iter_tol else 1e-4,  # Match MATLAB lineDefaults
            iter_max=options.iter_max if hasattr(options, 'iter_max') and options.iter_max else 100,  # Match MATLAB lineDefaults
        )

        # Get SCVchain and refstatchain
        SCVchain = chain_result.SCVchain if hasattr(chain_result, 'SCVchain') and chain_result.SCVchain is not None else np.ones((M, sn.nchains))
        refstatchain = chain_result.refstatchain if hasattr(chain_result, 'refstatchain') and chain_result.refstatchain is not None else np.zeros((sn.nchains, 1))

        amvald_result = solver_amvald(
            sn, Lchain, STchain, Vchain, alpha,
            Nchain_float, SCVchain, refstatchain, amvald_options
        )

        # Disaggregate to class level
        Xchain = amvald_result.X.reshape(1, -1) if amvald_result.X.ndim == 1 else amvald_result.X
        deagg_result = sn_deaggregate_chain_results(
            sn, Lchain, amvald_result.Q, STchain, Vchain, alpha,
            amvald_result.U, None, amvald_result.R, amvald_result.T, None, Xchain
        )

        result = SolverMVAReturn(
            Q=deagg_result.Q,
            U=deagg_result.U,
            R=deagg_result.R,
            T=deagg_result.T,
            C=deagg_result.C,
            X=deagg_result.X,
            lG=amvald_result.lG,
            runtime=time.time() - start_time,
            method=amvald_options.method,
            it=amvald_result.totiter
        )
        return result

    # For closed-only networks, use simpler methods
    Nchain_float = np.where(np.isfinite(Nchain_float), Nchain_float, 0)
    Nchain = Nchain_float.astype(int)

    # Get think times from delay stations
    sched_dict = sn.sched if sn.sched else {}
    Z = np.zeros(sn.nchains)
    for i in range(M):
        station_sched = sched_dict.get(i)
        if station_sched == SchedStrategy.INF:
            for c in range(sn.nchains):
                Z[c] += Lchain[i, c]

    # Get demands at queueing stations
    L_queue = []
    for i in range(M):
        station_sched = sched_dict.get(i)
        if station_sched not in [SchedStrategy.INF, SchedStrategy.EXT]:
            L_queue.append(Lchain[i, :])

    if len(L_queue) == 0:
        # All delay stations
        L = np.zeros((1, sn.nchains))
    else:
        L = np.array(L_queue)

    # Choose algorithm
    if method == 'bs':
        XN, QN, UN, RN, _ = pfqn_bs(L, Nchain, Z)
        TN = np.tile(XN, (QN.shape[0], 1))
        AN = TN.copy()
    elif method == 'sqni' and L.shape[0] == 1:
        Q, U, X = pfqn_sqni(L.flatten(), Nchain, Z)
        XN = X
        QN = Q[:1, :]
        UN = U[:1, :]
        CN = np.zeros((1, sn.nchains))
        RN = np.zeros_like(QN)
        TN = np.tile(XN, (QN.shape[0], 1))
        AN = TN.copy()
    else:
        # Default to AQL (Schweitzer approximation)
        XN, CN, QN, UN, RN, TN, AN = pfqn_aql(L, Nchain, Z)

    # Disaggregate to class level
    Rchain = np.zeros((M, sn.nchains))
    Tchain = np.zeros((M, sn.nchains))
    q_idx = 0
    for i in range(M):
        station_sched = sched_dict.get(i)
        if station_sched not in [SchedStrategy.INF, SchedStrategy.EXT]:
            if q_idx < RN.shape[0]:
                Rchain[i, :] = RN[q_idx, :]
                Tchain[i, :] = TN[q_idx, :]
            q_idx += 1
        elif station_sched == SchedStrategy.INF:
            Tchain[i, :] = XN.flatten()
            Rchain[i, :] = STchain[i, :] * Vchain[i, :]

    Xchain = XN.reshape(1, -1) if XN.ndim == 1 else XN

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
        lG=np.nan,
        runtime=time.time() - start_time,
        method=method if method != 'amva' else 'aql',
        it=0
    )

    return result


def solver_qna(
    sn: NetworkStruct,
    options: Optional[SolverMVAOptions] = None
) -> SolverMVAReturn:
    """
    Queueing Network Analyzer.

    Provides approximate analysis for general queueing networks
    including those with non-product-form characteristics.

    Implementation based on N. Gautaum's "Analysis of Queues" (CRC Press, 2012),
    Section 7.2.3 with minor corrections.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverMVAReturn with performance metrics

    References:
        MATLAB: matlab/src/solvers/MVA/solver_qna.m
    """
    start_time = time.time()

    if options is None:
        options = SolverMVAOptions()

    K = sn.nclasses
    M = sn.nstations
    C = sn.nchains

    # Extract parameters from network
    rt = sn.rt if sn.rt is not None else np.zeros((M * K, M * K))
    S = 1.0 / (sn.rates + 1e-10)  # Service times (M, K)
    scv = sn.scv.copy() if sn.scv is not None else np.ones((M, K))
    scv[np.isnan(scv)] = 0

    # Initialize results
    Q = np.zeros((M, K))
    Q_prev = Q + np.inf
    U = np.zeros((M, K))
    R = np.zeros((M, K))
    T = np.zeros((M, K))
    X = np.zeros((1, K))

    # Compute aggregate visit ratios
    V = np.zeros((M, K))
    if sn.visits:
        for c in sn.visits:
            V += sn.visits[c]

    # Chain information
    lambda_chain = np.zeros(C)
    inchain = sn.inchain if sn.inchain else {}

    tol = options.iter_tol if hasattr(options, 'iter_tol') else 1e-6
    max_iter = options.iter_max if hasattr(options, 'iter_max') else 1000

    # Compute departure process at source
    a1 = np.zeros((M, K))  # Arrival rates
    a2 = np.zeros((M, K))  # SCVs of arrivals
    d2 = np.zeros(M)  # SCVs of departure processes
    f2 = np.ones((M * K, M * K))  # SCV of each flow pair

    # Initialize throughputs at source
    for c in range(C):
        if c in inchain:
            classes_c = inchain[c]
            refstat = int(sn.refstat.flatten()[classes_c[0]])
            if refstat < M:
                for k in classes_c:
                    k = int(k)
                    if k < K and sn.njobs[k] == np.inf:  # Open class
                        T[refstat, k] = sn.rates[refstat, k]
                        lambda_chain[c] += sn.rates[refstat, k]

    # Compute initial SCVs at source
    for c in range(C):
        if c in inchain:
            classes_c = inchain[c]
            refstat = int(sn.refstat.flatten()[classes_c[0]])
            if refstat < M and lambda_chain[c] > 0:
                d2_c = 0
                for k in classes_c:
                    k = int(k)
                    if k < K and sn.njobs[k] == np.inf:
                        d2_c += scv[refstat, k] * sn.rates[refstat, k]
                d2[refstat] = d2_c / lambda_chain[c]

    # Main QNA iteration loop
    iteration = 0
    while np.max(np.abs(Q - Q_prev)) > tol and iteration < max_iter:
        iteration += 1
        Q_prev = Q.copy()

        # Normalize queue lengths for closed classes
        for c in range(C):
            if c in inchain:
                classes_c = inchain[c]
                njobs_c = np.sum([sn.njobs[int(k)] for k in classes_c if int(k) < K])
                if np.isfinite(njobs_c) and njobs_c > 0:
                    Q_sum = np.sum(Q[:, [int(k) for k in classes_c if int(k) < K]])
                    if Q_sum > 0:
                        Q[:, [int(k) for k in classes_c if int(k) < K]] *= njobs_c / Q_sum

        # Update throughputs
        if iteration == 1:
            for c in range(C):
                if c in inchain:
                    classes_c = inchain[c]
                    if lambda_chain[c] > 0:
                        for m in range(M):
                            for k in classes_c:
                                k = int(k)
                                if k < K:
                                    T[m, k] = V[m, k] * lambda_chain[c]

        # Superposition: compute arrival process parameters at each station
        for ist in range(M):
            a1[ist, :] = 0
            a2[ist, :] = 0
            lambda_i = np.sum(T[ist, :])

            for jst in range(M):
                for r in range(K):
                    for s in range(K):
                        idx_from = (jst) * K + r
                        idx_to = (ist) * K + s
                        if idx_from < rt.shape[0] and idx_to < rt.shape[1] and rt[idx_from, idx_to] > 0:
                            a1[ist, s] += T[jst, r] * rt[idx_from, idx_to]
                            if lambda_i > 0:
                                a2[ist, s] += (1.0 / lambda_i) * f2[idx_from, idx_to] * T[jst, r] * rt[idx_from, idx_to]

        # Handle different scheduling strategies
        sched_dict = sn.sched if sn.sched else {}
        for ist in range(M):
            station_sched = sched_dict.get(ist, SchedStrategy.FCFS)

            if station_sched == SchedStrategy.INF:  # Delay station
                for k in range(K):
                    T[ist, k] = a1[ist, k]
                    Q[ist, k] = T[ist, k] * S[ist, k] * V[ist, k]
                    U[ist, k] = Q[ist, k]
                    R[ist, k] = S[ist, k] * V[ist, k]
                d2[ist] = np.sum(a2[ist, :]) / np.sum(a1[ist, :]) if np.sum(a1[ist, :]) > 0 else 1

            elif station_sched == SchedStrategy.PS:  # Processor Sharing
                for c in range(C):
                    if c in inchain:
                        classes_c = inchain[c]
                        for k in classes_c:
                            k = int(k)
                            if k < K and lambda_chain[c] > 0:
                                T[ist, k] = lambda_chain[c] * V[ist, k]
                                U[ist, k] = S[ist, k] * T[ist, k]

                        # Compute queue length using approximation
                        Nc = np.sum([sn.njobs[int(k)] for k in classes_c if int(k) < K and np.isfinite(sn.njobs[int(k)])])
                        Uden = min(1.0 - options.tol if hasattr(options, 'tol') else 1e-6, np.sum(U[ist, :]))

                        for k in classes_c:
                            k = int(k)
                            if k < K:
                                if Uden < 1.0:
                                    Q[ist, k] = (U[ist, k] - U[ist, k] ** (Nc + 1)) / (1.0 - Uden)
                                else:
                                    Q[ist, k] = sn.njobs[k] if np.isfinite(sn.njobs[k]) else 1
                                R[ist, k] = Q[ist, k] / (T[ist, k] + 1e-10)

            elif station_sched == SchedStrategy.FCFS:  # FCFS queue
                mu_ist = sn.rates[ist, :]
                rho_ist_class = a1[ist, :] / (1e-10 + mu_ist)
                lambda_ist = np.sum(a1[ist, :])
                mi = sn.nservers[ist] if ist < len(sn.nservers) else 1
                rho_ist = np.sum(rho_ist_class) / mi

                if rho_ist < 1.0 - (options.tol if hasattr(options, 'tol') else 1e-6):
                    # Compute waiting time using diffusion approximation
                    mubar = lambda_ist / (rho_ist + 1e-10)
                    alpha_mi = (rho_ist ** mi + rho_ist) / 2.0 if rho_ist > 0.7 else rho_ist ** ((mi + 1) / 2.0)
                    c2 = -1

                    for r in range(K):
                        if mu_ist[r] > 0:
                            c2 += (a1[ist, r] / (lambda_ist + 1e-10)) * (mubar / (mi * mu_ist[r])) ** 2 * (scv[ist, r] + 1)

                    Wiq = (alpha_mi / mubar) * (1.0 / (1.0 - rho_ist + 1e-10)) * (np.sum(a2[ist, :]) + c2) / 2.0

                    for k in range(K):
                        Q[ist, k] = a1[ist, k] / (mu_ist[k] + 1e-10) + a1[ist, k] * Wiq
                        T[ist, k] = a1[ist, k]
                        U[ist, k] = T[ist, k] * S[ist, k] / mi
                        R[ist, k] = Q[ist, k] / (T[ist, k] + 1e-10)

                    d2[ist] = 1 + rho_ist ** 2 * (c2 - 1) / np.sqrt(mi) + (1 - rho_ist ** 2) * (np.sum(a2[ist, :]) - 1)
                else:
                    # System overloaded
                    for k in range(K):
                        Q[ist, k] = sn.njobs[k] if np.isfinite(sn.njobs[k]) else 1
                    d2[ist] = 1

        # Update flow SCVs for splitting
        for ist in range(M):
            for jst in range(M):
                for r in range(K):
                    for s in range(K):
                        idx_from = ist * K + r
                        idx_to = jst * K + s
                        if idx_from < rt.shape[0] and idx_to < rt.shape[1] and rt[idx_from, idx_to] > 0:
                            f2[idx_from, idx_to] = 1 + rt[idx_from, idx_to] * (d2[ist] - 1)

    # Final cleanup and normalization
    Q = np.abs(Q)
    Q[np.isnan(Q)] = 0
    U[np.isnan(U)] = 0
    R[np.isnan(R)] = 0
    T[np.isnan(T)] = 0

    C_result = np.sum(R, axis=0, keepdims=True)

    result = SolverMVAReturn(
        Q=Q,
        U=U,
        R=R,
        T=T,
        C=C_result,
        X=X,
        lG=0.0,  # QNA does not compute normalizing constant
        runtime=time.time() - start_time,
        method='qna',
        it=iteration
    )

    return result


__all__ = [
    'MVAResult',
    'solver_mva_analyzer',
    'solver_amva',
    'solver_qna',
]
