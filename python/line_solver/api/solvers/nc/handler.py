"""
NC Solver main handler.

Native Python implementation of the main NC solver handler that orchestrates
normalizing constant computation for product-form queueing networks.

Port from:

"""

import numpy as np
from math import exp, log
from typing import Optional, Dict, Any
from dataclasses import dataclass
import time

from ...sn import (
    NetworkStruct,
    SchedStrategy,
    sn_get_demands_chain,
    sn_deaggregate_chain_results,
)
from ...pfqn import (
    pfqn_nc,
)


# Constants
FINE_TOL = 1e-12
ZERO = 1e-10


@dataclass
class SolverNCReturn:
    """Result of NC solver analysis."""
    Q: Optional[np.ndarray]  # Mean queue lengths (M x K)
    U: Optional[np.ndarray]  # Utilizations (M x K)
    R: Optional[np.ndarray]  # Response times (M x K)
    T: Optional[np.ndarray]  # Throughputs (M x K)
    nchains: int             # Number of chains
    X: Optional[np.ndarray]  # System throughputs (1 x K)
    lG: float                # Log normalizing constant
    STeff: Optional[np.ndarray]  # Effective service times
    it: int                  # Number of iterations
    runtime: float           # Runtime in seconds
    method: str              # Method used


@dataclass
class SolverOptions:
    """Solver options for NC analysis."""
    method: str = 'default'
    tol: float = 1e-6
    iter_max: int = 1000
    iter_tol: float = 1e-4
    verbose: int = 0
    highvar: str = 'interp'  # Non-exponential approximation method ('interp', 'default', 'none')


def solver_nc(sn: NetworkStruct, options: Optional[SolverOptions] = None) -> SolverNCReturn:
    """
    Main NC solver handler.

    Performs normalizing constant analysis of a product-form queueing network.
    Supports standard closed/mixed networks, multiserver stations, and
    non-exponential service time distributions via approximations.

    Args:
        sn: NetworkStruct describing the queueing network
        options: Solver options

    Returns:
        SolverNCReturn with performance metrics

    Raises:
        RuntimeError: For unsupported network configurations
    """
    if options is None:
        options = SolverOptions()

    start_time = time.time()

    M = sn.nstations
    K = sn.nclasses
    C = sn.nchains

    nservers = sn.nservers
    if nservers is None:
        nservers = np.ones((M, 1))
    nservers = nservers.flatten()

    NK = sn.njobs.T if sn.njobs is not None else np.zeros((K, 1))
    if NK.ndim == 1:
        NK = NK.reshape(-1, 1)

    sched = sn.sched
    SCV = sn.scv

    # Check for LCFS scheduling - not supported in this version
    if sched is not None:
        for i in range(M):
            station_key = i
            if station_key in sched:
                if sched[station_key] in (SchedStrategy.LCFS, SchedStrategy.LCFSPR):
                    raise RuntimeError("LCFS queueing networks are not supported in this version.")

    # Compute visits matrix
    V = np.zeros((M, K))
    if sn.visits is not None:
        for key, visit_matrix in sn.visits.items():
            if visit_matrix is not None:
                if visit_matrix.shape == V.shape:
                    V = V + visit_matrix
                elif visit_matrix.shape[0] == V.shape[0]:
                    # Handle different column counts
                    min_cols = min(visit_matrix.shape[1], V.shape[1])
                    V[:, :min_cols] += visit_matrix[:, :min_cols]

    # Compute service times from rates
    rates = sn.rates
    if rates is None:
        rates = np.ones((M, K))
    with np.errstate(divide='ignore', invalid='ignore'):
        ST = np.where(rates > 0, 1.0 / rates, 0.0)
        ST = np.where(np.isnan(ST), 0.0, ST)
    ST0 = ST.copy()

    # Compute chain populations
    Nchain = np.zeros(C)
    for c in range(C):
        if c not in sn.inchain:
            continue
        inchain = sn.inchain[c].flatten().astype(int)
        Nchain[c] = sum(NK[k, 0] for k in inchain if k < NK.shape[0])

    # Identify open and closed chains
    open_chains = [c for c in range(C) if np.isinf(Nchain[c])]
    closed_chains = [c for c in range(C) if np.isfinite(Nchain[c])]

    # Check if iteration is needed (FCFS with non-exponential)
    has_fcfs = False
    if sched is not None:
        for station_sched in sched.values():
            if station_sched == SchedStrategy.FCFS:
                has_fcfs = True
                break

    if not has_fcfs:
        options.iter_max = 1

    # Initialize iteration variables
    gamma = np.zeros(M)
    eta_1 = np.zeros(M)
    eta = np.ones(M)
    it = 0
    tmp_eta = np.abs(1 - eta / (eta_1 + FINE_TOL))

    # Initialize output variables
    lmbda = None
    Lchain = None
    STchain = None
    Vchain = None
    alpha = None
    Q = None
    U = None
    R = None
    T = None
    X = None
    STeff = None
    lG = np.nan
    method = options.method

    # Main iteration loop
    while np.max(tmp_eta) > options.iter_tol and it < options.iter_max:
        it += 1
        eta_1 = eta.copy()

        if it == 1:
            # First iteration: compute chain-level demands
            lmbda = np.zeros(C)

            result = sn_get_demands_chain(sn)
            Lchain = result.Lchain
            STchain = result.STchain
            Vchain = result.Vchain
            alpha = result.alpha
            Nchain = result.Nchain.flatten() if result.Nchain is not None else np.zeros(C)

            # Set arrival rates for open chains
            # For open chains, find Source station and get its rate
            source_station = -1
            if sched is not None:
                for ist, sched_val in sched.items():
                    # Check for EXT scheduling (Source station)
                    if sched_val == SchedStrategy.EXT or (isinstance(sched_val, int) and sched_val == 16):
                        source_station = ist
                        break

            for c in range(C):
                if c not in sn.inchain:
                    continue
                inchain = sn.inchain[c].flatten().astype(int)

                is_open_chain = any(np.isinf(NK[k, 0]) for k in inchain if k < NK.shape[0])

                if is_open_chain:
                    # For open chains, arrival rate is 1/STchain at reference station
                    # STchain already correctly aggregates: STchain = 1/sum(arrival_rates_in_chain)
                    # This handles chains with multiple open classes correctly
                    first_class = inchain[0] if len(inchain) > 0 else 0
                    if sn.refstat is not None and first_class < len(sn.refstat):
                        ref_idx = int(sn.refstat[first_class])
                        if ref_idx < M and STchain[ref_idx, c] > 0:
                            lmbda[c] = 1.0 / STchain[ref_idx, c]
        else:
            # Subsequent iterations: update chain-level service times
            for c in range(C):
                if c not in sn.inchain:
                    continue
                inchain = sn.inchain[c].flatten().astype(int)

                for i in range(M):
                    ST_tmp = np.array([ST[i, k] for k in inchain if k < K])
                    alpha_tmp = np.array([alpha[i, k] for k in inchain if k < K])
                    if len(ST_tmp) > 0:
                        STchain[i, c] = np.dot(ST_tmp, alpha_tmp)
                        Lchain[i, c] = Vchain[i, c] * STchain[i, c]

        # Remove infinities
        STchain = np.where(np.isinf(STchain), 0.0, STchain)
        Lchain = np.where(np.isinf(Lchain), 0.0, Lchain)

        # Separate delay and queueing stations
        Lms = np.zeros((M, C))
        Z = np.zeros((M, C))
        Zms = np.zeros((M, C))

        inf_servers = []
        for i in range(M):
            if np.isinf(nservers[i]):
                inf_servers.append(i)
                Z[i, :] = Lchain[i, :]
            else:
                if options.method == 'exact' and nservers[i] > 1:
                    if options.verbose > 0:
                        print("SolverNC does not support exact multiserver yet. "
                              "Switching to approximate method.")

                for c in range(C):
                    Lms[i, c] = Lchain[i, c] / nservers[i]
                    Zms[i, c] = Lchain[i, c] * (nservers[i] - 1) / nservers[i]

        # Aggregate think times
        Z_new = np.sum(Z, axis=0) + np.sum(Zms, axis=0)

        # For mixed networks: scale demands by remaining capacity after open class load
        # This follows MATLAB's pfqn_nc: L(i,:) = L(i,:) / (1 - lambda*L(i,:)')
        Lms_scaled = Lms.copy()
        Qopen = np.zeros((M, C))
        Ut = np.ones(M)  # Remaining capacity at each station

        if open_chains and lmbda is not None:
            for i in range(M):
                # Compute open class utilization at station i
                open_util = 0.0
                for c in open_chains:
                    if c < len(lmbda):
                        open_util += lmbda[c] * Lchain[i, c] / nservers[i]

                Ut[i] = 1.0 - open_util
                if np.isnan(Ut[i]) or Ut[i] <= 0:
                    Ut[i] = FINE_TOL  # Avoid division by zero

                # Scale demands by remaining capacity
                Lms_scaled[i, :] = Lms[i, :] / Ut[i]

                # Compute open class queue lengths
                for c in open_chains:
                    if c < len(lmbda):
                        Qopen[i, c] = lmbda[c] * Lms[i, c] / Ut[i]

            Qopen = np.where(np.isnan(Qopen), 0.0, Qopen)

        # Also scale think times
        Z_scaled = Z_new.copy()

        # Compute normalizing constant (only for closed chains)
        pfqn_options = {'method': options.method, 'tol': options.tol}
        method = options.method

        if closed_chains:
            # Mixed or closed network: compute normalizing constant for closed chains only
            Nchain_closed = np.array([Nchain[c] if c in closed_chains else 0 for c in range(C)])
            G, lG = pfqn_nc(Lms_scaled, Nchain_closed, Z_scaled, method=options.method)
        else:
            # Purely open network: no normalizing constant needed
            G, lG = 1.0, 0.0

        # Compute throughputs
        Xchain = np.zeros(C)

        if np.sum(Zms) > FINE_TOL:
            # Multi-server case: need more sophisticated computation
            Qchain = np.zeros((M, C))
        else:
            Qchain = np.zeros((M, C))

        # Compute throughputs for closed chains
        for c in closed_chains:
            Nchain_c = float(Nchain[c]) if hasattr(Nchain[c], '__float__') else Nchain[c]
            if Nchain_c <= 0:
                continue

            Nchain_tmp = Nchain.copy()
            Nchain_tmp[c] -= 1
            # Replace infinities with 0 for open classes (NC only handles closed classes)
            Nchain_tmp = np.where(np.isinf(Nchain_tmp), 0.0, Nchain_tmp)

            # Use scaled demands for NC computation
            G_tmp, lG_tmp = pfqn_nc(Lms_scaled, Nchain_tmp, Z_scaled, method=options.method)

            if np.isfinite(lG) and np.isfinite(lG_tmp):
                Xchain[c] = exp(lG_tmp - lG)
            else:
                Xchain[c] = 0.0

            # Compute queue lengths using marginal analysis
            for i in range(M):
                if Lchain[i, c] > FINE_TOL:
                    if np.isinf(nservers[i]):
                        # Delay station: Q = L * X
                        Qchain[i, c] = Lchain[i, c] * Xchain[c]
                    else:
                        # Queueing station: use marginal analysis
                        # Build L_tmp by removing station i and adding it at end with extra class
                        Nchain_aug = np.append(Nchain_tmp, 1.0)
                        Z_aug = np.append(Z_scaled, 0.0)

                        # Create Lms_scaled with station i moved to end and augmented
                        L_tmp = []
                        for j in range(M):
                            if j != i:
                                L_tmp.append(np.append(Lms_scaled[j, :], 0.0))

                        # Add station i at end with demand in augmented class
                        L_i_aug = np.append(Lms_scaled[i, :], 1.0)
                        L_tmp.append(L_i_aug)
                        L_tmp = np.array(L_tmp)

                        # Compute marginal NC
                        _, lG_marg = pfqn_nc(L_tmp, Nchain_aug, Z_aug, method=options.method)

                        if np.isfinite(lG_marg) and np.isfinite(lG):
                            Qchain[i, c] = Zms[i, c] * Xchain[c] + Lms_scaled[i, c] * exp(lG_marg - lG)
                        else:
                            Qchain[i, c] = Zms[i, c] * Xchain[c] + Lms_scaled[i, c] * Xchain[c]

        # Compute queue lengths and throughputs for open chains
        # This follows MATLAB's formula for mixed networks:
        # Qchain(ist,r) = lambda(r)*Lchain(ist,r)/(1-lambda(openChains)*Lchain(ist,openChains)'/nservers(ist))*(1+sum(Qchain(ist,closedChains)))
        for c in open_chains:
            # For open chains, throughput equals arrival rate
            Xchain[c] = lmbda[c] if lmbda is not None and c < len(lmbda) else 0.0

        for c in open_chains:
            for i in range(M):
                if Lchain[i, c] > FINE_TOL:
                    lambda_c = lmbda[c] if lmbda is not None and c < len(lmbda) else 0.0

                    # Sum of queue lengths from closed chains at this station
                    closed_qlen_sum = sum(Qchain[i, cc] for cc in closed_chains)

                    # Mixed network formula from MATLAB:
                    # Qchain(ist,r) = lambda(r)*Lchain(ist,r)/(1-lambda(openChains)*Lchain(ist,openChains)'/nservers(ist))*(1+sum(Qchain(ist,closedChains)))
                    # Ut[i] = 1 - sum(lambda*Lchain/nservers) already incorporates nservers
                    if Ut[i] > FINE_TOL:
                        Qchain[i, c] = (lambda_c * Lchain[i, c] / Ut[i]) * (1 + closed_qlen_sum)
                    else:
                        Qchain[i, c] = np.inf

        # Remove NaN values
        Qchain = np.where(np.isnan(Qchain), 0.0, Qchain)

        # Compute response times
        Rchain = np.zeros((M, C))
        for i in range(M):
            for c in range(C):
                if Qchain[i, c] > ZERO and Xchain[c] > ZERO:
                    Vchain_val = Vchain[i, c] if Vchain[i, c] > 0 else 1.0
                    Rchain[i, c] = Qchain[i, c] / Xchain[c] / Vchain_val
                else:
                    Rchain[i, c] = 0.0

        # Fix response times for delay stations
        for i in inf_servers:
            for c in range(C):
                Vchain_val = Vchain[i, c] if Vchain[i, c] > 0 else 1.0
                Rchain[i, c] = Lchain[i, c] / Vchain_val

        # Compute chain throughputs
        Tchain = Vchain * Xchain

        # Disaggregate to class level
        # Pass None for Qchain to use Rchain-based formula which properly scales
        # by ST/STchain to account for different class service times
        deagg_result = sn_deaggregate_chain_results(
            sn, Lchain, ST, STchain, Vchain, alpha,
            None, None, Rchain, Tchain, None, Xchain.reshape(1, -1)
        )
        Q = deagg_result.Q
        U = deagg_result.U
        R = deagg_result.R
        T = deagg_result.T
        X = deagg_result.X
        STeff = ST.copy()

        # Non-exponential approximation
        from ...npfqn import npfqn_nonexp_approx
        nonexp_result = npfqn_nonexp_approx(
            method=options.highvar,
            sn=sn,
            ST=ST0,  # Use original service times
            V=V,
            SCV=SCV if SCV is not None else np.ones((M, K)),
            Tin=T,
            Uin=U,
            gamma=gamma,
            nservers=nservers
        )
        ST = nonexp_result.ST
        gamma = nonexp_result.gamma
        eta = nonexp_result.eta

        tmp_eta = np.abs(1 - eta / (eta_1 + FINE_TOL))

    # Final cleanup
    if Q is not None:
        Q = np.abs(Q)
        Q = np.where(~np.isfinite(Q), 0.0, Q)
    if R is not None:
        R = np.abs(R)
        R = np.where(~np.isfinite(R), 0.0, R)
    if X is not None:
        X = np.abs(X)
        X = np.where(~np.isfinite(X), 0.0, X)
    if U is not None:
        U = np.abs(U)
        U = np.where(~np.isfinite(U), 0.0, U)
        # Cap utilization at 1.0 for finite-server stations (matches Java/Kotlin Solver_mva.kt)
        # For INF servers, utilization represents mean number of busy servers (not capped)
        for i in range(M):
            if np.isfinite(nservers[i]):
                row_sum = np.sum(U[i, :])
                if row_sum > 1.0:
                    # Rescale to cap at 1.0
                    U[i, :] = U[i, :] / row_sum
    if T is not None:
        T = np.where(~np.isfinite(T), 0.0, T)

    # Chain population corrections
    for c in range(C):
        if c not in sn.inchain:
            continue
        inchain = sn.inchain[c].flatten().astype(int)
        Nchain_c = Nchain[c]

        if np.isfinite(Nchain_c) and Q is not None:
            sum_Q = sum(np.sum(Q[:, k]) for k in inchain if k < K)

            if sum_Q > 0:
                ratio = Nchain_c / sum_Q
                for k in inchain:
                    if k >= K:
                        continue
                    Q[:, k] *= ratio
                    if X is not None:
                        X[0, k] *= ratio
                    if T is not None:
                        T[:, k] *= ratio
                    if U is not None:
                        U[:, k] *= ratio
                    # Recalculate R = Q / T (MATLAB line 270)
                    if R is not None and T is not None:
                        for i in range(M):
                            if T[i, k] > 0:
                                R[i, k] = Q[i, k] / T[i, k]

    runtime = time.time() - start_time

    return SolverNCReturn(
        Q=Q, U=U, R=R, T=T,
        nchains=C, X=X, lG=lG,
        STeff=STeff, it=it, runtime=runtime,
        method=method
    )


@dataclass
class SolverNCLDReturn:
    """Result of load-dependent NC solver analysis."""
    Q: Optional[np.ndarray]  # Mean queue lengths (M x K)
    U: Optional[np.ndarray]  # Utilizations (M x K)
    R: Optional[np.ndarray]  # Response times (M x K)
    T: Optional[np.ndarray]  # Throughputs (M x K)
    C: Optional[np.ndarray]  # Covariances (optional)
    X: Optional[np.ndarray]  # System throughputs (1 x K)
    lG: float                # Log normalizing constant
    runtime: float           # Runtime in seconds
    it: int                  # Number of iterations
    method: str              # Method used


def solver_ncld(sn: NetworkStruct, options: Optional[SolverOptions] = None) -> SolverNCLDReturn:
    """
    Load-dependent NC solver handler.

    Performs normalizing constant analysis for load-dependent product-form
    queueing networks, including multi-server stations with limited load
    dependence.

    Args:
        sn: NetworkStruct describing the queueing network
        options: Solver options

    Returns:
        SolverNCLDReturn with performance metrics

    Raises:
        RuntimeError: For unsupported network configurations
    """
    if options is None:
        options = SolverOptions()

    start_time = time.time()

    M = sn.nstations
    K = sn.nclasses
    C = sn.nchains

    nservers = sn.nservers
    if nservers is None:
        nservers = np.ones((M, 1))
    nservers = nservers.flatten()

    nservers_finite = np.where(np.isfinite(nservers), nservers, 0)
    min_finite_server = np.inf
    for val in nservers_finite:
        if val > 0 and val < min_finite_server:
            min_finite_server = val

    NK = sn.njobs.T if sn.njobs is not None else np.zeros((K, 1))
    if NK.ndim == 1:
        NK = NK.reshape(-1, 1)

    # Check for open classes (not supported)
    if np.any(np.isinf(NK)):
        raise RuntimeError("The load-dependent solver does not support open classes yet.")

    # Handle multi-server stations via lldscaling
    lldscaling = sn.lldscaling if hasattr(sn, 'lldscaling') and sn.lldscaling is not None else np.array([])

    # Only check multiserver if there are finite servers
    if np.isfinite(min_finite_server) and min_finite_server > 1:
        if (lldscaling.size == 0 and M == 2 and np.isfinite(np.max(np.abs(NK)))):
            # Auto-generate lldscaling for simple 2-station networks
            Nt = int(np.sum(NK))
            lldscaling = np.zeros((M, Nt))
            for i in range(M):
                for j in range(Nt):
                    lldscaling[i, j] = min(j + 1, nservers[i])
        else:
            raise RuntimeError("The load-dependent solver does not support multi-server stations yet. "
                             "Specify multi-server stations via limited load-dependence.")

    SCV = sn.scv

    # Compute visits matrix
    V = np.zeros((M, K))
    if sn.visits is not None:
        for key, visit_matrix in sn.visits.items():
            if visit_matrix is not None:
                if visit_matrix.shape == V.shape:
                    V = V + visit_matrix
                elif visit_matrix.shape[0] == V.shape[0]:
                    min_cols = min(visit_matrix.shape[1], V.shape[1])
                    V[:, :min_cols] += visit_matrix[:, :min_cols]

    # Compute service times from rates
    rates = sn.rates
    if rates is None:
        rates = np.ones((M, K))
    with np.errstate(divide='ignore', invalid='ignore'):
        ST = np.where(rates > 0, 1.0 / rates, 0.0)
        ST = np.where(np.isnan(ST), 0.0, ST)
    ST0 = ST.copy()

    # Initialize lldscaling if empty
    NK_finite = NK.copy()
    NK_finite = np.where(np.isfinite(NK_finite), NK_finite, 0)
    Nt = int(np.sum(NK_finite))

    if lldscaling.size == 0:
        lldscaling = np.ones((M, max(Nt, 1)))

    # Get chain demands
    from ...sn import sn_get_demands_chain
    demandsChainReturn = sn_get_demands_chain(sn)
    Vchain = demandsChainReturn.Vchain
    alpha = demandsChainReturn.alpha

    # Iteration control
    eta_1 = np.zeros(M)
    eta = np.ones(M)
    gamma = np.zeros(M)  # Non-exponential correction factor
    sched = sn.sched

    has_fcfs = False
    if sched is not None:
        for station_sched in sched.values():
            if station_sched == SchedStrategy.FCFS:
                has_fcfs = True
                break

    if not has_fcfs:
        options.iter_max = 1

    # Initialize output variables
    it = 0
    Lchain = None
    STchain = None
    Q = None
    U = None
    R = None
    T = None
    X = None
    lG = np.nan
    method = options.method

    # Main iteration loop
    while np.max(np.abs(1 - eta / (eta_1 + FINE_TOL))) > options.iter_tol and it < options.iter_max:
        it += 1
        eta_1 = eta.copy()

        # Compute chain-level demands
        Lchain = np.zeros((M, C))
        STchain = np.zeros((M, C))
        SCVchain = np.zeros((M, C))
        Nchain = np.zeros(C)
        refstatchain = np.zeros(C)

        for c in range(C):
            if c not in sn.inchain:
                continue
            inchain = sn.inchain[c].flatten().astype(int)

            is_open_chain = any(np.isinf(NK[k, 0]) for k in inchain if k < NK.shape[0])

            for i in range(M):
                ST_tmp = np.array([ST[i, k] for k in inchain if k < K])
                alpha_tmp = np.array([alpha[i, k] for k in inchain if k < K])
                SCV_tmp = np.array([SCV[i, k] for k in inchain if k < K]) if SCV is not None else np.ones(len(ST_tmp))

                if len(ST_tmp) > 0:
                    STchain[i, c] = np.dot(ST_tmp, alpha_tmp)
                    Lchain[i, c] = Vchain[i, c] * STchain[i, c]
                    SCVchain[i, c] = np.dot(SCV_tmp, alpha_tmp)

                    if is_open_chain:
                        refstat = sn.refstat
                        if refstat is not None and len(inchain) > 0:
                            first_class = inchain[0]
                            if first_class < len(refstat) and i == int(refstat[first_class]):
                                ST_finite = np.where(np.isfinite(ST_tmp), ST_tmp, 0)
                                STchain[i, c] = np.sum(ST_finite)

            NK_inchain = np.array([NK[k, 0] for k in inchain if k < NK.shape[0]])
            Nchain[c] = np.sum(NK_inchain)

            if len(inchain) > 0:
                first_class = inchain[0]
                refstat = sn.refstat
                if refstat is not None and first_class < len(refstat):
                    refstatchain[c] = refstat[first_class]

        # Clean up infinities
        STchain = np.where(np.isfinite(STchain), STchain, 0.0)
        Lchain = np.where(np.isfinite(Lchain), Lchain, 0.0)

        # Compute normalizing constant
        Nchain_finite = np.where(np.isfinite(Nchain), Nchain, 0)
        Nt = int(np.sum(Nchain_finite))

        L = np.zeros((M, C))
        mu = np.ones((M, max(Nt, 1)))
        Z = np.zeros((M, C))
        inf_servers = []

        for i in range(M):
            if np.isinf(nservers[i]):
                inf_servers.append(i)
                L[i, :] = Lchain[i, :]
                Z[i, :] = Lchain[i, :]
                for j in range(max(Nt, 1)):
                    mu[i, j] = j + 1
            else:
                L[i, :] = Lchain[i, :]
                for j in range(max(Nt, 1)):
                    if j < lldscaling.shape[1]:
                        mu[i, j] = lldscaling[i, j]
                    else:
                        mu[i, j] = lldscaling[i, -1] if lldscaling.shape[1] > 0 else 1.0

        # Compute NC
        from ...pfqn import pfqn_ncld, pfqn_mushift, pfqn_fnc
        Nchain0 = np.zeros(C)
        pfqn_options = {'method': options.method, 'tol': options.tol}


        nc_result = pfqn_ncld(L, Nchain, Nchain0, mu, pfqn_options)
        lG = nc_result.lG
        method = nc_result.method

        # Compute throughputs
        Qchain = np.zeros((M, C))
        Xchain = np.zeros(C)

        for r in range(C):
            if Nchain[r] <= 0:
                continue

            Nchain_r = Nchain.copy()
            Nchain_r[r] -= 1

            lGr = pfqn_ncld(L, Nchain_r, Nchain0, mu, pfqn_options).lG
            if np.isfinite(lG) and np.isfinite(lGr):
                Xchain[r] = exp(lGr - lG)
            else:
                Xchain[r] = 0.0

            # Compute queue lengths
            if M == 2 and any(np.isinf(nservers)):
                # Simple 2-station case with delay
                first_delay = -1
                for i in range(M):
                    if np.isinf(nservers[i]):
                        first_delay = i
                        break

                if first_delay >= 0:
                    Qchain[first_delay, r] = Lchain[first_delay, r] * Xchain[r]
                    for i in range(M):
                        if i != first_delay:
                            Qchain[i, r] = Nchain[r] - Lchain[first_delay, r] * Xchain[r]
            else:
                # General case
                for i in range(M):
                    if Lchain[i, r] > 0:
                        if np.isinf(nservers[i]):
                            Qchain[i, r] = Lchain[i, r] * Xchain[r]
                        else:
                            # Check if this is the last finite server
                            nservers_finite_arr = [nservers[j] for j in range(M) if np.isfinite(nservers[j])]
                            if i == M - 1 and len(nservers_finite_arr) == 1:
                                # Population balance
                                Lchainsum = sum(Lchain[j, r] for j in range(M) if np.isinf(nservers[j]))
                                Qchainsum = sum(Qchain[j, r] for j in range(M) if not np.isinf(nservers[j]) and j != i)
                                Qchain[i, r] = max(0, Nchain[r] - Lchainsum * Xchain[r] - Qchainsum)
                            else:
                                # Use functional server approximation
                                muhati = pfqn_mushift(mu, i)
                                # Use muhati[i:i+1, :] for fnc to maintain column dimension consistency
                                fnc_result = pfqn_fnc(muhati[i:i+1, :] if i < muhati.shape[0] else muhati[:1, :])
                                muhati_f = fnc_result.mu
                                c_val = fnc_result.c[0, 0] if fnc_result.c.size > 0 else 0.0

                                # Compute L and mu with station i removed (for marginal analysis)
                                # Lms_i = L with row i removed
                                # mu_i = mu with row i removed
                                Lms_i = np.delete(L, i, axis=0)
                                mu_i = np.delete(mu, i, axis=0)


                                # Compute queue length using marginal analysis
                                lGhatir = pfqn_ncld(L, Nchain_r, Nchain0, muhati, pfqn_options).lG
                                lGr_i = pfqn_ncld(Lms_i, Nchain_r, Nchain0, mu_i, pfqn_options).lG

                                # Extended L and mu
                                L_ext = np.vstack([L, L[i:i+1, :]])
                                muhati_ext = np.vstack([muhati, muhati_f])


                                lGhat_fnci = pfqn_ncld(L_ext, Nchain_r, Nchain0, muhati_ext, pfqn_options).lG

                                dlGa = lGhat_fnci - lGhatir
                                dlG_i = lGr_i - lGhatir

                                CQchain_i = (exp(dlGa) - 1) + c_val * (exp(dlG_i) - 1) if np.isfinite(dlGa) and np.isfinite(dlG_i) else 0


                                if mu[i, 0] > 0:
                                    ldDemand = log(L[i, r]) + lGhatir - log(mu[i, 0]) - lGr if L[i, r] > 0 and np.isfinite(lGhatir) and np.isfinite(lGr) else -np.inf
                                    if np.isfinite(ldDemand):
                                        Qchain[i, r] = exp(ldDemand) * Xchain[r] * (1 + CQchain_i)
                                    else:
                                        Qchain[i, r] = 0


        # Compute response times
        Z_sum = np.sum(Z, axis=0)
        Rchain = np.zeros((M, C))
        for i in range(M):
            for c in range(C):
                if Xchain[c] > ZERO and Vchain[i, c] > ZERO:
                    Rchain[i, c] = Qchain[i, c] / (Xchain[c] * Vchain[i, c])
                else:
                    Rchain[i, c] = 0.0

        # Fix response times for delay stations
        for i in inf_servers:
            for c in range(C):
                if Vchain[i, c] > ZERO:
                    Rchain[i, c] = Lchain[i, c] / Vchain[i, c]

        # Compute throughputs per station
        Tchain = Vchain * Xchain


        # Disaggregate to class level
        # Note: MATLAB passes [] for Qchain, so Q is computed from Rchain
        # This matches MATLAB solver_ncld.m line 205
        from ...sn import sn_deaggregate_chain_results
        deagg_result = sn_deaggregate_chain_results(
            sn, Lchain, ST, STchain, Vchain, alpha,
            None, None, Rchain, Tchain, None, Xchain.reshape(1, -1)
        )
        Q = deagg_result.Q
        U = deagg_result.U
        R = deagg_result.R
        T = deagg_result.T
        X = deagg_result.X


        # Non-exponential approximation
        from ...npfqn import npfqn_nonexp_approx
        nonexp_result = npfqn_nonexp_approx(
            method=options.highvar,
            sn=sn,
            ST=ST0,  # Use original service times
            V=V,
            SCV=SCV if SCV is not None else np.ones((M, K)),
            Tin=T,
            Uin=U,
            gamma=gamma,
            nservers=nservers
        )
        ST = nonexp_result.ST
        gamma = nonexp_result.gamma
        eta = nonexp_result.eta

    # Final cleanup
    if Q is not None:
        Q = np.abs(Q)
        Q = np.where(np.isfinite(Q), Q, 0.0)
    if R is not None:
        R = np.abs(R)
        R = np.where(np.isfinite(R), R, 0.0)
    if U is not None:
        U = np.abs(U)
        U = np.where(np.isfinite(U), U, 0.0)
    if X is not None:
        X = np.abs(X)
        X = np.where(np.isfinite(X), X, 0.0)
    if T is not None:
        T = np.where(np.isfinite(T), T, 0.0)

    # Utilization correction for multi-server and infinite server stations
    for i in range(M):
        if U is None:
            continue
        if nservers[i] > 1 and np.isfinite(nservers[i]):
            # Multi-server utilization correction
            for r in range(K):
                if X is not None and X[0, r] > 0:
                    c_idx = None
                    if sn.chains is not None:
                        for c in range(C):
                            if c in sn.inchain:
                                inchain = sn.inchain[c].flatten().astype(int)
                                if r in inchain:
                                    c_idx = c
                                    break
                    if c_idx is not None and c_idx in sn.visits:
                        visit_c = sn.visits[c_idx]
                        refstat = sn.refstat
                        if refstat is not None and r < len(refstat):
                            ref_idx = int(refstat[r])
                            if visit_c[ref_idx, r] > 0:
                                U[i, r] = X[0, r] * visit_c[i, r] / visit_c[ref_idx, r] * ST[i, r] / nservers[i]
        elif np.isinf(nservers[i]):
            # Infinite server utilization correction
            # Match MATLAB solver_ncld.m lines 230-244
            for r in range(K):
                if NK[r, 0] == np.inf:
                    # Open class - use lambda
                    from ...sn import sn_get_product_form_params
                    lambda_arr, _ = sn_get_product_form_params(sn)
                    if lambda_arr[r] > 0:
                        c_idx = None
                        if sn.chains is not None:
                            for c in range(C):
                                if c in sn.inchain:
                                    inchain = sn.inchain[c].flatten().astype(int)
                                    if r in inchain:
                                        c_idx = c
                                        break
                        if c_idx is not None and c_idx in sn.visits:
                            visit_c = sn.visits[c_idx]
                            refstat = sn.refstat
                            if refstat is not None and r < len(refstat):
                                ref_idx = int(refstat[r])
                                if visit_c[ref_idx, r] > 0:
                                    U[i, r] = lambda_arr[r] * visit_c[i, r] / visit_c[ref_idx, r] * ST[i, r]
                else:
                    # Closed class - use X
                    if X is not None and X[0, r] > 0:
                        c_idx = None
                        if sn.chains is not None:
                            for c in range(C):
                                if c in sn.inchain:
                                    inchain = sn.inchain[c].flatten().astype(int)
                                    if r in inchain:
                                        c_idx = c
                                        break
                        if c_idx is not None and c_idx in sn.visits:
                            visit_c = sn.visits[c_idx]
                            refstat = sn.refstat
                            if refstat is not None and r < len(refstat):
                                ref_idx = int(refstat[r])
                                if visit_c[ref_idx, r] > 0:
                                    U[i, r] = X[0, r] * visit_c[i, r] / visit_c[ref_idx, r] * ST[i, r]
        else:
            # Single server with lldscaling
            if lldscaling.shape[0] > i:
                max_lld = np.max(lldscaling[i, :]) if lldscaling.shape[1] > 0 else 1.0
                if max_lld > 0:
                    for j in range(K):
                        U[i, j] = U[i, j] / max_lld
                # Ensure utilization sum <= 1
                U_sum = np.sum(U[i, :])
                if U_sum > 1:
                    U_notnan = U[i, :].copy()
                    U_notnan = np.where(np.isnan(U_notnan), 0, U_notnan)
                    if np.sum(U_notnan) > 0:
                        U[i, :] = U[i, :] / np.sum(U_notnan)

    # Clean up X and U
    if X is not None:
        X = np.where(np.isinf(X) | np.isnan(X), 0.0, X)
    if U is not None:
        U = np.where(np.isinf(U) | np.isnan(U), 0.0, U)

    # Chain population corrections
    Nchain_final = np.zeros(C)
    for c in range(C):
        if c not in sn.inchain:
            continue
        inchain = sn.inchain[c].flatten().astype(int)
        Nchain_final[c] = sum(NK[k, 0] for k in inchain if k < NK.shape[0])

    for c in range(C):
        if c not in sn.inchain:
            continue
        inchain = sn.inchain[c].flatten().astype(int)
        Nchain_c = Nchain_final[c]

        if np.isfinite(Nchain_c) and Q is not None:
            sum_Q = sum(np.sum(Q[:, k]) for k in inchain if k < K)

            if sum_Q > 0:
                ratio = Nchain_c / sum_Q
                for k in inchain:
                    if k >= K:
                        continue
                    Q[:, k] *= ratio
                    if X is not None:
                        X[0, k] *= ratio
                    if T is not None:
                        T[:, k] *= ratio
                    if U is not None:
                        U[:, k] *= ratio
                    # Recalculate R = Q / T (MATLAB line 270)
                    if R is not None and T is not None:
                        for i in range(M):
                            if T[i, k] > 0:
                                R[i, k] = Q[i, k] / T[i, k]

    # Handle self-looping classes
    if hasattr(sn, 'isslc') and sn.isslc is not None:
        isslc = sn.isslc.flatten() if sn.isslc.ndim > 1 else sn.isslc
        for k in range(K):
            if k < len(isslc) and isslc[k] == 1:
                if Q is not None:
                    Q[:, k] = 0.0
                    refstat = sn.refstat
                    if refstat is not None and k < len(refstat):
                        ist = int(refstat[k])
                        if ist < M:
                            Q[ist, k] = NK[k, 0]
                            if T is not None and rates is not None:
                                T[ist, k] = NK[k, 0] * rates[ist, k]
                            if R is not None and T is not None and T[ist, k] > 0:
                                R[ist, k] = Q[ist, k] / T[ist, k]
                            if U is not None:
                                U[ist, k] = ST[ist, k] * T[ist, k] if T is not None else 0

    runtime = time.time() - start_time

    return SolverNCLDReturn(
        Q=Q, U=U, R=R, T=T, C=None, X=X, lG=lG,
        runtime=runtime, it=it, method=method
    )


@dataclass
class SolverNCCacheReturn:
    """Result of NC cache analyzer."""
    QN: Optional[np.ndarray]  # Mean queue lengths (M x K)
    UN: Optional[np.ndarray]  # Utilizations (M x K)
    RN: Optional[np.ndarray]  # Response times (M x K)
    TN: Optional[np.ndarray]  # Throughputs (M x K)
    CN: Optional[np.ndarray]  # Covariances (optional)
    XN: Optional[np.ndarray]  # System throughputs (1 x K)
    lG: float                # Log normalizing constant
    pij: Optional[np.ndarray]  # Cache hit probabilities
    runtime: float           # Runtime in seconds
    method: str              # Method used


def solver_nc_cache_analyzer(sn: NetworkStruct, options: Optional[SolverOptions] = None) -> SolverNCCacheReturn:
    """
    NC solver for cache-only networks.

    Analyzes a network containing only a Source, Cache, and Sink using
    product-form analysis for cache systems.

    Args:
        sn: NetworkStruct describing the network with a Cache node
        options: Solver options (method: 'exact', 'sampling', 'spm'/default)

    Returns:
        SolverNCCacheReturn with cache performance metrics

    References:
        Original MATLAB: matlab/src/solvers/NC/solver_nc_cache_analyzer.m
    """
    from ...cache import (
        cache_gamma_lp,
        cache_prob_erec,
        cache_prob_spm,
        cache_miss_spm,
        cache_miss_is,
        cache_prob_is,
    )

    if options is None:
        options = SolverOptions()

    start_time = time.time()

    M = sn.nstations
    K = sn.nclasses

    # Initialize outputs
    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))
    CN = None
    XN = np.zeros(K)
    lG = np.nan
    pij = None
    method = options.method

    # Find Source station and get source rates
    source_ist = -1
    if sn.nodetype is not None:
        from ...sn import NodeType
        for node_idx, nt in enumerate(sn.nodetype):
            if nt == NodeType.SOURCE or nt == int(NodeType.SOURCE):
                source_ist = int(sn.nodeToStation[node_idx]) if hasattr(sn, 'nodeToStation') else node_idx
                break

    if source_ist < 0:
        raise RuntimeError("NC cache analyzer requires a Source station")

    # Get source rates
    rates = sn.rates if sn.rates is not None else np.zeros((M, K))
    sourceRate = rates[source_ist, :].flatten()
    sourceRate = np.where(np.isnan(sourceRate), 0.0, sourceRate)

    # Set throughput at source
    TN[source_ist, :] = sourceRate

    # Find Cache node and get parameters
    cache_idx = -1
    if sn.nodetype is not None:
        from ...sn import NodeType
        for node_idx, nt in enumerate(sn.nodetype):
            if nt == NodeType.CACHE or nt == int(NodeType.CACHE):
                cache_idx = node_idx
                break

    if cache_idx < 0:
        raise RuntimeError("NC cache analyzer requires a Cache node")

    # Get cache parameters from nodeparam
    if sn.nodeparam is None or cache_idx not in sn.nodeparam:
        raise RuntimeError("Cache node parameters not found in nodeparam")

    ch = sn.nodeparam[cache_idx]

    # Check replacement strategy - NC only supports RR (0) and FIFO (1)
    replstrat = getattr(ch, 'replacestrat', 0)
    # Handle enum types
    if hasattr(replstrat, 'value'):
        replstrat = replstrat.value
    replstrat = int(replstrat)

    # RR = 0, FIFO = 1; LRU = 3 and others are not supported
    if replstrat not in (0, 1):
        raise RuntimeError("[solver_nc_cache_analyzer] NC does not support approximate solution of the specified cache replacement policy.")

    # Get capacity and items
    m = np.atleast_1d(getattr(ch, 'itemcap', np.array([]))).flatten()
    n = int(getattr(ch, 'nitems', 0))

    if n == 0 or len(m) == 0:
        raise RuntimeError("Cache has invalid capacity or items count")

    total_cap = int(np.sum(m))
    if n < total_cap + 2:
        raise RuntimeError(f"NC requires the number of items ({n}) to exceed the cache capacity ({total_cap}) at least by 2")

    h = len(m)  # Number of cache levels
    u = K  # Number of classes (users)

    # Get pread (item popularity per class)
    pread_data = getattr(ch, 'pread', None)

    # Build lambda array (u x n x h+1)
    # lambda[v,k,l] = sourceRate[v] * pread[v][k]
    lambd = np.zeros((u, n, h + 1))

    for v in range(u):
        # Get pread for this class
        pread_v = None
        if pread_data is not None:
            if isinstance(pread_data, dict):
                pread_v = pread_data.get(v)
            elif isinstance(pread_data, (list, tuple)) and v < len(pread_data):
                pread_v = pread_data[v]
            elif isinstance(pread_data, np.ndarray):
                if pread_data.ndim == 1:
                    pread_v = pread_data
                elif v < pread_data.shape[0]:
                    pread_v = pread_data[v, :]

        if pread_v is not None:
            pread_v = np.atleast_1d(pread_v).flatten()
            for k in range(min(n, len(pread_v))):
                for l in range(h + 1):
                    lambd[v, k, l] = sourceRate[v] * pread_v[k]

    # Get routing structure (accost)
    R_data = getattr(ch, 'accost', None)

    # Build R structure for cache_gamma_lp
    # R is a list of lists: R[v][k] is the routing matrix for user v, item k
    # For simple cache models, we use a default diagonal structure
    R = [[None for _ in range(n)] for _ in range(u)]

    if R_data is not None:
        # R_data can be:
        # - A cell array (dict/list) with R[v][k] being a (h+1 x h+1) matrix
        # - A single matrix to be used for all users/items
        if isinstance(R_data, dict):
            for v in range(u):
                for k in range(n):
                    key = (v, k)
                    if key in R_data:
                        R[v][k] = np.asarray(R_data[key])
        elif isinstance(R_data, (list, tuple)):
            # Assume it's organized as R[v][k]
            for v in range(min(u, len(R_data))):
                if R_data[v] is not None:
                    if isinstance(R_data[v], (list, tuple)):
                        for k in range(min(n, len(R_data[v]))):
                            if R_data[v][k] is not None:
                                R[v][k] = np.asarray(R_data[v][k])
                    else:
                        # Single matrix for all items of this user
                        for k in range(n):
                            R[v][k] = np.asarray(R_data[v])
        elif isinstance(R_data, np.ndarray):
            # Single matrix for all users/items
            for v in range(u):
                for k in range(n):
                    R[v][k] = R_data.copy()

    # Fill in default routing if not specified
    # Default: sequential access from level 0 to level h
    for v in range(u):
        for k in range(n):
            if R[v][k] is None:
                # Default routing: diag(ones(1,h),1) with self-loop at end
                R_default = np.zeros((h + 1, h + 1))
                for l in range(h):
                    R_default[l, l + 1] = 1.0
                R_default[h, h] = 1.0  # Self-loop at last level
                R[v][k] = R_default

    # Compute gamma using cache_gamma_lp
    gamma, _, _, _ = cache_gamma_lp(lambd, R)

    # Compute miss rates based on method
    if options.method == 'exact':
        # Exact method using recursive computation
        pij = cache_prob_erec(gamma, m)
        missRate = np.zeros(u)
        for v in range(u):
            # missRate[v] = sum_k lambda[v,k,0] * pij[k,0]
            missRate[v] = np.dot(lambd[v, :, 0], pij[:, 0])
        method = 'exact'

    elif options.method == 'sampling':
        # Sampling method
        samples = getattr(options, 'samples', 10000)
        if samples is None:
            samples = 10000
        try:
            _, missRate_arr, _, _, lE = cache_miss_is(gamma, m, lambd, samples)
            missRate = missRate_arr.flatten() if missRate_arr is not None else np.zeros(u)
            pij = cache_prob_is(gamma, m, samples)
        except Exception:
            # Fallback to SPM
            _, missRate_arr, _, _, lE = cache_miss_spm(gamma, m, lambd)
            missRate = missRate_arr.flatten() if missRate_arr is not None else np.zeros(u)
            pij = cache_prob_spm(gamma, m)
        method = 'sampling'

    else:
        # Default: SPM (singular perturbation method)
        _, missRate_arr, _, _, lE = cache_miss_spm(gamma, m, lambd)
        missRate = missRate_arr.flatten() if missRate_arr is not None else np.zeros(u)
        pij = cache_prob_spm(gamma, m)
        method = 'spm'

    # Get hitclass and missclass
    hitclass = np.atleast_1d(getattr(ch, 'hitclass', np.array([]))).flatten().astype(int)
    missclass = np.atleast_1d(getattr(ch, 'missclass', np.array([]))).flatten().astype(int)

    # Set throughputs for hit and miss classes
    # Also find the Cache station index
    cache_ist = -1
    if hasattr(sn, 'nodeToStation') and cache_idx >= 0:
        cache_ist = int(sn.nodeToStation[cache_idx])

    # Compute and store actual hit/miss probabilities
    actualhitprob = np.zeros(K)
    actualmissprob = np.zeros(K)

    for r in range(K):
        if r < len(hitclass) and r < len(missclass):
            mc = int(missclass[r])
            hc = int(hitclass[r])
            # Python uses 0-based indexing (hitclass/missclass are already 0-based)
            # Valid indices are >= 0 and < K
            if mc >= 0 and hc >= 0 and mc < K and hc < K:
                miss_tput = missRate[r]
                hit_tput = sourceRate[r] - missRate[r]

                XN[mc] = XN[mc] + miss_tput
                XN[hc] = XN[hc] + hit_tput

                # Also set TN at the Cache station for hit/miss classes
                if cache_ist >= 0:
                    TN[cache_ist, mc] = miss_tput
                    TN[cache_ist, hc] = hit_tput

                # Compute hit/miss probabilities for this requesting class
                if sourceRate[r] > 0:
                    actualhitprob[r] = hit_tput / sourceRate[r]
                    actualmissprob[r] = miss_tput / sourceRate[r]

    # Store hit/miss probabilities in nodeparam for use by sn_get_node_arvr_from_tput
    ch.actualhitprob = actualhitprob
    ch.actualmissprob = actualmissprob

    # Compute throughputs for downstream stations by following routing from Cache
    # This handles routing through non-station nodes like Router
    if hasattr(sn, 'rt') and sn.rt is not None and cache_ist >= 0:
        rt_full = np.asarray(sn.rt)
        nstateful = sn.nstateful

        # Build station_to_stateful mapping
        station_to_stateful = np.zeros(M, dtype=int)
        if hasattr(sn, 'stationToNode') and sn.stationToNode is not None and \
           hasattr(sn, 'nodeToStateful') and sn.nodeToStateful is not None:
            for ist in range(M):
                node_idx = int(sn.stationToNode[ist]) if ist < len(sn.stationToNode) else ist
                stateful_idx = int(sn.nodeToStateful[node_idx]) if node_idx < len(sn.nodeToStateful) else ist
                station_to_stateful[ist] = stateful_idx

        # Build stateful_to_station mapping: -1 for non-station stateful nodes
        stateful_to_station = np.full(nstateful, -1, dtype=int)
        for ist in range(M):
            sf = station_to_stateful[ist]
            if sf < nstateful:
                stateful_to_station[sf] = ist

        def resolve_routing_to_stations(src_sf: int, src_class: int):
            """Resolve routing through non-station nodes to find downstream stations."""
            result = {}
            queue = []
            src_idx = src_sf * K + src_class
            for dst_sf in range(nstateful):
                for dst_class in range(K):
                    dst_idx = dst_sf * K + dst_class
                    if src_idx < rt_full.shape[0] and dst_idx < rt_full.shape[1]:
                        prob = rt_full[src_idx, dst_idx]
                        if prob > 1e-10:
                            queue.append((dst_sf, dst_class, prob))

            max_iters = 100
            iters = 0
            while queue and iters < max_iters:
                iters += 1
                sf, cls, prob = queue.pop(0)
                station = stateful_to_station[sf] if sf < len(stateful_to_station) else -1
                if station >= 0:
                    key = (station, cls)
                    result[key] = result.get(key, 0.0) + prob
                else:
                    idx = sf * K + cls
                    for next_sf in range(nstateful):
                        for next_cls in range(K):
                            next_idx = next_sf * K + next_cls
                            if idx < rt_full.shape[0] and next_idx < rt_full.shape[1]:
                                next_prob = rt_full[idx, next_idx]
                                if next_prob > 1e-10:
                                    queue.append((next_sf, next_cls, prob * next_prob))
            return result

        # Compute throughputs for hit/miss classes at downstream stations
        cache_sf = station_to_stateful[cache_ist]
        for r in range(K):
            if r < len(hitclass) and r < len(missclass):
                hc = int(hitclass[r])
                mc = int(missclass[r])
                if hc >= 0 and mc >= 0 and hc < K and mc < K:
                    # Route hit class from cache to downstream
                    hit_routing = resolve_routing_to_stations(cache_sf, hc)
                    for (dst_ist, dst_cls), prob in hit_routing.items():
                        if dst_ist != cache_ist:  # Don't overwrite cache throughput
                            TN[dst_ist, dst_cls] += TN[cache_ist, hc] * prob

                    # Route miss class from cache to downstream
                    miss_routing = resolve_routing_to_stations(cache_sf, mc)
                    for (dst_ist, dst_cls), prob in miss_routing.items():
                        if dst_ist != cache_ist:  # Don't overwrite cache throughput
                            TN[dst_ist, dst_cls] += TN[cache_ist, mc] * prob

    # MATLAB: QN/UN/RN stay empty for standalone cache (Source→Cache→Sink)
    # No real queueing stations exist, so leave QN/UN/RN as zeros.

    runtime = time.time() - start_time

    return SolverNCCacheReturn(
        QN=QN,
        UN=UN,
        RN=RN,
        TN=TN,
        CN=CN,
        XN=XN.reshape(1, -1),
        lG=lG,
        pij=pij,
        runtime=runtime,
        method=method
    )


@dataclass
class SolverNCCacheQNReturn:
    """Result of NC cache+queueing analyzer."""
    QN: Optional[np.ndarray]  # Mean queue lengths (M x K)
    UN: Optional[np.ndarray]  # Utilizations (M x K)
    RN: Optional[np.ndarray]  # Response times (M x K)
    TN: Optional[np.ndarray]  # Throughputs (M x K)
    CN: Optional[np.ndarray]  # Covariances (optional)
    XN: Optional[np.ndarray]  # System throughputs (1 x K)
    lG: float                # Log normalizing constant
    hitprob: Optional[np.ndarray]  # Hit probabilities per cache per class
    missprob: Optional[np.ndarray]  # Miss probabilities per cache per class
    runtime: float           # Runtime in seconds
    it: int                  # Number of iterations
    method: str              # Method used
    visits: Optional[dict] = None       # Updated visit ratios (stateful)
    nodevisits: Optional[dict] = None   # Updated visit ratios (all nodes)


def solver_nc_cacheqn_analyzer(
    sn,
    options: Optional[SolverOptions] = None
) -> SolverNCCacheQNReturn:
    """
    NC solver for queueing networks with integrated cache nodes.

    This implements the iterative cache-QN analyzer that solves networks
    containing both Cache nodes and queueing stations. Unlike the standalone
    cache analyzer (solver_nc_cache_analyzer), this handles mixed networks
    where cache hit/miss affects routing to downstream queues.

    Algorithm:
    1. Find Cache nodes and extract parameters
    2. Convert Cache nodes to ClassSwitch (routing handled by rtnodes)
    3. Initialize random lambda (arrival rates at caches)
    4. Iterate until convergence:
       a. Solve isolated cache: build lambda_cache, compute gamma via cache_gamma_lp
       b. For 'exact' method: compute pij = cache_prob_erec, derive miss rates
       c. Compute hit/miss probabilities from miss rates and lambda
       d. Update rtnodes with hit/miss routing probabilities
       e. Recompute rt = dtmc_stochcomp(rtnodes, statefulNodesClasses)
       f. Refresh visits: sn_refresh_visits(sn, chains, rt, rtnodes)
       g. Solve QN: solver_nc_analyzer or solver_ncld_analyzer → get XN
       h. Update lambda from XN and nodevisits
       i. Check convergence: norm(lambda - lambda_prev) < iter_tol

    Args:
        sn: NetworkStruct describing the queueing network with Cache nodes
        options: Solver options (method: 'exact' for cache_prob_erec, else SPM)

    Returns:
        SolverNCCacheQNReturn with performance metrics and hit/miss probabilities

    References:
        MATLAB: matlab/src/solvers/NC/solver_nc_cacheqn_analyzer.m
    """
    from ...cache import (
        cache_gamma_lp,
        cache_prob_erec,
        cache_miss_spm,
    )
    from ...mc.dtmc import dtmc_stochcomp
    from ...sn.transforms import sn_refresh_visits
    from ...sn.network_struct import NodeType

    if options is None:
        options = SolverOptions()

    start_time = time.time()

    # Make a working copy of sn to avoid modifying original
    # (MATLAB: snorig = self.model.getStruct; sn = snorig;)
    import copy
    sn_work = copy.deepcopy(sn)

    I = sn_work.nnodes
    K = sn_work.nclasses
    M = sn_work.nstations

    # Find stateful nodes and build index mapping
    stateful_nodes = []
    if sn_work.isstateful is not None:
        for i in range(I):
            if i < len(sn_work.isstateful) and sn_work.isstateful[i]:
                stateful_nodes.append(i)
    else:
        stateful_nodes = list(range(I))

    stateful_nodes_classes = []
    for ind in stateful_nodes:
        for k in range(K):
            stateful_nodes_classes.append(ind * K + k)
    stateful_nodes_classes = np.array(stateful_nodes_classes, dtype=int)

    # Find cache nodes
    caches = []
    if sn_work.nodetype is not None:
        for ind in range(I):
            if ind < len(sn_work.nodetype):
                if sn_work.nodetype[ind] == NodeType.CACHE:
                    caches.append(ind)

    if not caches:
        # No cache nodes - fall back to standard NC
        nc_result = solver_nc(sn, options)
        return SolverNCCacheQNReturn(
            QN=nc_result.Q,
            UN=nc_result.U,
            RN=nc_result.R,
            TN=nc_result.T,
            CN=None,
            XN=nc_result.X,
            lG=nc_result.lG,
            hitprob=None,
            missprob=None,
            runtime=nc_result.runtime,
            it=nc_result.it,
            method=nc_result.method
        )

    # Initialize arrays
    lambd = np.zeros(K)
    lambd_prev = np.zeros(K)
    hitprob_all = np.zeros((len(caches), K))
    missprob_all = np.zeros((len(caches), K))
    missrate = np.zeros((len(caches), K))

    # Main iteration loop
    QN, UN, RN, TN, XN = None, None, None, None, None
    lG = np.nan
    method = options.method

    iter_max = options.iter_max if hasattr(options, 'iter_max') else 10
    iter_tol = options.iter_tol if hasattr(options, 'iter_tol') else 1e-4

    for it in range(1, iter_max + 1):
        for cache_idx, ind in enumerate(caches):
            # Get cache parameters from nodeparam
            ch = sn_work.nodeparam.get(ind) if sn_work.nodeparam else None
            if ch is None:
                continue

            hitclass = getattr(ch, 'hitclass', None)
            missclass = getattr(ch, 'missclass', None)
            if hitclass is None or missclass is None:
                continue

            # Convert to numpy arrays
            hitclass = np.atleast_1d(hitclass).flatten().astype(int)
            missclass = np.atleast_1d(missclass).flatten().astype(int)

            # Find input classes (classes with valid hit/miss mappings)
            input_classes = []
            for r in range(K):
                if r < len(hitclass) and hitclass[r] >= 0:
                    input_classes.append(r)

            if not input_classes:
                continue

            m = np.atleast_1d(getattr(ch, 'itemcap', np.array([1]))).flatten().astype(int)
            n = int(getattr(ch, 'nitems', 1))

            if it == 1:
                # Check capacity constraint
                total_cap = int(np.sum(m))
                if n < total_cap + 2:
                    raise RuntimeError(
                        f"NC requires the number of items ({n}) to exceed "
                        f"the cache capacity ({total_cap}) at least by 2."
                    )

                # Initialize random lambda for input classes
                rng = np.random.default_rng(42)
                lambd_prev[input_classes] = rng.uniform(0.1, 1.0, len(input_classes))
                lambd = lambd_prev.copy()

                # Convert Cache node to ClassSwitch for QN solving
                sn_work.nodetype[ind] = NodeType.CLASSSWITCH

            # Build lambda_cache array: lambda[v,k,l] = lambda[v] * pread[v][k]
            h = len(m)  # Number of cache levels
            u = K  # Number of classes (users)
            lambda_cache = np.zeros((u, n, h + 1))

            pread_data = getattr(ch, 'pread', None)
            for v in range(u):
                pread_v = None
                if pread_data is not None:
                    if isinstance(pread_data, dict):
                        pread_v = pread_data.get(v)
                    elif isinstance(pread_data, (list, tuple)) and v < len(pread_data):
                        pread_v = pread_data[v]
                    elif isinstance(pread_data, np.ndarray):
                        if pread_data.ndim == 1:
                            pread_v = pread_data
                        elif v < pread_data.shape[0]:
                            pread_v = pread_data[v, :]

                if pread_v is not None:
                    pread_v = np.atleast_1d(pread_v).flatten()
                    for k in range(min(n, len(pread_v))):
                        if not np.isnan(pread_v[k]):
                            for l in range(h + 1):
                                lambda_cache[v, k, l] = lambd[v] * pread_v[k]

            # Get routing cost structure
            R_cost = getattr(ch, 'accost', None)

            # Build R structure for cache_gamma_lp
            R = [[None for _ in range(n)] for _ in range(u)]
            if R_cost is not None:
                if isinstance(R_cost, (list, tuple)):
                    for v in range(min(u, len(R_cost))):
                        if R_cost[v] is not None:
                            if isinstance(R_cost[v], (list, tuple)):
                                for k in range(min(n, len(R_cost[v]))):
                                    if R_cost[v][k] is not None:
                                        R[v][k] = np.asarray(R_cost[v][k])
                            else:
                                for k in range(n):
                                    R[v][k] = np.asarray(R_cost[v])
                elif isinstance(R_cost, np.ndarray):
                    for v in range(u):
                        for k in range(n):
                            R[v][k] = R_cost.copy()

            # Fill default routing if not specified
            for v in range(u):
                for k in range(n):
                    if R[v][k] is None:
                        R_default = np.zeros((h + 1, h + 1))
                        for l in range(h):
                            R_default[l, l + 1] = 1.0
                        R_default[h, h] = 1.0
                        R[v][k] = R_default

            # Compute gamma using cache_gamma_lp
            try:
                gamma, _, _, _ = cache_gamma_lp(lambda_cache, R)
            except Exception:
                # Fallback: uniform gamma
                gamma = np.ones((n, h + 1)) / n

            # Compute miss rates based on method
            if options.method == 'exact':
                # Exact method using recursive computation
                try:
                    pij = cache_prob_erec(gamma, m)
                    for v in range(u):
                        # missRate[v] = sum_k lambda_cache[v,k,0] * pij[k,0]
                        missrate[cache_idx, v] = np.dot(lambda_cache[v, :, 0], pij[:, 0])
                    method = 'exact'
                except Exception:
                    # Fallback to SPM
                    try:
                        _, missrate_arr, _, _, _ = cache_miss_spm(gamma, m, lambda_cache)
                        missrate[cache_idx, :] = missrate_arr.flatten() if missrate_arr is not None else np.zeros(u)
                    except Exception:
                        missrate[cache_idx, :] = 0.5 * lambd
                    method = 'spm'
            else:
                # Default: SPM method
                try:
                    _, missrate_arr, _, _, _ = cache_miss_spm(gamma, m, lambda_cache)
                    missrate[cache_idx, :] = missrate_arr.flatten() if missrate_arr is not None else np.zeros(u)
                except Exception as e:
                    missrate[cache_idx, :] = 0.5 * lambd
                method = 'spm'

            # Compute hit/miss probabilities
            for v in range(K):
                if lambd[v] > 0:
                    missprob_all[cache_idx, v] = missrate[cache_idx, v] / lambd[v]
                    hitprob_all[cache_idx, v] = 1.0 - missprob_all[cache_idx, v]
                else:
                    missprob_all[cache_idx, v] = 0.0
                    hitprob_all[cache_idx, v] = 0.0

            # Clean up NaN values
            hitprob_all = np.where(np.isnan(hitprob_all), 0.0, hitprob_all)
            missprob_all = np.where(np.isnan(missprob_all), 0.0, missprob_all)

            # Update routing matrix (MATLAB lines 82-89)
            # For each input class, route from Cache to connected nodes with hit/miss probabilities
            if sn_work.rtnodes is not None:
                for r in input_classes:
                    # Zero out the row for input class at cache
                    sn_work.rtnodes[(ind) * K + r, :] = 0

                    # Find connected nodes and set hit/miss routing
                    for jnd in range(I):
                        if sn_work.connmatrix is not None:
                            if ind < sn_work.connmatrix.shape[0] and jnd < sn_work.connmatrix.shape[1]:
                                if sn_work.connmatrix[ind, jnd]:
                                    hc = hitclass[r] if r < len(hitclass) else -1
                                    mc = missclass[r] if r < len(missclass) else -1
                                    if hc >= 0 and hc < K:
                                        sn_work.rtnodes[(ind) * K + r, (jnd) * K + hc] = hitprob_all[cache_idx, r]
                                    if mc >= 0 and mc < K:
                                        sn_work.rtnodes[(ind) * K + r, (jnd) * K + mc] = missprob_all[cache_idx, r]

                # Recompute rt using stochastic complement (MATLAB line 91)
                try:
                    new_rt = dtmc_stochcomp(sn_work.rtnodes, stateful_nodes_classes)
                    sn_work.rt = new_rt
                    # CRITICAL: Also update rt_visits since sn_refresh_visits uses it
                    if hasattr(sn_work, 'rt_visits') and sn_work.rt_visits is not None:
                        sn_work.rt_visits = new_rt.copy()
                except Exception:
                    pass

        # Refresh visits (MATLAB line 93-95)
        # Python sn_refresh_visits modifies sn in place and doesn't take extra arguments
        try:
            sn_refresh_visits(sn_work)
        except Exception:
            pass

        # Solve the QN (MATLAB lines 97-101)
        # Check for load-dependent scaling
        has_lld = hasattr(sn_work, 'lldscaling') and sn_work.lldscaling is not None and sn_work.lldscaling.size > 0
        has_cd = hasattr(sn_work, 'cdscaling') and sn_work.cdscaling is not None

        try:
            if has_lld or has_cd:
                nc_result = solver_ncld(sn_work, options)
                QN = nc_result.Q
                UN = nc_result.U
                RN = nc_result.R
                TN = nc_result.T
                XN = nc_result.X
                lG = nc_result.lG
            else:
                nc_result = solver_nc(sn_work, options)
                QN = nc_result.Q
                UN = nc_result.U
                RN = nc_result.R
                TN = nc_result.T
                XN = nc_result.X
                lG = nc_result.lG
        except Exception as e:
            # If QN solving fails, return partial results
            runtime = time.time() - start_time
            return SolverNCCacheQNReturn(
                QN=np.zeros((M, K)),
                UN=np.zeros((M, K)),
                RN=np.zeros((M, K)),
                TN=np.zeros((M, K)),
                CN=None,
                XN=np.zeros((1, K)),
                lG=np.nan,
                hitprob=hitprob_all,
                missprob=missprob_all,
                runtime=runtime,
                it=it,
                method=method
            )

        # Aggregate nodevisits across chains (MATLAB line 103: nodevisits = cellsum(nodevisits))
        # Note: nodevisits can be indexed by nodes or stateful nodes depending on how visits were computed
        if sn_work.nodevisits is not None:
            # Determine the shape from first non-None entry
            nodevisits_shape = None
            for chain_id, nv in sn_work.nodevisits.items():
                if nv is not None:
                    nodevisits_shape = nv.shape
                    break
            if nodevisits_shape is None:
                nodevisits_shape = (I, K)  # Default to nnodes
            nodevisits_sum = np.zeros(nodevisits_shape)
            for chain_id, nv in sn_work.nodevisits.items():
                if nv is not None and nv.shape == nodevisits_shape:
                    nodevisits_sum += nv
        else:
            nodevisits_sum = np.ones((I, K))

        # Update lambda from throughputs (MATLAB lines 104-114)
        for cache_idx, ind in enumerate(caches):
            ch = sn_work.nodeparam.get(ind) if sn_work.nodeparam else None
            if ch is None:
                continue

            hitclass = getattr(ch, 'hitclass', None)
            if hitclass is None:
                continue
            hitclass = np.atleast_1d(hitclass).flatten().astype(int)

            input_classes = [r for r in range(K) if r < len(hitclass) and hitclass[r] >= 0]

            for r in input_classes:
                # Find chain for class r
                c = -1
                if sn_work.chains is not None:
                    for chain_id in range(sn_work.nchains):
                        if chain_id in sn_work.inchain:
                            inchain = sn_work.inchain[chain_id].flatten().astype(int)
                            if r in inchain:
                                c = chain_id
                                break

                if c < 0:
                    continue

                inchain = sn_work.inchain[c].flatten().astype(int) if c in sn_work.inchain else [r]

                # Sum of throughputs in chain (MATLAB: sum(XN(inchain)))
                if XN is not None:
                    X_flat = XN.flatten()
                    Xchain_sum = sum(X_flat[k] for k in inchain if k < len(X_flat))
                else:
                    Xchain_sum = 0.0

                # Get reference station and class for this chain
                refstat_r = int(sn_work.refstat[r]) if sn_work.refstat is not None and r < len(sn_work.refstat) else 0
                refclass = int(sn_work.refclass[c]) if sn_work.refclass is not None and c < len(sn_work.refclass) else r

                # nodevisits_sum is indexed by node (not stateful), matching MATLAB's nodevisits(ind,r)
                refstat_node = int(sn_work.stationToNode[refstat_r]) if hasattr(sn_work, 'stationToNode') and refstat_r < len(sn_work.stationToNode) else refstat_r

                # Update lambda (MATLAB line 109-112)
                nv_cache_r = nodevisits_sum[ind, r] if ind < nodevisits_sum.shape[0] and r < nodevisits_sum.shape[1] else 1.0

                if refclass >= 0 and refclass < K:
                    nv_ref = nodevisits_sum[refstat_node, refclass] if refstat_node < nodevisits_sum.shape[0] and refclass < nodevisits_sum.shape[1] else 1.0
                else:
                    nv_ref = nodevisits_sum[refstat_node, r] if refstat_node < nodevisits_sum.shape[0] and r < nodevisits_sum.shape[1] else 1.0

                if nv_ref > 0:
                    lambd[r] = Xchain_sum * nv_cache_r / nv_ref
                else:
                    lambd[r] = Xchain_sum * nv_cache_r

        # Check convergence (MATLAB lines 115-118)
        if np.linalg.norm(lambd - lambd_prev, 1) < iter_tol:
            break

        lambd_prev = lambd.copy()

    runtime = time.time() - start_time

    # Ensure output dimensions match stations (M x K)
    if QN is not None and QN.shape[0] != M:
        # Results are in stateful node dimensions, need to extract station results
        QN_stations = np.zeros((M, K))
        UN_stations = np.zeros((M, K))
        RN_stations = np.zeros((M, K))
        TN_stations = np.zeros((M, K))
        for ist in range(M):
            if hasattr(sn_work, 'stationToStateful'):
                sf_idx = int(sn_work.stationToStateful[ist])
                if sf_idx < QN.shape[0]:
                    QN_stations[ist, :] = QN[sf_idx, :]
                    UN_stations[ist, :] = UN[sf_idx, :] if UN is not None else 0
                    RN_stations[ist, :] = RN[sf_idx, :] if RN is not None else 0
                    TN_stations[ist, :] = TN[sf_idx, :] if TN is not None else 0
        QN = QN_stations
        UN = UN_stations
        RN = RN_stations
        TN = TN_stations

    return SolverNCCacheQNReturn(
        QN=QN,
        UN=UN,
        RN=RN,
        TN=TN,
        CN=None,
        XN=XN,
        lG=lG,
        hitprob=hitprob_all,
        missprob=missprob_all,
        runtime=runtime,
        it=it,
        method=method,
        visits=sn_work.visits,
        nodevisits=sn_work.nodevisits,
    )


__all__ = [
    'solver_nc',
    'solver_ncld',
    'solver_nc_cache_analyzer',
    'solver_nc_cacheqn_analyzer',
    'SolverNCReturn',
    'SolverNCLDReturn',
    'SolverNCCacheReturn',
    'SolverNCCacheQNReturn',
    'SolverOptions',
]
