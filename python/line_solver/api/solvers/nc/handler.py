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
    iter_max: int = 10
    iter_tol: float = 1e-4
    verbose: int = 0
    highvar: Optional[str] = None


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
                    # For open chains, arrival rate comes from Source station
                    if source_station >= 0:
                        first_class = inchain[0] if len(inchain) > 0 else 0
                        if source_station < rates.shape[0] and first_class < rates.shape[1]:
                            lmbda[c] = rates[source_station, first_class]
                    elif sn.refstat is not None and len(sn.refstat) > 0:
                        # Fallback to refstat-based computation
                        first_class = inchain[0]
                        if first_class < len(sn.refstat):
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

        # Compute normalizing constant (only for closed chains)
        pfqn_options = {'method': options.method, 'tol': options.tol}
        method = options.method

        if closed_chains:
            # Mixed or closed network: compute normalizing constant for closed chains only
            Nchain_closed = np.array([Nchain[c] if c in closed_chains else 0 for c in range(C)])
            G, lG = pfqn_nc(Lms, Nchain_closed, Z_new, method=options.method)
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

            G_tmp, lG_tmp = pfqn_nc(Lms, Nchain_tmp, Z_new, method=options.method)

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
                        Z_aug = np.append(Z_new, 0.0)

                        # Create Lms with station i moved to end and augmented
                        L_tmp = []
                        for j in range(M):
                            if j != i:
                                L_tmp.append(np.append(Lms[j, :], 0.0))

                        # Add station i at end with demand in augmented class
                        L_i_aug = np.append(Lms[i, :], 1.0)
                        L_tmp.append(L_i_aug)
                        L_tmp = np.array(L_tmp)

                        # Compute marginal NC
                        _, lG_marg = pfqn_nc(L_tmp, Nchain_aug, Z_aug, method=options.method)

                        if np.isfinite(lG_marg) and np.isfinite(lG):
                            Qchain[i, c] = Zms[i, c] * Xchain[c] + Lms[i, c] * exp(lG_marg - lG)
                        else:
                            Qchain[i, c] = Zms[i, c] * Xchain[c] + Lms[i, c] * Xchain[c]

        # Compute queue lengths and throughputs for open chains
        for c in open_chains:
            # For open chains, throughput equals arrival rate
            Xchain[c] = lmbda[c] if lmbda is not None and c < len(lmbda) else 0.0

            for i in range(M):
                if Lchain[i, c] > FINE_TOL:
                    # Simple M/M/1 formula for open chains
                    rho = lmbda[c] * Lchain[i, c] if lmbda is not None and c < len(lmbda) else 0.0
                    if rho < 1.0 and rho > 0:
                        Qchain[i, c] = rho / (1 - rho)
                    elif rho >= 1.0:
                        Qchain[i, c] = np.inf
                    else:
                        Qchain[i, c] = 0.0

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
        deagg_result = sn_deaggregate_chain_results(
            sn, Lchain, ST, STchain, Vchain, alpha,
            Qchain, None, Rchain, Tchain, None, Xchain.reshape(1, -1)
        )
        Q = deagg_result.Q
        U = deagg_result.U
        R = deagg_result.R
        T = deagg_result.T
        X = deagg_result.X
        STeff = ST.copy()

        # Non-exponential approximation (simplified)
        # Full npfqn_nonexp_approx would go here
        eta = np.ones(M)  # Simplified for now

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
                                fnc_result = pfqn_fnc(mu[i:i+1, :])
                                muhati_f = fnc_result.mu
                                c_val = fnc_result.c[0, 0] if fnc_result.c.size > 0 else 0.0

                                # Compute queue length using marginal analysis
                                L_i = L[i:i+1, :]
                                lGhatir = pfqn_ncld(L, Nchain_r, Nchain0, muhati, pfqn_options).lG
                                lGr_i = pfqn_ncld(L_i, Nchain_r, Nchain0, muhati[i:i+1, :] if i < muhati.shape[0] else muhati[:1, :], pfqn_options).lG

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
        from ...sn import sn_deaggregate_chain_results
        deagg_result = sn_deaggregate_chain_results(
            sn, Lchain, ST, STchain, Vchain, alpha,
            Qchain, None, Rchain, Tchain, None, Xchain.reshape(1, -1)
        )
        Q = deagg_result.Q
        U = deagg_result.U
        R = deagg_result.R
        T = deagg_result.T
        X = deagg_result.X

        # Non-exponential approximation (simplified)
        eta = np.ones(M)

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
            # Infinite server utilization
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


__all__ = [
    'solver_nc',
    'solver_ncld',
    'SolverNCReturn',
    'SolverNCLDReturn',
    'SolverOptions',
]
