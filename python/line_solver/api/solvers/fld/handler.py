"""
FLD Solver handler.

Native Python implementation of FLD (Fluid/Mean-Field Approximation) solver handler
that orchestrates ODE-based fluid analysis of queueing networks.

Port from:


"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple, Any
import time
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag

from ...sn import (
    NetworkStruct,
    SchedStrategy,
    NodeType,
)


@dataclass
class SolverFLDOptions:
    """Options for FLD solver."""
    method: str = 'matrix'
    tol: float = 1e-6
    verbose: bool = False
    stiff: bool = True
    iter_max: int = 100
    timespan: Tuple[float, float] = (0.0, float('inf'))
    init_sol: Optional[np.ndarray] = None
    pstar: Optional[List[float]] = None  # P-norm smoothing values per station
    num_cdf_pts: int = 200  # Number of points for CDF computation


@dataclass
class SolverFLDReturn:
    """
    Result of FLD solver handler.

    Attributes:
        Q: Mean queue lengths (M x K)
        U: Utilizations (M x K)
        R: Response times (M x K)
        T: Throughputs (M x K)
        C: Cycle times (1 x K)
        X: System throughputs (1 x K)
        Qt: Transient queue lengths (list of (M x K) arrays per time point)
        Ut: Transient utilizations
        Tt: Transient throughputs
        t: Time vector
        odeStateVec: Final ODE state vector
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
    Qt: Optional[List[np.ndarray]] = None
    Ut: Optional[List[np.ndarray]] = None
    Tt: Optional[List[np.ndarray]] = None
    t: Optional[np.ndarray] = None
    odeStateVec: Optional[np.ndarray] = None
    runtime: float = 0.0
    method: str = "matrix"
    it: int = 0


def _build_transition_matrix_W(
    sn: NetworkStruct,
    proc: Dict,
    pie: Dict,
    rt: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Build transition rate matrix W as per Ruuskanen et al., PEVA 151 (2021).

    W = psi + B * P * A' where:
    - psi: Block diagonal of internal phase transitions
    - A: Block diagonal of initial phase probabilities (transposed)
    - B: Block diagonal of completion rate column vectors
    - P: Routing probability matrix (rt)

    Args:
        sn: Network structure
        proc: Process information (station -> class -> [psi_matrix, completion_matrix])
        pie: Initial phase probabilities (station -> class -> probability vector)
        rt: Routing probability matrix

    Returns:
        Tuple of (W, A, B, psi matrices)
    """
    M = sn.nstations
    K = sn.nclasses

    # Build block diagonal matrices psi, A, B
    psi_blocks = []
    A_blocks = []
    B_blocks = []

    stations = list(range(M))
    jobclasses = list(range(K))

    # Get phases matrix
    phases = sn.phases if sn.phases is not None else np.ones((M, K))

    for i in stations:
        for r in jobclasses:
            nphases = int(phases[i, r]) if phases[i, r] > 0 else 0

            if nphases == 0 or phases[i, r] == 0:
                # Disabled class-station pair
                psi_blocks.append(np.zeros((1, 1)))
                A_blocks.append(np.array([[np.nan]]))
                B_blocks.append(np.zeros((1, 1)))
            else:
                # Get process info
                if proc and i in proc and r in proc[i]:
                    proc_ir = proc[i][r]
                    if isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 2:
                        # proc[i][r] = [psi_matrix, completion_matrix]
                        psi_ir = np.asarray(proc_ir[0])
                        completion_ir = np.asarray(proc_ir[1])
                    else:
                        # Single matrix - use as psi
                        psi_ir = np.asarray(proc_ir)
                        completion_ir = -np.sum(psi_ir, axis=1, keepdims=True)
                else:
                    # Default: simple exponential service
                    rate = sn.rates[i, r] if sn.rates is not None else 1.0
                    psi_ir = np.array([[-rate]])
                    completion_ir = np.array([[rate]])

                psi_blocks.append(psi_ir)

                # Completion rate column vector - sum of completion rates per phase
                B_ir = completion_ir.sum(axis=1, keepdims=True) if completion_ir.ndim > 1 else completion_ir.reshape(-1, 1)
                B_blocks.append(B_ir)

                # Initial phase probability
                if pie and i in pie and r in pie[i]:
                    pie_ir = np.asarray(pie[i][r]).flatten()
                else:
                    # Default: start in first phase
                    pie_ir = np.zeros(nphases)
                    pie_ir[0] = 1.0
                A_blocks.append(pie_ir.reshape(-1, 1))

    # Build block diagonal matrices
    psi = block_diag(*psi_blocks)

    # A is block diagonal of initial phase probabilities (transposed)
    A = block_diag(*[block.T for block in A_blocks])

    # B is block diagonal of completion rate vectors
    B = block_diag(*[block for block in B_blocks])

    # Compute W = psi + (B * P * A')^T
    # W represents the complete graph of transition rates
    # W[j,i] = rate of fluid flowing from state i to state j
    # The routing term B @ rt @ A.T gives entry [i,j] for rate from i to j,
    # so we transpose to get W[j,i] = rate from i to j
    if rt is not None:
        W = psi + (B @ rt @ A.T).T
    else:
        W = psi

    # Remove disabled transitions (NaN columns/rows)
    if np.any(np.isnan(W)):
        valid_cols = ~np.isnan(W.sum(axis=0))
        valid_rows = ~np.isnan(W.sum(axis=1))
        valid = valid_cols & valid_rows
        W = W[np.ix_(valid, valid)]
        A = A[np.ix_(valid, valid)]
        B = B[np.ix_(valid, valid)]
        psi = psi[np.ix_(valid, valid)]

    return W, A, B, psi


def _build_state_mappings(
    sn: NetworkStruct,
    phases: np.ndarray,
    nservers: np.ndarray,
    nservers_orig: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Build state mapping matrices for converting ODE state to performance metrics.

    Args:
        sn: Network structure
        phases: (M x K) number of phases per station-class
        nservers: (M,) number of servers per station (with inf replaced by population)
        nservers_orig: (M,) original number of servers (may contain inf for IS nodes)

    Returns:
        Tuple of (Qa, SQC, SUC, STC, SQ) where:
        - Qa: (1, total_phases) state -> station mapping
        - SQC: (M*K, total_phases) state -> queue length
        - SUC: (M*K, total_phases) state -> utilization
        - STC: (M*K, total_phases) state -> throughput
        - SQ: (total_phases, total_phases) state -> total queue at station
    """
    M = sn.nstations
    K = sn.nclasses

    total_phases = int(np.sum(phases))

    Qa = np.zeros((1, total_phases))
    SQC = np.zeros((M * K, total_phases))
    SUC = np.zeros((M * K, total_phases))
    STC = np.zeros((M * K, total_phases))
    SQ = np.zeros((total_phases, total_phases))

    state = 0
    for i in range(M):
        for r in range(K):
            nphases = int(phases[i, r])

            # Get completion rates for throughput calculation
            # For exponential service, completion_rate = μ (service rate)
            completion_rate = 1.0
            if sn.proc and i in sn.proc and r in sn.proc[i]:
                proc_ir = sn.proc[i][r]
                if isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 2:
                    completion_ir = np.asarray(proc_ir[1])
                    if completion_ir.size > 0:
                        completion_rate = np.sum(completion_ir, axis=1) if completion_ir.ndim > 1 else completion_ir
            elif sn.rates is not None and i < sn.rates.shape[0] and r < sn.rates.shape[1]:
                # Use service rate from rates matrix for exponential service
                completion_rate = sn.rates[i, r]

            for k in range(nphases):
                Qa[0, state] = i
                SQC[i * K + r, state] = 1.0

                # For IS (Infinite Server) nodes, utilization = queue length (all jobs in service)
                # For finite-server nodes, utilization = jobs in service / servers
                if np.isinf(nservers_orig[i]):
                    # IS node: U = Q (utilization coefficient is 1.0)
                    SUC[i * K + r, state] = 1.0
                else:
                    SUC[i * K + r, state] = 1.0 / nservers[i] if nservers[i] > 0 else 0.0

                # Throughput contribution from this phase: T = μ * theta
                if isinstance(completion_rate, np.ndarray):
                    STC[i * K + r, state] = completion_rate[k] if k < len(completion_rate) else completion_rate[-1]
                else:
                    STC[i * K + r, state] = completion_rate

                state += 1

    # Build SQ matrix - maps state to total queue at each station
    state = 0
    for i in range(M):
        for r in range(K):
            nphases = int(phases[i, r])
            for k in range(nphases):
                # Mark all phases at the same station
                for col in range(total_phases):
                    if Qa[0, col] == i:
                        SQ[state, col] = 1.0
                state += 1

    return Qa, SQC, SUC, STC, SQ


def _fluid_ode(
    t: float,
    x: np.ndarray,
    W: np.ndarray,
    SQ: np.ndarray,
    Sa: np.ndarray,
    ALambda: np.ndarray,
    pstar: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Fluid ODE right-hand side for queueing network fluid analysis.

    dx/dt = W * (x .* g(x)) + ALambda

    where g(x) is the server constraint:
    - min(S, sum_station(x)) / sum_station(x) (without smoothing)
    - p-norm smoothed constraint function (with smoothing)

    The W matrix encodes both service rates and routing:
    - W[i,i] = -mu_i (departure rate from station i)
    - W[i,j] = mu_j * P_{j->i} (arrival rate at i from j)

    Args:
        t: Time
        x: State vector (queue lengths)
        W: Transition rate matrix
        SQ: State-to-station queue mapping
        Sa: Server capacity per state
        ALambda: External arrival rates
        pstar: P-norm smoothing parameters

    Returns:
        dx/dt state derivative
    """
    x = np.maximum(x, 0)  # Ensure non-negative

    # Compute total queue at each state's station
    sum_x_Qa = SQ @ x + 1e-14  # Add small value for numerical stability

    if pstar is not None and len(pstar) > 0:
        # P-norm smoothed constraint as per Ruuskanen et al.
        ghat = np.zeros_like(x)
        for i in range(len(x)):
            x_val = sum_x_Qa[i]
            c_val = Sa[i]
            p_val = pstar[i] if i < len(pstar) else pstar[-1]

            if p_val > 0 and c_val > 0:
                ghat_val = 1.0 / np.power(1 + np.power(x_val / c_val, p_val), 1.0 / p_val)
                if np.isnan(ghat_val):
                    ghat[i] = 0.0
                else:
                    ghat[i] = ghat_val
            else:
                ghat[i] = 1.0

        dxdt = W @ (x * ghat) + ALambda.flatten()
    else:
        # Standard fluid constraint
        min_vals = np.minimum(sum_x_Qa, Sa.flatten())
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.where(sum_x_Qa > 1e-14, min_vals / sum_x_Qa, 1.0)
        ratio = np.nan_to_num(ratio, nan=1.0, posinf=1.0, neginf=0.0)

        dxdt = W @ (x * ratio) + ALambda.flatten()

    return dxdt


def solver_fld(
    sn: NetworkStruct,
    options: Optional[SolverFLDOptions] = None
) -> SolverFLDReturn:
    """
    FLD solver handler using matrix method.

    Performs Fluid/Mean-Field analysis by:
    1. Building transition rate matrix W from process representations
    2. Setting up initial state from network configuration
    3. Integrating ODEs to steady state or specified timespan
    4. Extracting performance metrics from ODE solution

    Args:
        sn: Network structure with proc, pie, rt fields
        options: Solver options

    Returns:
        SolverFLDReturn with all performance metrics

    Raises:
        RuntimeError: For unsupported configurations
    """
    start_time = time.time()

    if options is None:
        options = SolverFLDOptions()

    M = sn.nstations
    K = sn.nclasses

    # Get phases matrix
    if sn.phases is not None:
        phases = sn.phases.copy()
    else:
        # Default to 1 phase (exponential)
        phases = np.ones((M, K))

    # Get server counts
    if sn.nservers is not None and len(sn.nservers.flatten()) > 0:
        nservers_orig = sn.nservers.flatten().copy()  # Keep original for IS detection
        nservers = nservers_orig.copy()
        # Replace inf with total population for delay stations (for ODE numerics)
        njobs_flat = np.asarray(sn.njobs).flatten()
        total_pop = sn.nclosedjobs if sn.nclosedjobs > 0 else np.sum(njobs_flat[np.isfinite(njobs_flat)])
        if total_pop == 0:
            total_pop = 1000  # Default for open networks
        nservers = np.where(np.isinf(nservers), total_pop, nservers)
    else:
        nservers = np.ones(M)
        nservers_orig = np.ones(M)

    # Build process-related structures
    proc = sn.proc if hasattr(sn, 'proc') and sn.proc else {}
    pie = sn.pie if hasattr(sn, 'pie') and sn.pie else {}
    rt = sn.rt if hasattr(sn, 'rt') and sn.rt is not None else None

    # Determine if this is a closed network
    is_closed = np.all(np.isfinite(sn.njobs))

    # Build transition matrix W
    if not proc and not pie:
        # No process info - use simple model
        W, A, B, psi = _build_simple_W(sn)
    else:
        try:
            W, A, B, psi = _build_transition_matrix_W(sn, proc, pie, rt)
        except Exception as e:
            # Fall back to simple exponential model
            W, A, B, psi = _build_simple_W(sn)

    total_phases = int(np.sum(phases))

    # Handle dimension mismatch (when W was reduced due to disabled classes)
    if W.shape[0] != total_phases:
        # Rebuild with only valid phases
        valid_phases = W.shape[0]
        total_phases = valid_phases

    # Build state mappings
    Qa, SQC, SUC, STC, SQ = _build_state_mappings(sn, phases, nservers, nservers_orig)

    # Ensure matrices match ODE dimension
    if SQ.shape[0] != W.shape[0]:
        # Resize to match
        n = W.shape[0]
        SQ = np.eye(n)
        SQC = np.eye(n)[:M * K, :] if M * K <= n else np.zeros((M * K, n))
        SUC = SQC / nservers.reshape(-1, 1)[:M * K, :1] if SQC.shape[0] > 0 else SQC
        STC = SQC.copy()
        Qa = np.zeros((1, n))

    # Server capacity per state
    Sa = np.array([nservers[int(Qa[0, i]) if i < Qa.shape[1] else 0] for i in range(W.shape[0])])

    # External arrival rates
    lambda_arr = np.zeros(M * K)
    sched_dict = sn.sched if sn.sched else {}

    for i in range(M):
        station_sched = sched_dict.get(i)
        node_idx = int(sn.stationToNode[i]) if len(sn.stationToNode) > i else i

        if node_idx < len(sn.nodetype) and sn.nodetype[node_idx] == NodeType.SOURCE:
            for r in range(K):
                if sn.rates is not None and sn.rates[i, r] > 0:
                    lambda_arr[i * K + r] = sn.rates[i, r]

    # A * lambda for arrivals
    if A.shape[0] > 0 and len(lambda_arr) > 0:
        ALambda = A @ lambda_arr.reshape(-1, 1) if A.shape[1] == len(lambda_arr) else np.zeros((A.shape[0], 1))
    else:
        ALambda = np.zeros((W.shape[0], 1))

    # Initial state - distribute jobs across stations
    if options.init_sol is not None and len(options.init_sol) > 0:
        x0 = options.init_sol.flatten()
        if len(x0) != W.shape[0]:
            x0 = np.zeros(W.shape[0])
            # Distribute initial jobs evenly
            njobs_flat_init = np.asarray(sn.njobs).flatten()
            total_jobs = np.sum(njobs_flat_init[np.isfinite(njobs_flat_init)])
            if total_jobs > 0:
                x0[:] = total_jobs / len(x0)
    else:
        # Default: place jobs at reference station or first station
        x0 = np.zeros(W.shape[0])

        # For each class, place jobs at reference station
        for r in range(K):
            if r < len(sn.njobs.flatten()) and np.isfinite(sn.njobs.flatten()[r]):
                n_jobs = sn.njobs.flatten()[r]
                # Determine which station to place jobs
                ref_i = int(sn.refstat[r]) if r < len(sn.refstat) else 0
                ref_i = min(ref_i, M - 1)  # Clamp to valid range

                # Find the state index for this station-class
                idx = 0
                for i in range(M):
                    for s in range(K):
                        if i == ref_i and s == r:
                            x0[idx] = n_jobs
                        idx += int(phases[i, s])

    # P-star values for smoothing
    pstar = None
    if options.pstar and len(options.pstar) > 0:
        # Expand pstar to match state dimension
        pstar = np.zeros(W.shape[0])
        idx = 0
        for i in range(M):
            pstar_i = options.pstar[i] if i < len(options.pstar) else 10.0
            for r in range(K):
                nphases = int(phases[i, r])
                for k in range(nphases):
                    pstar[idx] = pstar_i
                    idx += 1

    # Time span
    min_rate = np.abs(W[W != 0]).min() if np.any(W != 0) else 1.0
    T_end = min(options.timespan[1], abs(10 * options.iter_max / min_rate))
    T_start = options.timespan[0] if np.isfinite(options.timespan[0]) else 0.0

    # Check if W is essentially zero (equilibrium case)
    W_norm = np.linalg.norm(W)
    if W_norm < 1e-10:
        # W is essentially zero - system is at equilibrium
        # Return initial state as the solution
        t_vec = np.array([T_start, T_end])
        x_vec = np.vstack([x0, x0])  # Constant solution
    else:
        # Solve ODE
        try:
            if options.stiff:
                method = 'LSODA'
            else:
                method = 'RK45'

            sol = solve_ivp(
                lambda t, x: _fluid_ode(t, x, W, SQ, Sa, ALambda, pstar),
                [T_start, T_end],
                x0,
                method=method,
                rtol=options.tol,
                atol=options.tol * 1e-3,
                dense_output=True,
            )

            t_vec = sol.t
            x_vec = sol.y.T  # Shape: (n_times, n_states)

        except Exception as e:
            # Return empty result on failure
            result = SolverFLDReturn(
                Q=np.full((M, K), np.nan),
                U=np.full((M, K), np.nan),
                R=np.full((M, K), np.nan),
                T=np.full((M, K), np.nan),
                C=np.full((1, K), np.nan),
                X=np.full((1, K), np.nan),
                t=np.array([]),
                odeStateVec=np.array([]),
                runtime=time.time() - start_time,
                method=options.method,
                it=0
            )
            return result

    # Extract final state
    x_final = x_vec[-1, :]
    x_final = np.maximum(x_final, 0)  # Ensure non-negative

    # Identify Source and Sink stations (they don't hold jobs)
    source_stations = set()
    sink_stations = set()
    for i in range(M):
        node_idx = int(sn.stationToNode[i]) if i < len(sn.stationToNode) else i
        if node_idx < len(sn.nodetype):
            if sn.nodetype[node_idx] == NodeType.SOURCE:
                source_stations.add(i)
            elif sn.nodetype[node_idx] == NodeType.SINK:
                sink_stations.add(i)

    # Compute performance metrics from final state
    Q = np.zeros((M, K))
    U = np.zeros((M, K))
    T = np.zeros((M, K))
    R = np.zeros((M, K))

    # Compute theta (effective service rate fraction)
    sum_x_Qa = SQ @ x_final + 1e-14
    theta = x_final.copy()
    for phase in range(len(x_final)):
        station = int(Qa[0, phase]) if phase < Qa.shape[1] else 0
        theta[phase] = x_final[phase] / sum_x_Qa[phase] * min(Sa[phase], sum_x_Qa[phase])

    # Queue lengths
    if SQC.shape[1] == len(x_final):
        QN_flat = SQC @ x_final
        for i in range(M):
            for r in range(K):
                if i * K + r < len(QN_flat):
                    Q[i, r] = QN_flat[i * K + r]

    # Utilizations
    if SUC.shape[1] == len(theta):
        UN_flat = SUC @ theta
        for i in range(M):
            for r in range(K):
                if i * K + r < len(UN_flat):
                    U[i, r] = UN_flat[i * K + r]

    # Throughputs
    if STC.shape[1] == len(theta):
        TN_flat = STC @ theta
        for i in range(M):
            for r in range(K):
                if i * K + r < len(TN_flat):
                    T[i, r] = TN_flat[i * K + r]

    # Response times via Little's Law
    with np.errstate(divide='ignore', invalid='ignore'):
        R = np.where(T > 1e-14, Q / T, 0.0)
    R = np.nan_to_num(R, nan=0.0, posinf=0.0, neginf=0.0)

    # Override Source and Sink station metrics
    # Source: Q=0, U=0, R=0, T=arrival_rate
    # Sink: Q=0, U=0, R=0, T=throughput flowing into sink
    for i in source_stations:
        for r in range(K):
            Q[i, r] = 0.0
            U[i, r] = 0.0
            R[i, r] = 0.0
            # Throughput is the arrival rate
            if sn.rates is not None and i < sn.rates.shape[0] and r < sn.rates.shape[1]:
                T[i, r] = sn.rates[i, r]

    for i in sink_stations:
        for r in range(K):
            Q[i, r] = 0.0
            U[i, r] = 0.0
            R[i, r] = 0.0
            # Throughput at sink is the arrival rate (in steady state, same as source)

    # System throughput (per class)
    X = np.zeros((1, K))
    for r in range(K):
        ref_stat = int(sn.refstat[r]) if r < len(sn.refstat) else 0
        if ref_stat < M:
            X[0, r] = T[ref_stat, r]

    # Cycle times
    C = np.zeros((1, K))
    for r in range(K):
        if X[0, r] > 1e-14:
            C[0, r] = np.sum(R[:, r])

    # Build transient data
    Qt = [[np.zeros(len(t_vec)) for _ in range(K)] for _ in range(M)]
    Ut = [[np.zeros(len(t_vec)) for _ in range(K)] for _ in range(M)]
    Tt = [[np.zeros(len(t_vec)) for _ in range(K)] for _ in range(M)]

    for step in range(len(t_vec)):
        x_step = np.maximum(x_vec[step, :], 0)
        sum_x_step = SQ @ x_step + 1e-14
        theta_step = x_step.copy()
        for phase in range(len(x_step)):
            station = int(Qa[0, phase]) if phase < Qa.shape[1] else 0
            theta_step[phase] = x_step[phase] / sum_x_step[phase] * min(Sa[phase], sum_x_step[phase])

        if SQC.shape[1] == len(x_step):
            QN_step = SQC @ x_step
            UN_step = SUC @ theta_step
            TN_step = STC @ theta_step

            for i in range(M):
                for r in range(K):
                    idx = i * K + r
                    if idx < len(QN_step):
                        Qt[i][r][step] = QN_step[idx]
                        Ut[i][r][step] = UN_step[idx]
                        Tt[i][r][step] = TN_step[idx]

    result = SolverFLDReturn(
        Q=Q,
        U=U,
        R=R,
        T=T,
        C=C,
        X=X,
        Qt=Qt,
        Ut=Ut,
        Tt=Tt,
        t=t_vec,
        odeStateVec=x_final,
        runtime=time.time() - start_time,
        method=options.method,
        it=len(t_vec)
    )

    return result


def _build_simple_W(sn: NetworkStruct) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Build simple transition matrix W for exponential service.

    Fallback when proc/pie structures are not available.
    Uses standard fluid equations for closed queueing networks.
    """
    M = sn.nstations
    K = sn.nclasses

    # For fluid analysis of closed networks, we build a simple generator
    # W[i,i] = -mu_i (departure rate from state i)
    # W[i,j] = mu_j * V_j / sum_k(V_k) (arrival rate to state i from state j)
    n = M * K
    W = np.zeros((n, n))
    A = np.eye(n)
    B = np.eye(n)
    psi = np.zeros((n, n))

    # Get visits for routing
    V = np.ones((M, K))
    if sn.visits:
        for chain_id, chain_visits in sn.visits.items():
            if chain_visits is not None and chain_visits.shape == V.shape:
                V = chain_visits.copy()
                break

    for i in range(M):
        for r in range(K):
            idx = i * K + r
            rate = sn.rates[i, r] if sn.rates is not None and i < sn.rates.shape[0] else 0.0

            if rate > 0:
                # Departure rate
                psi[idx, idx] = -rate

                # Routing - use rt if available, otherwise use visits
                if sn.rt is not None:
                    for j in range(M):
                        for s in range(K):
                            idx_j = j * K + s
                            src_idx = i * K + r
                            dst_idx = j * K + s
                            if src_idx < sn.rt.shape[0] and dst_idx < sn.rt.shape[1]:
                                p_ij = sn.rt[src_idx, dst_idx]
                                if p_ij > 0:
                                    W[idx_j, idx] += rate * p_ij  # Note: W[j,i] for arrival at j from i
                else:
                    # Use visits-based routing (cyclic)
                    # Next station is (i+1) mod M
                    j = (i + 1) % M
                    idx_j = j * K + r
                    W[idx_j, idx] += rate  # All jobs go to next station

    W = psi + W

    return W, A, B, psi


def _build_closed_network_W(sn: NetworkStruct) -> np.ndarray:
    """
    Build W matrix specifically for closed networks using standard fluid model.

    For a closed network with M stations and K classes, the fluid equations are:
    dx_i/dt = -mu_i * min(1, S_i / sum_x_i) * x_i + sum_j(mu_j * P_{ji} * x_j * min(1, S_j / sum_x_j))

    This function returns W for: dx/dt = W * theta(x)
    where theta(x) = x * min(1, S / sum_x)
    """
    M = sn.nstations
    K = sn.nclasses
    n = M * K

    W = np.zeros((n, n))

    # Get visits for routing
    V = np.ones((M, K))
    if sn.visits:
        for chain_id, chain_visits in sn.visits.items():
            if chain_visits is not None and chain_visits.shape == V.shape:
                V = chain_visits.copy()
                break

    for i in range(M):
        for r in range(K):
            idx = i * K + r
            rate = sn.rates[i, r] if sn.rates is not None and i < sn.rates.shape[0] else 0.0

            if rate > 0:
                # Diagonal: departure rate from this station-class
                W[idx, idx] = -rate

                # Off-diagonal: arrivals from other stations
                if sn.rt is not None and sn.rt.shape[0] > idx and sn.rt.shape[1] > idx:
                    for j in range(M):
                        for s in range(K):
                            src_idx = j * K + s
                            if src_idx < sn.rt.shape[0] and idx < sn.rt.shape[1]:
                                p_ji = sn.rt[src_idx, idx]  # Prob from j,s to i,r
                                rate_j = sn.rates[j, s] if j < sn.rates.shape[0] else 0.0
                                if p_ji > 0 and rate_j > 0:
                                    W[idx, src_idx] += rate_j * p_ji
                else:
                    # Use cyclic routing based on visits
                    # Previous station sends jobs here
                    j = (i - 1) % M
                    src_idx = j * K + r
                    rate_j = sn.rates[j, r] if j < sn.rates.shape[0] else 0.0
                    if rate_j > 0:
                        W[idx, src_idx] += rate_j

    return W


__all__ = [
    'solver_fld',
    'SolverFLDReturn',
    'SolverFLDOptions',
]
