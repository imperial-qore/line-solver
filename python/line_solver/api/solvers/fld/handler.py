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
    Build transition rate matrix W for phase-type fluid ODE.

    W encodes:
    1. Internal phase transitions (psi): from phase k to phase k' within same station-class
    2. Service completions with routing: from station i to station j via routing matrix

    For phase-type distributions:
    - D0 (psi block): internal transitions within phases
    - D1: completion rates (absorption) to exit the station
    - pie: initial phase probabilities upon arrival

    Args:
        sn: Network structure
        proc: Process information (station -> class -> [D0, D1])
        pie: Initial phase probabilities (station -> class -> probability vector)
        rt: Routing probability matrix (M*K x M*K), class-level routing

    Returns:
        Tuple of (W, A, B, psi matrices)
    """
    M = sn.nstations
    K = sn.nclasses

    # Get phases matrix
    phases = sn.phases if sn.phases is not None else np.ones((M, K))
    total_phases = int(np.sum(phases))

    # Compute starting index for each station-class
    q_indices = np.zeros((M, K), dtype=int)
    idx = 0
    for i in range(M):
        for r in range(K):
            q_indices[i, r] = idx
            idx += int(phases[i, r])

    # Initialize W matrix
    W = np.zeros((total_phases, total_phases))

    # Build W from D0 matrices (internal phase transitions)
    # D0[i,j] represents rate FROM phase i TO phase j
    # For ODE dx/dt = W@x, W[j,i] is rate from state i to state j
    # So we need to TRANSPOSE D0 for the off-diagonal terms
    # Diagonal terms stay as-is (departure rates are negative)
    for i in range(M):
        for r in range(K):
            nphases = int(phases[i, r])
            if nphases == 0:
                continue

            base_idx = q_indices[i, r]

            # Get rate for this station-class
            rate = sn.rates[i, r] if sn.rates is not None else 1.0

            # Skip disabled classes (NaN rates)
            if np.isnan(rate):
                continue

            # Check if proc has explicit D0 matrix
            # proc can be a dict or list
            has_proc_ir = False
            proc_ir = None
            if isinstance(proc, dict):
                if i in proc and r in proc[i]:
                    proc_ir = proc[i][r]
                    has_proc_ir = True
            elif isinstance(proc, (list, np.ndarray)):
                if i < len(proc) and proc[i] is not None:
                    if isinstance(proc[i], dict):
                        if r in proc[i]:
                            proc_ir = proc[i][r]
                            has_proc_ir = True
                    elif isinstance(proc[i], (list, np.ndarray)) and r < len(proc[i]):
                        proc_ir = proc[i][r]
                        has_proc_ir = True

            if has_proc_ir and proc_ir is not None:
                if isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 1:
                    D0 = np.asarray(proc_ir[0])
                elif isinstance(proc_ir, dict) and 'D0' in proc_ir:
                    D0 = np.asarray(proc_ir['D0'])
                else:
                    D0 = np.array([[-rate]])
            else:
                D0 = np.diag([-rate] * nphases)

            # Copy D0 with correct orientation:
            # - Diagonal: W[k,k] = D0[k,k] (departure rate, negative)
            # - Off-diagonal: W[j,i] = D0[i,j] (arrival rate at j from i)
            for k_from in range(min(nphases, D0.shape[0])):
                for k_to in range(min(nphases, D0.shape[1])):
                    if k_from == k_to:
                        # Diagonal - keep as-is
                        W[base_idx + k_from, base_idx + k_to] = D0[k_from, k_to]
                    else:
                        # Off-diagonal - transpose: D0[from,to] goes to W[to,from]
                        W[base_idx + k_to, base_idx + k_from] = D0[k_from, k_to]

    # Build psi for return (not actually used in ODE, just for compatibility)
    psi_blocks = []
    for i in range(M):
        for r in range(K):
            nphases = int(phases[i, r])
            rate = sn.rates[i, r] if sn.rates is not None else 1.0

            # Handle disabled classes (NaN rates) - use zero block
            if np.isnan(rate):
                psi_blocks.append(np.zeros((max(nphases, 1), max(nphases, 1))))
                continue

            # Check if proc has explicit D0 matrix
            has_proc_ir = False
            proc_ir = None
            if isinstance(proc, dict):
                if i in proc and r in proc[i]:
                    proc_ir = proc[i][r]
                    has_proc_ir = True
            elif isinstance(proc, (list, np.ndarray)):
                if i < len(proc) and proc[i] is not None:
                    if isinstance(proc[i], dict):
                        if r in proc[i]:
                            proc_ir = proc[i][r]
                            has_proc_ir = True
                    elif isinstance(proc[i], (list, np.ndarray)) and r < len(proc[i]):
                        proc_ir = proc[i][r]
                        has_proc_ir = True

            if nphases > 0 and has_proc_ir and proc_ir is not None:
                if isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 1:
                    D0 = np.asarray(proc_ir[0])
                    psi_blocks.append(D0)
                elif isinstance(proc_ir, dict) and 'D0' in proc_ir:
                    psi_blocks.append(np.asarray(proc_ir['D0']))
                else:
                    psi_blocks.append(np.array([[-rate]]))
            elif nphases > 0:
                psi_blocks.append(np.diag([-rate] * nphases))
            else:
                psi_blocks.append(np.zeros((1, 1)))

    psi = block_diag(*psi_blocks)

    # Add routing transitions: completion at station i routes to station j
    # For each source (i, r) and destination (j, s):
    #   W[dest_phase, src_phase] += completion_rate[src_phase] * P[i,r -> j,s] * pie[j,s][dest_phase]

    if rt is not None:
        for i in range(M):
            for r in range(K):
                nphases_src = int(phases[i, r])
                if nphases_src == 0:
                    continue

                # Skip disabled classes (NaN rates)
                src_rate = sn.rates[i, r] if sn.rates is not None else 1.0
                if np.isnan(src_rate):
                    continue

                src_idx = q_indices[i, r]

                # Get completion rates for source
                completion_rates = np.ones(nphases_src)
                # Check proc structure
                has_proc_ir = False
                proc_ir = None
                if isinstance(proc, dict):
                    if i in proc and r in proc[i]:
                        proc_ir = proc[i][r]
                        has_proc_ir = True
                elif isinstance(proc, (list, np.ndarray)):
                    if i < len(proc) and proc[i] is not None:
                        if isinstance(proc[i], dict):
                            if r in proc[i]:
                                proc_ir = proc[i][r]
                                has_proc_ir = True
                        elif isinstance(proc[i], (list, np.ndarray)) and r < len(proc[i]):
                            proc_ir = proc[i][r]
                            has_proc_ir = True

                if has_proc_ir and proc_ir is not None:
                    if isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 2:
                        D1 = np.asarray(proc_ir[1])
                        if D1.ndim == 1:
                            completion_rates = D1
                        else:
                            completion_rates = np.sum(D1, axis=1)
                    elif isinstance(proc_ir, dict) and 'D1' in proc_ir:
                        D1 = np.asarray(proc_ir['D1'])
                        completion_rates = np.sum(D1, axis=1) if D1.ndim == 2 else D1
                else:
                    completion_rates = np.full(nphases_src, src_rate)

                # Route to all destinations
                for j in range(M):
                    for s in range(K):
                        nphases_dst = int(phases[j, s])
                        if nphases_dst == 0:
                            continue

                        dst_idx = q_indices[j, s]

                        # Get routing probability
                        rt_idx_src = i * K + r
                        rt_idx_dst = j * K + s
                        if rt_idx_src < rt.shape[0] and rt_idx_dst < rt.shape[1]:
                            p_route = rt[rt_idx_src, rt_idx_dst]
                        else:
                            p_route = 0.0

                        if p_route <= 0:
                            continue

                        # Get initial phase probabilities at destination
                        pie_dst = np.zeros(nphases_dst)
                        pie_dst[0] = 1.0
                        if pie and j in pie and s in pie[j]:
                            pie_dst = np.asarray(pie[j][s]).flatten()

                        # Add routing contributions:
                        # W[dst_phase, src_phase] += completion_rate[src] * P * pie[dst]
                        for k_src in range(nphases_src):
                            for k_dst in range(nphases_dst):
                                rate = completion_rates[k_src] * p_route * pie_dst[k_dst]
                                W[dst_idx + k_dst, src_idx + k_src] += rate

    # Create placeholder A and B for compatibility
    A = np.eye(total_phases)
    B = np.eye(total_phases)

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
            # For multi-phase distributions, completion_rate[k] = D1 row sum for phase k
            completion_rates = np.ones(max(nphases, 1))
            if sn.proc and i in sn.proc and r in sn.proc[i]:
                proc_ir = sn.proc[i][r]
                if isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 2:
                    D1 = np.asarray(proc_ir[1])
                    if D1.ndim == 1:
                        completion_rates = D1
                    elif D1.ndim == 2:
                        completion_rates = np.sum(D1, axis=1)
                    if len(completion_rates) < nphases:
                        # Pad with last value
                        completion_rates = np.pad(
                            completion_rates,
                            (0, nphases - len(completion_rates)),
                            mode='edge'
                        )
            elif sn.rates is not None and i < sn.rates.shape[0] and r < sn.rates.shape[1]:
                # Use service rate from rates matrix for exponential service
                completion_rates = np.full(nphases, sn.rates[i, r])

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

                # Throughput contribution from this phase: T = completion_rate[k] * x[k]
                # This captures that only completion from certain phases contributes to throughput
                STC[i * K + r, state] = completion_rates[k] if k < len(completion_rates) else 0.0

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
    isSourceState: Optional[np.ndarray] = None,
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
        isSourceState: Boolean array indicating which states belong to Source stations

    Returns:
        dx/dt state derivative
    """
    x = np.maximum(x, 0)  # Ensure non-negative

    # Compute total queue at each state's station
    sum_x_Qa = SQ @ x + 1e-8  # Add FineTol for numerical stability (matching MATLAB GlobalConstants.FineTol)

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

        theta = x * ghat
        # For Source stations, theta = 0 to bypass Source in dynamics
        # (matching MATLAB where Source is excluded from state space)
        if isSourceState is not None:
            theta[isSourceState] = 0.0
        dxdt = W @ theta + ALambda.flatten()
    else:
        # Standard fluid constraint
        min_vals = np.minimum(sum_x_Qa, Sa.flatten())
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.where(sum_x_Qa > 1e-8, min_vals / sum_x_Qa, 1.0)
        ratio = np.nan_to_num(ratio, nan=1.0, posinf=1.0, neginf=0.0)

        theta = x * ratio
        # For Source stations, theta = 0 to bypass Source in dynamics
        # (matching MATLAB where Source is excluded from state space)
        if isSourceState is not None:
            theta[isSourceState] = 0.0
        dxdt = W @ theta + ALambda.flatten()

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
    from ....solvers.solver_fld.utils.phase_type import (
        prepare_phase_type_structures,
        extract_mu_phi_from_phase_type,
        compute_q_indices,
    )

    start_time = time.time()

    if options is None:
        options = SolverFLDOptions()

    M = sn.nstations
    K = sn.nclasses

    # Convert dict-based proc to phase-type matrices and compute phases
    proc_matrix, pie_dict, phases = prepare_phase_type_structures(sn)

    # Update sn with computed phases and proc for downstream use
    sn.phases = phases
    sn.proc = proc_matrix
    sn.pie = pie_dict

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

    # Build process-related structures (already converted to matrix form above)
    proc = proc_matrix
    pie = pie_dict
    # Extract station-to-station routing matrix P from stateful-indexed sn.rt
    # sn.rt is (nstateful*K x nstateful*K), indexed by stateful node, not station.
    # We need a (M*K x M*K) matrix indexed by station.
    # (Matching MATLAB solver_fluid_matrix.m lines 22-30)
    P = None
    if hasattr(sn, 'rt') and sn.rt is not None:
        station_indices = []
        for ist in range(M):
            isf = int(sn.stationToStateful[ist]) if hasattr(sn, 'stationToStateful') and sn.stationToStateful is not None else ist
            for r in range(K):
                station_indices.append(isf * K + r)
        P = sn.rt[np.ix_(station_indices, station_indices)].copy()

        # Remove Sink->Source feedback routing for open classes
        # In open networks, jobs exit at Sink and should not recirculate back to Source.
        # The routing matrix includes this feedback (added by getRoutingMatrix for
        # pseudo-closed network analysis), but it causes incorrect flow balance in the
        # fluid ODE formulation because arrivals are already accounted for via ALambda.
        # (Following MATLAB solver_fluid_matrix.m lines 37-56)
        for src_ist in range(M):
            # Check if this is a Source station
            node_idx = int(sn.stationToNode[src_ist]) if src_ist < len(sn.stationToNode) else src_ist
            is_source = False
            if node_idx < len(sn.nodetype) and sn.nodetype[node_idx] == NodeType.SOURCE:
                is_source = True

            if is_source:
                # This is a Source station - remove feedback routing TO it for open classes
                for r in range(K):
                    # Check if this is an open class with external arrivals
                    if sn.rates is not None and sn.rates[src_ist, r] > 0:
                        # Zero out routing TO this Source for this class from all other stations
                        src_col = src_ist * K + r  # Column index in P for (Source, class r)
                        for from_ist in range(M):
                            if from_ist != src_ist:  # Don't modify Source's own outgoing routing
                                for from_r in range(K):
                                    from_row = from_ist * K + from_r
                                    P[from_row, src_col] = 0  # Remove feedback to Source

    # Determine if this is a closed network
    is_closed = np.all(np.isfinite(sn.njobs))

    # Build transition matrix W
    if not proc and not pie:
        # No process info - use simple model (single phase per station-class)
        W, A, B, psi = _build_simple_W(sn, P)
        # _build_simple_W uses M*K dimensions (one phase per station-class)
        phases = np.ones((M, K))
        sn.phases = phases
    else:
        try:
            W, A, B, psi = _build_transition_matrix_W(sn, proc, pie, P)
        except Exception as e:
            # Fall back to simple exponential model (single phase per station-class)
            W, A, B, psi = _build_simple_W(sn, P)
            phases = np.ones((M, K))
            sn.phases = phases

    total_phases = int(np.sum(phases))

    # Handle dimension mismatch (when W was reduced due to disabled classes)
    if W.shape[0] != total_phases:
        # Rebuild phases to match W dimension (one phase per station-class)
        phases = np.ones((M, K))
        sn.phases = phases
        total_phases = int(np.sum(phases))

    # Build state mappings
    Qa, SQC, SUC, STC, SQ = _build_state_mappings(sn, phases, nservers, nservers_orig)

    # Ensure matrices match ODE dimension
    if SQ.shape[0] != W.shape[0]:
        # Resize to match
        n = W.shape[0]
        SQ = np.eye(n)
        SQC = np.eye(n)[:M * K, :] if M * K <= n else np.zeros((M * K, n))
        # Expand nservers from (M,) to (M*K, 1) by repeating each station's value K times
        nservers_expanded = np.repeat(nservers, K).reshape(-1, 1)
        SUC = SQC / nservers_expanded if SQC.shape[0] > 0 else SQC
        STC = SQC.copy()
        Qa = np.zeros((1, n))

    # Server capacity per state
    Sa = np.array([nservers[int(Qa[0, i]) if i < Qa.shape[1] else 0] for i in range(W.shape[0])])

    # External arrival rates - following MATLAB solver_fluid_matrix.m approach:
    # Build Alambda: arrivals go to QUEUE phases (where jobs route from Source), not Source phases
    # This matches dx = W' * theta + A * lambda where lambda represents arrivals INTO queues

    sched_dict = sn.sched if sn.sched else {}

    def _is_source_station(station_idx: int) -> bool:
        """Check if a station is a Source station."""
        node_idx = int(sn.stationToNode[station_idx]) if len(sn.stationToNode) > station_idx else station_idx
        if node_idx < len(sn.nodetype) and sn.nodetype[node_idx] == NodeType.SOURCE:
            return True
        if station_idx in sched_dict:
            sched_val = sched_dict[station_idx]
            if sched_val == SchedStrategy.EXT:
                return True
            if isinstance(sched_val, int) and sched_val == SchedStrategy.EXT.value:
                return True
        return False

    # First, identify Source stations and their arrival rates per class
    source_arrivals = np.zeros((M, K))
    for src_ist in range(M):
        if _is_source_station(src_ist):
            for r in range(K):
                if sn.rates is not None and not np.isnan(sn.rates[src_ist, r]) and sn.rates[src_ist, r] > 0:
                    source_arrivals[src_ist, r] = sn.rates[src_ist, r]

    # Compute phase indices for each (station, class)
    q_indices = np.zeros((M, K), dtype=int)
    idx = 0
    for i in range(M):
        for r in range(K):
            q_indices[i, r] = idx
            idx += int(phases[i, r])

    # Build ALambda: arrivals go to QUEUE phases (not Source phases)
    ALambda = np.zeros((W.shape[0], 1))
    state = 0
    for ist in range(M):
        for r in range(K):
            nphases_ir = int(phases[ist, r])
            if nphases_ir > 0:
                if _is_source_station(ist):
                    # Source station: do NOT add arrivals here (arrivals go to downstream queues)
                    state += nphases_ir
                else:
                    # Queue station: check if it receives arrivals from any Source
                    arrival_rate_to_queue = 0.0
                    for src_ist in range(M):
                        if source_arrivals[src_ist, r] > 0:
                            # Get routing probability from Source to this queue for this class
                            # P is station-indexed (M*K x M*K)
                            src_row = src_ist * K + r
                            queue_col = ist * K + r
                            if P is not None and src_row < P.shape[0] and queue_col < P.shape[1]:
                                routing_prob = P[src_row, queue_col]
                                arrival_rate_to_queue += source_arrivals[src_ist, r] * routing_prob

                    if arrival_rate_to_queue > 0:
                        # Apply arrivals according to entrance probability pie
                        pie_ir = np.zeros(nphases_ir)
                        pie_ir[0] = 1.0  # Default: all arrivals go to first phase
                        if pie and ist in pie and r in pie[ist]:
                            pie_arr = np.asarray(pie[ist][r]).flatten()
                            pie_ir[:min(len(pie_arr), nphases_ir)] = pie_arr[:min(len(pie_arr), nphases_ir)]

                        for k in range(nphases_ir):
                            if state < ALambda.shape[0]:
                                ALambda[state, 0] = pie_ir[k] * arrival_rate_to_queue
                            state += 1
                    else:
                        state += nphases_ir
            else:
                # Disabled class - add placeholder state
                state += 1

    # Initial state - use sn.state to initialize ODE state vector
    # This matches MATLAB's solver_fluid_initsol which reads model state
    # (typically all closed-class jobs at reference station, zero elsewhere)

    x0 = np.zeros(W.shape[0])

    if options.init_sol is not None and len(options.init_sol) == W.shape[0]:
        # Use explicitly provided initial solution
        x0 = np.array(options.init_sol, dtype=float)
    elif hasattr(sn, 'state') and sn.state is not None and len(sn.state) > 0:
        # Build x0 from sn.state (matching MATLAB solver_fluid_initsol)
        # sn.state[isf] contains per-class marginal job counts
        idx = 0
        for ist in range(M):
            isf = int(sn.stationToStateful[ist]) if hasattr(sn, 'stationToStateful') and sn.stationToStateful is not None else ist
            for k in range(K):
                nphases_ik = int(phases[ist, k])
                if nphases_ik == 0:
                    continue

                # Get number of jobs of class k at station ist from sn.state
                n_jobs_at_station = 0.0
                if isf < len(sn.state) and sn.state[isf] is not None:
                    state_vec = np.asarray(sn.state[isf]).flatten()
                    if k < len(state_vec):
                        n_jobs_at_station = float(state_vec[k])
                    elif len(state_vec) == 1 and K == 1:
                        n_jobs_at_station = float(state_vec[0])

                # Check if this is a Source station - set state to 0
                node_idx_ist = int(sn.stationToNode[ist]) if ist < len(sn.stationToNode) else ist
                is_source = (node_idx_ist < len(sn.nodetype) and sn.nodetype[node_idx_ist] == NodeType.SOURCE)
                sched_ist = sched_dict.get(ist) if sched_dict else None
                if sched_ist == SchedStrategy.EXT:
                    is_source = True

                if is_source:
                    # Source stations: no mass (arrivals via ALambda)
                    idx += nphases_ik
                    continue

                if np.isnan(n_jobs_at_station):
                    n_jobs_at_station = 0.0

                # Place all jobs in first phase (matching MATLAB: jobs in waiting buffer
                # are re-started from phase 1)
                x0[idx] = n_jobs_at_station
                idx += nphases_ik
    else:
        # Fallback: distribute evenly (legacy behavior)
        match = np.zeros((M, K))
        if P is not None:
            for ist in range(M):
                for k in range(K):
                    dst_idx = ist * K + k
                    if dst_idx < P.shape[1]:
                        match[ist, k] = 1 if np.sum(P[:, dst_idx]) > 0 else 0

        njobs_flat = sn.njobs.flatten() if sn.njobs is not None else np.zeros(K)
        assigned = np.zeros(K)

        idx = 0
        for ist in range(M):
            is_source_station = False
            node_idx = int(sn.stationToNode[ist]) if ist < len(sn.stationToNode) else ist
            if node_idx < len(sn.nodetype) and sn.nodetype[node_idx] == NodeType.SOURCE:
                is_source_station = True
            sched_ist = sched_dict.get(ist) if sched_dict else None
            if sched_ist == SchedStrategy.EXT:
                is_source_station = True

            for k in range(K):
                nphases_ik = int(phases[ist, k])
                if nphases_ik == 0:
                    continue

                if match[ist, k] > 0:
                    if k < len(njobs_flat) and np.isinf(njobs_flat[k]):
                        if is_source_station:
                            to_assign = 1
                        else:
                            to_assign = 0
                    elif k < len(njobs_flat) and np.isfinite(njobs_flat[k]):
                        n_jobs = njobs_flat[k]
                        num_stations_serving = int(np.sum(match[:, k]))
                        if num_stations_serving > 0:
                            to_assign = int(n_jobs // num_stations_serving)
                            remaining_stations = int(np.sum(match[ist+1:, k]))
                            if remaining_stations == 0:
                                to_assign = int(n_jobs - assigned[k])
                        else:
                            to_assign = 0
                    else:
                        to_assign = 0

                    x0[idx] = to_assign
                    assigned[k] += to_assign

                idx += nphases_ik

    # P-star values for smoothing
    pstar = None
    if options.pstar and len(options.pstar) > 0:
        # Expand pstar to match state dimension
        pstar = np.zeros(W.shape[0])
        idx = 0
        n_states = W.shape[0]
        for i in range(M):
            pstar_i = options.pstar[i] if i < len(options.pstar) else 10.0
            for r in range(K):
                nphases = int(phases[i, r])
                for k in range(nphases):
                    if idx < n_states:
                        pstar[idx] = pstar_i
                    idx += 1

    # Identify Source station states (EXT scheduler)
    # For Source stations, theta should be 0.0 to effectively bypass Source in dynamics.
    # This matches the MATLAB implementation where Source is excluded from state space.
    # Arrivals are injected directly into queue phases via ALambda.
    isSourceState = np.zeros(W.shape[0], dtype=bool)
    state_idx = 0
    for i in range(M):
        node_idx = int(sn.stationToNode[i]) if i < len(sn.stationToNode) else i
        is_source = (node_idx < len(sn.nodetype) and sn.nodetype[node_idx] == NodeType.SOURCE)
        for r in range(K):
            nphases_ir = int(phases[i, r])
            for k in range(nphases_ir):
                if is_source:
                    isSourceState[state_idx] = True
                    x0[state_idx] = 0.0  # Initialize Source phases to 0 (no mass at Source)
                state_idx += 1

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
                lambda t, x: _fluid_ode(t, x, W, SQ, Sa, ALambda, pstar, isSourceState),
                [T_start, T_end],
                x0,
                method=method,
                rtol=options.tol,
                atol=options.tol,
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
    sum_x_Qa = SQ @ x_final + 1e-8
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
        sum_x_step = SQ @ x_step + 1e-8
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


def _build_simple_W(sn: NetworkStruct, P: np.ndarray = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Build simple transition matrix W for exponential service.

    Fallback when proc/pie structures are not available.
    Uses standard fluid equations for closed queueing networks.

    Args:
        sn: Network structure
        P: Station-indexed routing matrix (M*K x M*K). If None, uses visits-based routing.
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

                # Routing - use station-indexed P if available, otherwise use visits
                if P is not None:
                    for j in range(M):
                        for s in range(K):
                            idx_j = j * K + s
                            src_idx = i * K + r
                            dst_idx = j * K + s
                            if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                                p_ij = P[src_idx, dst_idx]
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
