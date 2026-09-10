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
from ....constants import GlobalConstants


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
    odemaxstep: Optional[float] = None  # ODE solver max step size override (None = unbounded)
    # Integrator override; see SolverFLDOptions.odesolver in
    # solvers/solver_fld/options.py. None keeps the method chosen below.
    odesolver: Optional[Any] = None
    # Instants the caller wants the trajectory AT, increasing. When set, the ODE
    # is evaluated on exactly these (plus the endpoints) through the
    # integrator's own continuous extension, instead of on its step grid. Port
    # of MATLAB's options.tranpoints in solver_fluid_iteration.m; SolverENV sets
    # it to the sojourn quadrature grid, because reading that grid off a LINEAR
    # interpolation of the step grid is a first-order error the dense output
    # does not have. Unset, nothing changes for any other caller.
    tranpoints: Optional[np.ndarray] = None


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


def _fluid_theta(
    x: np.ndarray,
    SQ: np.ndarray,
    Sa: np.ndarray,
    pstar: Optional[np.ndarray] = None,
    isSourceState: Optional[np.ndarray] = None,
    isInfState: Optional[np.ndarray] = None,
    varclosure: bool = False,
) -> np.ndarray:
    """
    Mass in service, i.e. theta(x) of Ruuskanen et al., PEVA 151 (2021).

    Without smoothing, theta = x * min(S, sum_station(x)) / sum_station(x),
    eq. (12). With pstar set, theta = x * ghat with the p-norm ghat of
    eq. (26). An INF station has no min() to smooth, so its share stays 1;
    Sa holds the total population there, which would smooth it as a k = N
    queue.

    The same theta feeds the drift and the metrics: eq. (23) reads the
    utilization off the share the ODE integrated, so U and T must not
    revert to min() after a smoothed solve.

    Args:
        x: State vector (queue lengths)
        SQ: State-to-station queue mapping
        Sa: Server capacity per state
        pstar: P-norm smoothing parameters, None or empty for the hard min
        isSourceState: States belonging to Source stations
        isInfState: States belonging to INF (delay) stations
        varclosure: replace min(E[n], c) by E[min(n, c)] under the equilibrium
            geometric marginal of the station. Set ONLY by the degeneracy repair
            in solver_fld; see _fluid_fixed_point_is_degenerate for why.

    Returns:
        theta, the per-state mass in service
    """
    x = np.maximum(x, 0)  # Ensure non-negative

    # Compute total queue at each state's station
    sum_x_Qa = SQ @ x + 1e-8  # Add FineTol for numerical stability (matching MATLAB GlobalConstants.FineTol)

    if pstar is not None and len(pstar) > 0:
        # P-norm smoothed constraint as per Ruuskanen et al.
        ghat = np.ones_like(x)
        for i in range(len(x)):
            if isInfState is not None and isInfState[i]:
                continue
            x_val = sum_x_Qa[i]
            c_val = Sa[i]
            p_val = pstar[i] if i < len(pstar) else pstar[-1]

            if p_val > 0 and c_val > 0:
                ghat_val = 1.0 / np.power(1 + np.power(x_val / c_val, p_val), 1.0 / p_val)
                if np.isnan(ghat_val):
                    ghat[i] = 0.0
                else:
                    ghat[i] = ghat_val

        theta = x * ghat
    elif varclosure:
        # E[min(n, c)] UNDER A GEOMETRIC MARGINAL, not min(E[n], c).
        #
        #   n ~ Geometric(mean m)  =>  P(n >= k) = p^k with p = m/(1+m), and
        #   E[min(n,c)] = sum_{k=1..c} p^k = m * (1 - p^c).
        #
        # It has the two properties the hard min lacks and the repair needs:
        # it is STRICTLY INCREASING in m everywhere (slope 1/(1+m)^2 at c = 1,
        # so still 1e-2 at m = 9 -- a restoring force the integrator can follow
        # inside its horizon), and it carries the same asymptote, -> c as
        # m -> inf and -> m as m -> 0. It is the first-order face of what the
        # `dae` rung does by seeding the variance positive, which is why both
        # isolate the same fixed point.
        m = sum_x_Qa
        c = Sa.flatten()
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            p = np.where(m > 0, m / (1.0 + m), 0.0)
            min_vals = m * (1.0 - np.power(p, np.maximum(c, 0.0)))
        min_vals = np.nan_to_num(min_vals, nan=0.0, posinf=0.0, neginf=0.0)
        if isInfState is not None:
            # An INF station has a server per job: there is no min() to close,
            # and Sa holds the whole population there, which the closure would
            # otherwise read as a finite queue of that many servers.
            min_vals = np.where(isInfState, m, min_vals)
        min_vals = np.minimum(min_vals, sum_x_Qa)
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.where(sum_x_Qa > 1e-8, min_vals / sum_x_Qa, 1.0)
        ratio = np.nan_to_num(ratio, nan=1.0, posinf=1.0, neginf=0.0)
        theta = x * ratio
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
    return theta


def _fluid_ode(
    t: float,
    x: np.ndarray,
    W: np.ndarray,
    SQ: np.ndarray,
    Sa: np.ndarray,
    ALambda: np.ndarray,
    pstar: Optional[np.ndarray] = None,
    isSourceState: Optional[np.ndarray] = None,
    isInfState: Optional[np.ndarray] = None,
    varclosure: bool = False,
) -> np.ndarray:
    """
    Fluid ODE right-hand side for queueing network fluid analysis.

    dx/dt = W * theta(x) + ALambda

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
        isInfState: Boolean array indicating which states belong to INF stations

    Returns:
        dx/dt state derivative
    """
    theta = _fluid_theta(x, SQ, Sa, pstar, isSourceState, isInfState, varclosure)
    return W @ theta + ALambda.flatten()


def _fluid_station_groups(SQC, K, isSourceState=None, isInfState=None):
    """State index groups, one per (station, CLASS), keyed by class.

    PER CLASS, NOT PER STATION, and that is the whole correctness of the probe.
    A direction that moves a station's mass proportionally across ALL its
    classes is not a direction the model can take: a SelfLoopingClass is pinned
    at one station and can never leave it, so such a direction is infeasible,
    the drift is trivially unchanged along it, and the fixed point reads as
    degenerate. That is what it did to
    `sanity_CQN_2q_psfcfs_1class_1slcateachqueue`, whose two queues each hold a
    self-looping job: RespT came back 1.4336 against a baseline of 0.726303, a
    97% error, on a model with nothing wrong with it.

    Moving ONE class between two stations it actually occupies IS feasible, and
    a self-looping class occupies exactly one station, so no pair exists for it
    and no direction is proposed.

    SQC[i*K + r, a] is 1 exactly when state a belongs to station i and class r.
    Source and INF states are dropped: a Source carries theta = 0 by
    construction and an INF station carries theta = x with no min() to pin.

    Returns {class r: [state-index arrays, one per station that class occupies]},
    with only the classes that occupy at least two stations.
    """
    SQC = np.asarray(SQC)
    if SQC.ndim != 2 or K <= 0:
        return {}
    n = SQC.shape[1]
    M = SQC.shape[0] // K
    by_class = {}
    for r in range(K):
        per_station = []
        for i in range(M):
            members = [a for a in range(n)
                       if SQC[i * K + r, a] > 0
                       and not (isSourceState is not None and isSourceState[a])
                       and not (isInfState is not None and isInfState[a])]
            if members:
                per_station.append(np.asarray(members, dtype=int))
        if len(per_station) >= 2:
            by_class[r] = per_station
    return by_class


def _fluid_fixed_point_is_degenerate(x, rhs, SQC, K, isSourceState=None,
                                     isInfState=None, rate_scale=1.0):
    """Is the returned point one of a CONTINUUM of fixed points?

    A station whose queue exceeds its server count has theta pinned at the
    server count: min(S, sum_x) stops depending on sum_x, so the drift cannot
    tell one split of the mass between two such stations from another. Every
    split is then an equilibrium and the first-order method returns whichever
    one the integrator happened to stop at -- [9 1] where the exact answer is
    [5 5], on two identical saturated stations in a closed cycle.

    The test is direct rather than structural: move a little mass of ONE CLASS
    from one station to another along a population-conserving direction and see
    whether the drift moves at all. Both directions are tried, because the
    integrator typically stops on the BOUNDARY of the degenerate set, where one
    of the two does change the drift. The direction must be FEASIBLE -- see
    _fluid_station_groups for why moving a station's whole mass is not.

    The DIRECTIONAL DERIVATIVE is the scale-free quantity to threshold: a live
    direction moves the drift at the station's own service rate, a null one only
    by the FineTol the share carries, which leaves four orders between them.

    Returns False whenever the point is not a fixed point in the first place --
    a transient or timespan-limited run -- so no repair is ever attempted on a
    trajectory the caller asked to see mid-flight.
    """
    x = np.asarray(x, dtype=float).ravel()
    d0 = np.asarray(rhs(x), dtype=float).ravel()
    if d0.size == 0 or float(np.max(np.abs(d0))) > 1e-6 * max(1.0, float(np.max(np.abs(x)))):
        return False
    by_class = _fluid_station_groups(SQC, K, isSourceState, isInfState)
    if not by_class:
        return False
    eps = 1e-3 * max(1.0, float(np.max(np.abs(x))))
    rate_scale = max(float(rate_scale), 1e-12)
    for groups in by_class.values():
        mass = [float(np.sum(x[g])) for g in groups]
        for a in range(len(groups)):
            if mass[a] <= eps:
                continue
            for b in range(len(groups)):
                if a == b:
                    continue
                ga, gb = groups[a], groups[b]
                d = np.zeros_like(x)
                d[ga] -= x[ga] / mass[a]                   # take, proportionally
                if mass[b] > 0:
                    d[gb] += x[gb] / mass[b]               # give, proportionally
                else:
                    d[gb] += 1.0 / len(gb)
                dd = np.asarray(rhs(x + eps * d), dtype=float).ravel() - d0
                if float(np.max(np.abs(dd))) / eps <= 1e-4 * rate_scale:
                    return True
    return False


def _expand_eliminated(x_reduced, state_map, n_original):
    """Reduced ODE state back to the pre-elimination layout, eliminated slots at 0."""
    x = np.asarray(x_reduced, dtype=float).ravel()
    if state_map is None or x.size != len(state_map):
        return x_reduced
    full = np.zeros(int(n_original))
    full[np.asarray(state_map, dtype=int)] = x
    return full


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
    # Use stochastic complementation to eliminate non-station stateful nodes
    # (e.g., Router, Cache-as-ClassSwitch), matching MATLAB solver_fluid_matrix.m line 31:
    #   P = dtmc_stochcomp(P_full, station_indices);
    P = None
    if hasattr(sn, 'rt') and sn.rt is not None:
        from ...mc.dtmc import dtmc_stochcomp as _dtmc_stochcomp
        station_indices = []
        for ist in range(M):
            isf = int(sn.stationToStateful[ist]) if hasattr(sn, 'stationToStateful') and sn.stationToStateful is not None else ist
            for r in range(K):
                station_indices.append(isf * K + r)
        station_indices_arr = np.array(station_indices, dtype=int)
        if len(station_indices_arr) < sn.rt.shape[0]:
            # Non-station stateful nodes exist — use stochcomp to eliminate them
            P = _dtmc_stochcomp(sn.rt, station_indices_arr)
        else:
            # All stateful nodes are stations — simple extraction suffices
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
        # Build x0 from sn.state, the port of MATLAB SOLVER_FLUID_INITSOL.
        #
        # sn.state[isf] IS NOT A PER-CLASS COUNT VECTOR, and reading it as one is
        # a silent wrong answer rather than an error. Its layout depends on the
        # station's scheduling: PS/INF carry per-class counts, but FCFS, HOL,
        # LCFS and SIRO carry the BUFFER ORDERING -- one entry per job holding
        # that job's class -- so entry r is the class of the r-th queued job, not
        # the number of class-r jobs. Taking state_vec[r] there put ONE job in
        # the state of a station holding N of them, and the ODE conserves
        # whatever it is handed: a closed cycle referencing an FCFS station
        # returned a total population of 1 against N=6.
        #
        # State.toMarginal is the decoder for every layout, exactly as the MATLAB
        # reference calls it, and nir is the per-class count regardless of
        # scheduling.
        from ....api.state.marginal import toMarginal

        idx = 0
        for ist in range(M):
            node_idx_ist = int(sn.stationToNode[ist]) if ist < len(sn.stationToNode) else ist
            isf = int(sn.stationToStateful[ist]) if hasattr(sn, 'stationToStateful') and sn.stationToStateful is not None else ist

            nir_i = None
            kir_i = None
            if isf < len(sn.state) and sn.state[isf] is not None:
                _, nir_i, _, kir_i = toMarginal(sn, node_idx_ist, np.asarray(sn.state[isf]))
                nir_i = np.atleast_2d(np.asarray(nir_i, dtype=float))[0]
                kir_i = np.asarray(kir_i, dtype=float)
                if kir_i.ndim == 3:
                    kir_i = kir_i[0]
                else:
                    kir_i = np.atleast_2d(kir_i)

            for k in range(K):
                nphases_ik = int(phases[ist, k])
                if nphases_ik == 0:
                    continue

                is_source = (node_idx_ist < len(sn.nodetype) and sn.nodetype[node_idx_ist] == NodeType.SOURCE)
                sched_ist = sched_dict.get(ist) if sched_dict else None
                if sched_ist == SchedStrategy.EXT:
                    is_source = True

                if is_source:
                    # Source stations: no mass (arrivals via ALambda)
                    idx += nphases_ik
                    continue

                if nir_i is None or k >= len(nir_i):
                    idx += nphases_ik
                    continue

                # A job in the WAITING BUFFER has no phase yet, so it is restarted
                # in phase 1; only the jobs actually in service carry kir.
                in_service = np.zeros(nphases_ik)
                for ph in range(nphases_ik):
                    if k < kir_i.shape[0] and ph < kir_i.shape[1]:
                        in_service[ph] = kir_i[k, ph]

                total_k = float(nir_i[k])
                if np.isnan(total_k):
                    total_k = 0.0
                x0[idx] = total_k - float(np.sum(in_service[1:]))
                for ph in range(1, nphases_ik):
                    x0[idx + ph] = in_service[ph]
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
    # An INF station carries no min() to smooth, and Sa holds the population
    # there, which the p-norm would otherwise read as a k = N queue
    isInfState = np.zeros(W.shape[0], dtype=bool)
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
                if i < len(nservers_orig) and np.isinf(nservers_orig[i]):
                    isInfState[state_idx] = True
                state_idx += 1

    # STOCHASTIC-COMPLEMENT THE INSTANTANEOUS STATES OUT OF THE LINEAR GENERATOR,
    # the same reduction the event-set routes take, so this route does not
    # integrate an InfRate mode either. Port of MATLAB solver_fluid_matrix.m /
    # eliminate_immediate_matrix.m. Here W is the TRANSPOSE of a generator --
    # W[j,i] is the rate i -> j, because the drift is W @ theta -- so the block
    # splits by COLUMN and every per-state object is projected alongside it.
    #
    # The read-off matrices are CORRECTED, not truncated: an eliminated state
    # holds O(1/InfRate) mass and yet carries a FINITE throughput, because the
    # rate read off it is InfRate itself, so dropping its column would silently
    # delete every completion the instantaneous phase makes.
    try:
        from ....solvers.solver_fld.immediate import fluid_hide_immediate
        _do_hide = fluid_hide_immediate(sn, options)
    except Exception:  # noqa: BLE001 - a handler-level options object may not carry the flag
        _do_hide = False
    imm_state_map = None
    n_pre_elim = W.shape[0]
    if _do_hide and W.shape[0] > 1:
        imm_tol = GlobalConstants.Immediate * (1.0 - 0.01)
        cfg = getattr(options, 'config', None) or {}
        if isinstance(cfg, dict) and cfg.get('immediate_tol') is not None:
            imm_tol = float(cfg['immediate_tol'])
        # outgoing rates of state i live in COLUMN i
        imm = np.nonzero(np.abs(W).max(axis=0) >= imm_tol)[0]
        timed = np.setdiff1d(np.arange(W.shape[0]), imm)
        if imm.size > 0 and timed.size > 1:
            Q = W.T  # generator in the usual row-source convention
            QTT = Q[np.ix_(timed, timed)]
            QTI = Q[np.ix_(timed, imm)]
            QIT = Q[np.ix_(imm, timed)]
            QII = Q[np.ix_(imm, imm)]
            try:
                neg_inv = np.linalg.inv(-QII)
                absorb_ii = neg_inv @ QIT      # absorption distribution
                sojourn_ti = QTI @ neg_inv     # mass held per unit timed mass
                Q_red = QTT + QTI @ absorb_ii
            except np.linalg.LinAlgError:
                Q_red = None
            if Q_red is not None and np.all(np.isfinite(Q_red)):
                x0 = x0[timed] + absorb_ii.T @ x0[imm]
                ALambda = (np.asarray(ALambda).ravel()[timed]
                           + absorb_ii.T @ np.asarray(ALambda).ravel()[imm]).reshape(-1, 1)
                SQC = SQC[:, timed] + SQC[:, imm] @ sojourn_ti.T
                SUC = SUC[:, timed] + SUC[:, imm] @ sojourn_ti.T
                STC = STC[:, timed] + STC[:, imm] @ sojourn_ti.T
                Qa = Qa[:, timed]
                SQ = SQ[np.ix_(timed, timed)]
                Sa = Sa[timed]
                pstar = pstar[timed]
                isSourceState = isSourceState[timed]
                isInfState = isInfState[timed]
                W = Q_red.T
                imm_state_map = timed

    # Time span
    min_rate = np.abs(W[W != 0]).min() if np.any(W != 0) else 1.0
    T_end = min(options.timespan[1], abs(10 * options.iter_max / min_rate))
    T_start = options.timespan[0] if np.isfinite(options.timespan[0]) else 0.0

    # Set by the degeneracy repair below when the drift had to be re-integrated
    # with a variance-carrying saturation term; the metrics must then be read
    # off the SAME share the drift used.
    varclosure_used = False

    # Check if W is essentially zero (equilibrium case)
    W_norm = np.linalg.norm(W)
    if W_norm < 1e-10:
        # W is essentially zero - system is at equilibrium
        # Return initial state as the solution
        t_vec = np.array([T_start, T_end])
        x_vec = np.vstack([x0, x0])  # Constant solution
    else:
        # Solve ODE
        def _rhs(t, x):
            return _fluid_ode(t, x, W, SQ, Sa, ALambda, pstar, isSourceState, isInfState)
        try:
            if not options.stiff:
                method = 'RK45'
            elif W.size and np.abs(W).max() >= 0.99 * GlobalConstants.Immediate:
                # Immediate transitions (rate ~GlobalConstants.Immediate=1e8, e.g.
                # LQN entry/activity layers) make W extremely stiff. LSODA's Fortran
                # core takes O(1e8) tiny steps and its callback aborts on the
                # transient overflow (returning a NaN solution). Use the implicit
                # BDF solver (scipy analogue of MATLAB's ode15s), which steps over
                # the fast modes once they reach quasi-equilibrium. Non-immediate
                # models keep LSODA unchanged.
                method = 'BDF'
            else:
                method = 'LSODA'
            # An explicit integrator wins over all three branches above, and is
            # the only way the in-tree line_solver.lib.lsoda is reached
            if getattr(options, 'odesolver', None) is not None:
                method = options.odesolver

            t_eval = None
            if options.tranpoints is not None and len(options.tranpoints) > 0:
                pts = np.asarray(options.tranpoints, dtype=float).ravel()
                pts = pts[(pts > T_start) & (pts < T_end)]
                if pts.size:
                    t_eval = np.unique(np.concatenate(([T_start], pts, [T_end])))

            # SUCCESSIVE TIME WINDOWS, not one integration to the horizon.
            # This port used to compute T_end = 10*iter_max/min_rate up front and
            # hand solve_ivp the whole span in a single call. MATLAB
            # (solver_fluid_iteration.m) and the JAR
            # (ClosingAndStateDepMethodsAnalyzer) both advance in consecutive
            # warm-started windows T_k = 10*k/min_rate over [T_{k-1}, T_k], and
            # stop as soon as the iteration's own geometric tail says the fixed
            # point is within tolerance. The total model time is the same either
            # way; what the windows buy is the EARLY STOP, so a model that settles
            # in one window pays one window instead of the full horizon, and what
            # they cost is an integrator restart per window. Aligning here keeps
            # the three codebases on one algorithm rather than three.
            _mx = options.odemaxstep if (options.odemaxstep is not None
                                         and np.isfinite(options.odemaxstep)) else np.inf

            def _integrate(t0, t1, y0):
                pts = None
                if options.tranpoints is not None and len(options.tranpoints) > 0:
                    q = np.asarray(options.tranpoints, dtype=float).ravel()
                    q = q[(q > t0) & (q < t1)]
                    if q.size:
                        pts = np.unique(np.concatenate(([t0], q, [t1])))
                r = solve_ivp(_rhs, [t0, t1], y0, method=method,
                              rtol=options.tol, atol=options.tol,
                              dense_output=True, t_eval=pts, max_step=_mx)
                return r.t, r.y.T

            earlystop = True
            cfg = getattr(options, 'config', None)
            if isinstance(cfg, dict) and cfg.get('fluid_earlystop') is not None:
                earlystop = bool(cfg['fluid_earlystop'])
            elif cfg is not None and getattr(cfg, 'fluid_earlystop', None) is not None:
                earlystop = bool(cfg.fluid_earlystop)

            # The residual cannot be driven below the error the integrator carries.
            drift_tol = max(getattr(options, 'iter_tol', options.tol), options.tol)
            slowest_rate = min_rate
            min_horizon = 10.0 / slowest_rate
            drift_safety = 0.01   # headroom on the tail, since rho is estimated
            rho_hist = [np.nan, np.nan, np.nan]
            drift_below = 0       # consecutive windows satisfying the residual test
            moved_prev = np.inf

            T0 = T_start
            T = 0.0
            it_w = 0
            goon = True
            x_prev = np.asarray(x0, dtype=float)
            t_chunks, x_chunks = [], []
            while ((np.isfinite(options.timespan[1]) and T < options.timespan[1])
                   or (goon and it_w < options.iter_max)):
                it_w += 1
                T = min(options.timespan[1], abs(10.0 / min_rate)) if it_w == 1 \
                    else min(options.timespan[1], abs(10.0 * it_w / min_rate))
                t_it, x_it = _integrate(T0, T, x_prev)
                if t_it.size == 0:
                    break
                t_chunks.append(t_it)
                x_chunks.append(x_it)
                x_end = x_it[-1, :]
                denom = np.sum(x_prev)
                moved = (np.linalg.norm(x_end - x_prev, 1) / 2.0 / denom) if denom > 0 else 0.0
                T0 = T
                # MOVED MASS IS NOT THE TERMINATION TEST: it is one window's motion
                # and drops exactly the geometric tail r*rho/(1-rho) still to come.
                # rho is read off the iteration itself, and the drift F(x) -- zero AT
                # a fixed point -- is an independent second bound; BOTH must hold on
                # two consecutive windows, past the slowest relaxation time.
                if earlystop and goon and T >= min_horizon and it_w > 1:
                    rho_hist[it_w % len(rho_hist)] = moved / max(moved_prev, GlobalConstants.Zero)
                    finite = [v for v in rho_hist if np.isfinite(v)]
                    rho = max(finite) if finite else np.inf
                    xe = np.asarray(x_end, dtype=float)
                    drift_displ = (np.linalg.norm(_rhs(T, xe), 1) / 2.0
                                   / max(np.sum(xe), GlobalConstants.Zero) / slowest_rate)
                    if rho < 1:
                        tail = moved * rho / (1.0 - rho)
                        if tail < drift_safety * drift_tol and drift_displ < drift_tol:
                            drift_below += 1
                            if drift_below >= 2:
                                goon = False
                        else:
                            drift_below = 0
                    else:
                        drift_below = 0
                moved_prev = moved
                x_prev = x_end
                if T >= options.timespan[1]:
                    goon = False

            if t_chunks:
                t_vec = np.concatenate(t_chunks)
                x_vec = np.vstack(x_chunks)
            else:
                t_vec = np.array([T_start])
                x_vec = np.asarray(x0, dtype=float).reshape(1, -1)

        except (NameError, AttributeError, TypeError, ImportError):
            # A PROGRAMMING ERROR IS NOT AN INTEGRATION FAILURE, and turning one into
            # a NaN table hides it: `GlobalConstants` shadowed by a function-local
            # import made every method='matrix' model return NaN, and forty fluid
            # tests failed with no message naming the cause. The NaN fallback below
            # stays for what it is for -- a drift the integrator cannot follow.
            raise
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

    # DEGENERATE DRIFT: re-integrate with a closed saturation term, do not touch
    # the answer that came back. min(E[n], c) is FLAT above the server count, so
    # a network of saturated stations has a CONTINUUM of fixed points and this
    # method returns whichever one the integrator stopped at -- [9 1] against an
    # exact [5 5] on two identical saturated stations in a closed cycle, and
    # [8 2] on the same pair with two servers each.
    #
    # The repair is applied to the DRIFT, not to the point: the same trajectory
    # is integrated again with E[min(n, c)] in place of min(E[n], c), which is
    # strictly increasing and therefore isolates one fixed point. A selection
    # rule imposed after the fact would not be a solution of anything.
    #
    # WHY A CLOSURE AND NOT A SMOOTHED min: any smoothing sharp enough to stay
    # faithful to min away from the kink is numerically FLAT far from it -- the
    # Boltzmann softmin at alpha = 20 carries a restoring force of exp(-160) at
    # the [9 1] point, and the p-norm trades the two off directly (pstar = 2
    # recovers [5 5], pstar = 8 gives [7.64 2.36], pstar = 128 gives
    # [8.94 1.06]). The closure has no such trade-off because its slope comes
    # from the VARIANCE of the marginal rather than from a smoothing width.
    #
    # Only a model that is ACTUALLY degenerate pays for it: the test is a
    # null-direction probe at the returned point, so a well-posed model
    # integrates once and is bit-for-bit unchanged.
    if W_norm >= 1e-10 and (pstar is None or len(pstar) == 0):
        try:
            if _fluid_fixed_point_is_degenerate(
                    x_final, lambda xx: _rhs(T_end, xx), SQC, K,
                    isSourceState, isInfState,
                    float(np.max(np.abs(W))) if W.size else 1.0):
                def _rhs_closed(t, x):
                    return _fluid_ode(t, x, W, SQ, Sa, ALambda, pstar,
                                      isSourceState, isInfState, True)
                sol_c = solve_ivp(
                    _rhs_closed, [T_start, T_end], x0, method=method,
                    rtol=options.tol, atol=options.tol, dense_output=True,
                    t_eval=t_eval,
                    max_step=options.odemaxstep if (options.odemaxstep is not None
                                                    and np.isfinite(options.odemaxstep))
                    else np.inf,
                )
                x_c = np.maximum(sol_c.y.T[-1, :], 0)
                if np.all(np.isfinite(x_c)):
                    t_vec = sol_c.t
                    x_vec = sol_c.y.T
                    x_final = x_c
                    varclosure_used = True
                    if options.verbose:
                        from ...io.logging import line_printf
                        line_printf(
                            'Fluid: the first-order fixed point is not isolated '
                            '(two or more saturated stations), so the drift was '
                            're-integrated with a closed saturation term.\n')
        except Exception:
            # A failed repair must leave the unrepaired answer standing rather
            # than turning a wrong number into no number.
            pass

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

    # The same share the drift used, smoothed, closed or neither
    theta = _fluid_theta(x_final, SQ, Sa, pstar, isSourceState, isInfState,
                         varclosure_used)

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

    # Disabled (station,class) pairs (NaN service rate) hold no jobs and carry no
    # throughput; report 0 (matching MVA/MATLAB) instead of the NaN inherited from
    # the STC/SQC completion rates (np.full(nphases, sn.rates[i,r]) with a NaN
    # rate). Otherwise the NaN leaks into the LN activity throughput, which is
    # summed across layers (an activity is a disabled class in the layers that
    # only call it), turning a whole activity/entry row to NaN.
    if sn.rates is not None:
        rates_arr = np.asarray(sn.rates, dtype=float)
        rr = min(M, rates_arr.shape[0])
        cc = min(K, rates_arr.shape[1])
        disabled = np.isnan(rates_arr[:rr, :cc])
        Q[:rr, :cc][disabled] = 0.0
        U[:rr, :cc][disabled] = 0.0
        T[:rr, :cc][disabled] = 0.0

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
        theta_step = _fluid_theta(x_step, SQ, Sa, pstar, isSourceState, isInfState,
                                  varclosure_used)

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

    # Disabled (station,class) pairs (NaN service rate) hold no jobs and carry no
    # throughput at every time point, exactly as in the steady-state block above;
    # without this the NaN inherited from the STC/SQC completion rates leaks into
    # the transient traces (and, through the LN coupled transient, into the
    # inter-layer demand trajectories).
    if sn.rates is not None:
        rates_arr = np.asarray(sn.rates, dtype=float)
        for i in range(min(M, rates_arr.shape[0])):
            for r in range(min(K, rates_arr.shape[1])):
                if np.isnan(rates_arr[i, r]):
                    Qt[i][r][:] = 0.0
                    Ut[i][r][:] = 0.0
                    Tt[i][r][:] = 0.0

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
        # Handed back in the PRE-ELIMINATION layout: the caller stores it as
        # options.init_sol for the next iterate, which rebuilds the state over
        # every phase, so a vector shortened by the immediate elimination would
        # run that iterate off the end. The eliminated slots come back empty,
        # which is what they hold.
        odeStateVec=(_expand_eliminated(x_final, imm_state_map, n_pre_elim)
                     if imm_state_map is not None else x_final),
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
