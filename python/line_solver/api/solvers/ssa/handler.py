"""
SSA Solver handler.

Native Python implementation of SSA (Stochastic Simulation Algorithm) solver
handler that analyzes queueing networks through discrete-event simulation.

The SSA solver uses Gillespie's algorithm to simulate sample paths and
estimate steady-state performance metrics.

Port from:

"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple, Any
import time

from ...sn import (
    NetworkStruct,
    SchedStrategy,
    NodeType,
    sn_is_open_model,
    sn_is_closed_model,
    sn_has_open_classes,
)


@dataclass
class SolverSSAOptions:
    """Options for SSA solver."""
    method: str = 'default'
    tol: float = 1e-6
    verbose: bool = False
    samples: int = 10000  # Number of simulation samples/events
    timespan: Tuple[float, float] = (0.0, float('inf'))  # Simulation time window
    seed: int = 0  # Random seed for reproducibility
    cutoff: int = 10  # Cutoff for open class populations
    confidence_level: float = 0.95  # Confidence level for CI


@dataclass
class SolverSSAReturn:
    """
    Result of SSA solver handler.

    Attributes:
        Q: Mean queue lengths (M x K)
        U: Utilizations (M x K)
        R: Response times (M x K)
        T: Throughputs (M x K)
        C: Cycle times (1 x K)
        X: System throughputs (1 x K)
        Q_ci: Queue length confidence intervals
        U_ci: Utilization confidence intervals
        R_ci: Response time confidence intervals
        T_ci: Throughput confidence intervals
        total_time: Total simulated time
        runtime: Runtime in seconds
        method: Method used
        samples: Number of samples collected
    """
    Q: Optional[np.ndarray] = None
    U: Optional[np.ndarray] = None
    R: Optional[np.ndarray] = None
    T: Optional[np.ndarray] = None
    C: Optional[np.ndarray] = None
    X: Optional[np.ndarray] = None
    Q_ci: Optional[np.ndarray] = None
    U_ci: Optional[np.ndarray] = None
    R_ci: Optional[np.ndarray] = None
    T_ci: Optional[np.ndarray] = None
    total_time: float = 0.0
    runtime: float = 0.0
    method: str = "default"
    samples: int = 0


def _init_state(
    sn: NetworkStruct,
    cutoff: int = 10
) -> np.ndarray:
    """
    Initialize the simulation state.

    For closed networks, places all jobs at reference stations.
    For open networks, starts with empty queues.

    Args:
        sn: Network structure
        cutoff: Population cutoff for open classes

    Returns:
        Initial state vector (M x K)
    """
    M = sn.nstations
    K = sn.nclasses
    state = np.zeros((M, K))

    if sn.njobs is not None:
        N = sn.njobs.flatten()
    else:
        N = np.zeros(K)

    for k in range(K):
        if np.isfinite(N[k]) and N[k] > 0:
            # Closed class: place jobs at reference station
            ref_stat = int(sn.refstat[k]) if hasattr(sn, 'refstat') and k < len(sn.refstat) else 0
            if ref_stat < M:
                state[ref_stat, k] = N[k]

    return state


def _get_enabled_transitions(
    sn: NetworkStruct,
    state: np.ndarray
) -> List[Tuple[int, int, int, int, float]]:
    """
    Find all enabled transitions in current state.

    Returns list of (source_station, source_class, dest_station, dest_class, rate)
    for each enabled transition.

    Args:
        sn: Network structure
        state: Current state (M x K)

    Returns:
        List of enabled transitions with their rates
    """
    M = sn.nstations
    K = sn.nclasses
    transitions = []

    # Get service rates
    if hasattr(sn, 'rates') and sn.rates is not None:
        rates = np.asarray(sn.rates)
    else:
        rates = np.ones((M, K))

    # Get routing probabilities
    if hasattr(sn, 'rt') and sn.rt is not None:
        P = np.asarray(sn.rt)
    else:
        # Default: uniform routing to next station
        P = np.zeros((M * K, M * K))
        for i in range(M):
            for k in range(K):
                next_i = (i + 1) % M
                P[i * K + k, next_i * K + k] = 1.0

    # Get number of servers
    if hasattr(sn, 'nservers') and sn.nservers is not None:
        nservers = np.asarray(sn.nservers).flatten()
    else:
        nservers = np.ones(M)

    # External arrivals (for open/mixed networks)
    # Use sn_has_open_classes instead of sn_is_open_model to properly handle
    # mixed networks (which have both open and closed classes)
    is_open = sn_has_open_classes(sn)
    if is_open:
        for k in range(K):
            # Check if class k has external arrivals
            if sn.njobs is not None:
                N = sn.njobs.flatten()
                if np.isinf(N[k]):
                    # Find arrival rate
                    for ist in range(M):
                        if hasattr(sn, 'sched') and sn.sched is not None:
                            sched = sn.sched.get(ist, SchedStrategy.FCFS)
                        else:
                            sched = SchedStrategy.FCFS

                        # Compare by name to handle enum type mismatches
                        sched_name = sched.name if hasattr(sched, 'name') else str(sched)
                        is_ext = sched_name == 'EXT' or (hasattr(sched, 'value') and sched.value == 11)

                        if is_ext:
                            arr_rate = rates[ist, k] if ist < rates.shape[0] and k < rates.shape[1] else 0
                            if arr_rate > 0:
                                # External arrival goes to first queue
                                for jst in range(M):
                                    src_idx = ist * K + k
                                    dst_idx = jst * K + k
                                    if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                                        p = P[src_idx, dst_idx]
                                        if p > 0:
                                            transitions.append((-1, k, jst, k, arr_rate * p))

    # Service completions
    for ist in range(M):
        total_at_station = np.sum(state[ist, :])
        if total_at_station <= 0:
            continue

        if hasattr(sn, 'sched') and sn.sched is not None:
            sched = sn.sched.get(ist, SchedStrategy.FCFS)
        else:
            sched = SchedStrategy.FCFS

        # Compare by name to handle enum type mismatches
        sched_name = sched.name if hasattr(sched, 'name') else str(sched)

        # Skip source stations
        if sched_name == 'EXT' or (hasattr(sched, 'value') and sched.value == 11):
            continue

        for k in range(K):
            n_ik = state[ist, k]
            if n_ik <= 0:
                continue

            mu = rates[ist, k] if ist < rates.shape[0] and k < rates.shape[1] else 1.0

            if sched_name == 'INF':
                # Infinite server: rate = n * mu
                rate = n_ik * mu
            else:
                # FCFS, PS, etc.: rate = min(n, c) * mu * fraction
                c = nservers[ist] if ist < len(nservers) else 1
                if np.isinf(c):
                    c = total_at_station
                active_servers = min(total_at_station, c)
                rate = active_servers * mu * (n_ik / total_at_station)

            if rate <= 0:
                continue

            # Find destinations
            src_idx = ist * K + k
            total_routing_prob = 0.0

            for jst in range(M):
                for r in range(K):
                    dst_idx = jst * K + r

                    if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                        p = P[src_idx, dst_idx]
                    else:
                        p = 0

                    if p > 0:
                        transitions.append((ist, k, jst, r, rate * p))
                        total_routing_prob += p

            # For open networks: if routing probs sum to < 1, remaining goes to sink
            # Use destination -2 to indicate departure from system
            if total_routing_prob < 1.0 - 1e-6:
                sink_prob = 1.0 - total_routing_prob
                transitions.append((ist, k, -2, k, rate * sink_prob))

    return transitions


def _fire_transition(
    state: np.ndarray,
    transition: Tuple[int, int, int, int, float]
) -> np.ndarray:
    """
    Fire a transition and return new state.

    Args:
        state: Current state
        transition: (src_station, src_class, dst_station, dst_class, rate)
            src_st = -1 means external arrival
            dst_st = -2 means departure to sink

    Returns:
        New state after transition
    """
    new_state = state.copy()
    src_st, src_k, dst_st, dst_k, _ = transition

    if src_st >= 0:
        new_state[src_st, src_k] -= 1

    # dst_st = -2 means departure to sink (job leaves system)
    if dst_st >= 0:
        new_state[dst_st, dst_k] += 1

    return new_state


def solver_ssa_basic(
    sn: NetworkStruct,
    options: Optional[SolverSSAOptions] = None
) -> SolverSSAReturn:
    """
    Basic SSA solver using Gillespie algorithm.

    Simulates the queueing network and collects statistics for
    performance metric estimation.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverSSAReturn with performance metrics
    """
    start_time = time.time()

    if options is None:
        options = SolverSSAOptions()

    # Set random seed
    if options.seed > 0:
        np.random.seed(options.seed)

    M = sn.nstations
    K = sn.nclasses

    # Initialize state
    state = _init_state(sn, options.cutoff)

    # Statistics accumulators
    total_time = 0.0
    time_weighted_queue = np.zeros((M, K))
    time_weighted_util = np.zeros((M, K))
    departure_counts = np.zeros((M, K))
    arrival_counts = np.zeros((M, K))

    # Get number of servers for utilization calculation
    if hasattr(sn, 'nservers') and sn.nservers is not None:
        nservers = np.asarray(sn.nservers).flatten()
    else:
        nservers = np.ones(M)

    samples_collected = 0
    max_time = options.timespan[1] if np.isfinite(options.timespan[1]) else 1e6

    while samples_collected < options.samples and total_time < max_time:
        # Find enabled transitions
        transitions = _get_enabled_transitions(sn, state)

        if not transitions:
            # Deadlock or absorbing state
            if options.verbose:
                print(f"SSA: No enabled transitions at sample {samples_collected}")
            break

        # Total rate
        total_rate = sum(t[4] for t in transitions)

        if total_rate <= 0:
            break

        # Sample time to next event (exponential)
        dt = np.random.exponential(1.0 / total_rate)

        # Update time-weighted statistics before transition
        time_weighted_queue += dt * state

        # Utilization: fraction of servers busy
        for ist in range(M):
            total_at_station = np.sum(state[ist, :])
            c = nservers[ist] if ist < len(nservers) else 1
            if np.isinf(c):
                c = total_at_station if total_at_station > 0 else 1

            if hasattr(sn, 'sched') and sn.sched is not None:
                sched = sn.sched.get(ist, SchedStrategy.FCFS)
            else:
                sched = SchedStrategy.FCFS

            # Compare by name to handle enum type mismatches
            sched_name = sched.name if hasattr(sched, 'name') else str(sched)

            if sched_name == 'INF':
                # Delay station: utilization = jobs in service
                for k in range(K):
                    time_weighted_util[ist, k] += dt * state[ist, k]
            else:
                # Queue: fraction of servers busy
                busy_servers = min(total_at_station, c)
                for k in range(K):
                    if total_at_station > 0:
                        time_weighted_util[ist, k] += dt * busy_servers * (state[ist, k] / total_at_station) / c

        total_time += dt

        # Select which transition fires
        rand = np.random.random() * total_rate
        cumsum = 0.0
        selected = None
        for t in transitions:
            cumsum += t[4]
            if cumsum >= rand:
                selected = t
                break

        if selected is None:
            selected = transitions[-1]

        # Record departure and arrival
        src_st, src_k, dst_st, dst_k, _ = selected
        if src_st >= 0:
            departure_counts[src_st, src_k] += 1
        if dst_st >= 0:
            arrival_counts[dst_st, dst_k] += 1

        # Fire transition
        state = _fire_transition(state, selected)

        samples_collected += 1

    # Compute average metrics
    if total_time > 0:
        QN = time_weighted_queue / total_time
        UN = time_weighted_util / total_time
        TN = departure_counts / total_time
    else:
        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        TN = np.zeros((M, K))

    # Response times via Little's law: R = Q / T
    RN = np.zeros((M, K))
    for ist in range(M):
        for k in range(K):
            if TN[ist, k] > 0:
                RN[ist, k] = QN[ist, k] / TN[ist, k]

    # Cycle times and system throughput
    CN = np.sum(RN, axis=0).reshape(1, -1)
    XN = np.zeros((1, K))
    for k in range(K):
        ref_stat = int(sn.refstat[k]) if hasattr(sn, 'refstat') and k < len(sn.refstat) else 0
        if ref_stat < M:
            XN[0, k] = TN[ref_stat, k]

    # Clean up NaN values
    QN = np.nan_to_num(QN, nan=0.0)
    UN = np.nan_to_num(UN, nan=0.0)
    RN = np.nan_to_num(RN, nan=0.0)
    TN = np.nan_to_num(TN, nan=0.0)
    CN = np.nan_to_num(CN, nan=0.0)
    XN = np.nan_to_num(XN, nan=0.0)

    result = SolverSSAReturn()
    result.Q = QN
    result.U = UN
    result.R = RN
    result.T = TN
    result.C = CN
    result.X = XN
    result.total_time = total_time
    result.runtime = time.time() - start_time
    result.method = "serial"
    result.samples = samples_collected

    return result


def _run_ssa_replica(args: Tuple) -> Dict:
    """
    Run a single SSA simulation replica for parallel execution.

    This function is designed to be picklable for multiprocessing.

    Args:
        args: Tuple of (sn, samples, seed, cutoff, timespan)

    Returns:
        Dict with raw metrics from this replica
    """
    sn, samples, seed, cutoff, timespan = args

    # Set seed for this replica
    np.random.seed(seed)

    M = sn.nstations
    K = sn.nclasses

    # Initialize state
    state = _init_state(sn, cutoff)

    # Statistics accumulators
    Q_accum = np.zeros((M, K))
    service_completions = np.zeros((M, K))
    arrivals = np.zeros((M, K))
    busy_time = np.zeros((M, K))
    response_time_sum = np.zeros((M, K))
    response_count = np.zeros((M, K))

    current_time = timespan[0] if timespan else 0.0
    max_time = timespan[1] if timespan and len(timespan) > 1 else float('inf')
    last_event_time = current_time

    # Get rates
    rates = _compute_rates(sn, state)
    total_rate = np.sum(rates)

    for _ in range(samples):
        if total_rate <= 0:
            break

        if current_time >= max_time:
            break

        # Time to next event
        dt = np.random.exponential(1.0 / total_rate) if total_rate > 0 else float('inf')

        if current_time + dt > max_time:
            dt = max_time - current_time

        # Accumulate time-weighted statistics
        Q_accum += state * dt
        for i in range(M):
            for r in range(K):
                if state[i, r] > 0:
                    busy_time[i, r] += dt

        current_time += dt

        if current_time >= max_time:
            break

        # Select and execute event
        flat_rates = rates.flatten()
        probs = flat_rates / total_rate
        event_idx = np.random.choice(len(probs), p=probs)
        station = event_idx // K
        job_class = event_idx % K

        # Departure from station
        if state[station, job_class] > 0:
            state[station, job_class] -= 1
            service_completions[station, job_class] += 1

            # Route to next station
            next_station = _get_next_station(sn, station, job_class)
            if next_station >= 0 and next_station < M:
                state[next_station, job_class] += 1
                arrivals[next_station, job_class] += 1

        # Update rates
        rates = _compute_rates(sn, state)
        total_rate = np.sum(rates)

        last_event_time = current_time

    # Compute metrics from this replica
    total_time = last_event_time - (timespan[0] if timespan else 0.0)
    if total_time <= 0:
        total_time = 1.0

    Q_mean = Q_accum / total_time
    U_mean = busy_time / total_time
    T_mean = service_completions / total_time

    # Response times via Little's Law: R = Q / T
    R_mean = np.zeros((M, K))
    for i in range(M):
        for r in range(K):
            if T_mean[i, r] > 0:
                R_mean[i, r] = Q_mean[i, r] / T_mean[i, r]

    return {
        'Q': Q_mean,
        'U': U_mean,
        'R': R_mean,
        'T': T_mean,
        'total_time': total_time,
        'completions': service_completions.sum()
    }


def solver_ssa_parallel(
    sn: NetworkStruct,
    options: Optional[SolverSSAOptions] = None,
    n_workers: int = None,
    n_replicas: int = 4
) -> SolverSSAReturn:
    """
    Parallel SSA solver using multiprocessing.

    Runs multiple independent simulation replicas in parallel and
    aggregates results with confidence intervals.

    Args:
        sn: Network structure
        options: Solver options
        n_workers: Number of worker processes (default: CPU count)
        n_replicas: Number of independent replicas (default: 4)

    Returns:
        SolverSSAReturn with aggregated performance metrics and CIs
    """
    import multiprocessing as mp
    from scipy import stats

    start_time = time.time()

    if options is None:
        options = SolverSSAOptions()

    if n_workers is None:
        n_workers = min(mp.cpu_count(), n_replicas)

    # Prepare arguments for each replica
    base_seed = options.seed if options.seed > 0 else int(time.time())
    samples_per_replica = max(1000, options.samples // n_replicas)

    replica_args = []
    for i in range(n_replicas):
        seed = base_seed + i * 12345
        replica_args.append((sn, samples_per_replica, seed, options.cutoff, options.timespan))

    # Run replicas in parallel
    if options.verbose:
        print(f"Running {n_replicas} SSA replicas with {n_workers} workers...")

    try:
        with mp.Pool(processes=n_workers) as pool:
            results = pool.map(_run_ssa_replica, replica_args)
    except Exception as e:
        # Fallback to serial if multiprocessing fails
        if options.verbose:
            print(f"Parallel execution failed ({e}), falling back to serial.")
        results = [_run_ssa_replica(args) for args in replica_args]

    # Aggregate results
    M = sn.nstations
    K = sn.nclasses

    Q_all = np.array([r['Q'] for r in results])
    U_all = np.array([r['U'] for r in results])
    R_all = np.array([r['R'] for r in results])
    T_all = np.array([r['T'] for r in results])

    # Mean across replicas
    Q_mean = np.mean(Q_all, axis=0)
    U_mean = np.mean(U_all, axis=0)
    R_mean = np.mean(R_all, axis=0)
    T_mean = np.mean(T_all, axis=0)

    # Confidence intervals
    alpha = 1.0 - options.confidence_level
    t_critical = stats.t.ppf(1 - alpha / 2, n_replicas - 1) if n_replicas > 1 else 1.96

    Q_std = np.std(Q_all, axis=0, ddof=1) if n_replicas > 1 else np.zeros((M, K))
    U_std = np.std(U_all, axis=0, ddof=1) if n_replicas > 1 else np.zeros((M, K))
    R_std = np.std(R_all, axis=0, ddof=1) if n_replicas > 1 else np.zeros((M, K))
    T_std = np.std(T_all, axis=0, ddof=1) if n_replicas > 1 else np.zeros((M, K))

    se_factor = t_critical / np.sqrt(n_replicas)

    Q_ci = Q_std * se_factor
    U_ci = U_std * se_factor
    R_ci = R_std * se_factor
    T_ci = T_std * se_factor

    # System throughput
    X = np.sum(T_mean, axis=0, keepdims=True)

    # Cycle times (for closed networks)
    C = np.zeros((1, K))
    if sn.njobs is not None:
        N = sn.njobs.flatten()
        for r in range(K):
            if X[0, r] > 0 and np.isfinite(N[r]) and N[r] > 0:
                C[0, r] = N[r] / X[0, r]

    total_time = sum(r['total_time'] for r in results) / n_replicas
    runtime = time.time() - start_time

    if options.verbose:
        print(f"Parallel SSA completed in {runtime:.3f}s ({n_replicas} replicas)")

    return SolverSSAReturn(
        Q=Q_mean,
        U=U_mean,
        R=R_mean,
        T=T_mean,
        C=C,
        X=X,
        Q_ci=Q_ci,
        U_ci=U_ci,
        R_ci=R_ci,
        T_ci=T_ci,
        total_time=total_time,
        runtime=runtime,
        method='parallel',
        samples=options.samples
    )


def solver_ssa(
    sn: NetworkStruct,
    options: Optional[SolverSSAOptions] = None
) -> SolverSSAReturn:
    """
    Main SSA solver handler.

    Routes to appropriate simulation method based on options.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverSSAReturn with performance metrics
    """
    if options is None:
        options = SolverSSAOptions()

    method = options.method.lower()

    if method in ['default', 'serial', 'ssa']:
        return solver_ssa_basic(sn, options)
    elif method in ['parallel', 'para', 'ssa.parallel']:
        return solver_ssa_parallel(sn, options)
    elif method == 'nrm':
        # Next Reaction Method - same algorithm, different selection
        return solver_ssa_basic(sn, options)
    else:
        # Unknown method
        if options.verbose:
            print(f"Warning: Unknown SSA method '{method}'. Using serial.")
        return solver_ssa_basic(sn, options)


__all__ = [
    'solver_ssa',
    'solver_ssa_basic',
    'solver_ssa_parallel',
    'SolverSSAReturn',
    'SolverSSAOptions',
]
