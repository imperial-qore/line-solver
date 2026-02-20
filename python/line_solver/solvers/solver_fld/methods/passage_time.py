"""
Passage Time Distribution Analysis for SolverFLD.

Computes response time CDFs (passage time distributions) for stations in queueing
networks via network augmentation and ODE integration.

Algorithm: Network Augmentation (transient class approach)
1. Run steady-state analysis to get ODE state vector
2. Add transient job class (K+1) with initial jobs at target station
3. Route transient jobs to return to original classes upon completion at target station
4. Solve augmented network ODE system using the same matrix formulation as steady-state
5. Track transient fluid over time: F(t) = 1 - sum(transient_fluid) / initial

Features:
- Supports both open and closed networks
- Reuses handler's W matrix construction and ODE function for accurate dynamics
- Adaptive refinement for large CDF jumps
- Automatic time extension when CDF(0) > 0.01
- Compatible with all SolverFLD methods (matrix, softmin, etc.)

Reference: Solver_Fluid_Passage_Time.m from MATLAB LINE implementation
"""

import numpy as np
from scipy.integrate import solve_ivp
from typing import Optional, Tuple, Dict, Any, List
import warnings

from ..options import SolverFLDOptions
from ....api.sn import SchedStrategy


def compute_passage_time_cdf(
    sn,
    station_idx: int,
    job_class: int,
    options: SolverFLDOptions,
    steady_state_vec: Optional[np.ndarray] = None,
    t_span: Optional[Tuple[float, float]] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute response time CDF for a station using transient fluid analysis.

    Uses network augmentation approach with the same W matrix formulation
    as the steady-state solver (Ruuskanen et al., PEVA 151, 2021):
    1. Create augmented network with transient class K+1
    2. Build augmented proc/pie/rt structures
    3. Construct W matrix using W = psi + B * P * A' formulation
    4. Solve augmented ODE system
    5. Track: F(t) = 1 - remaining_transient_fluid / initial_transient_fluid

    Parameters
    ----------
    sn : NetworkStruct
        Network structure with service rates, routing, phases
    station_idx : int
        Station index for response time CDF
    job_class : int
        Job class index
    options : SolverFLDOptions
        Solver configuration
    steady_state_vec : np.ndarray, optional
        ODE state vector from steady-state solution
    t_span : tuple, optional
        Time interval for integration (t_min, t_max)

    Returns
    -------
    t : np.ndarray
        Time points
    cdf : np.ndarray
        CDF values F(t) = P(response_time <= t)
    """
    M = sn.nstations
    K = sn.nclasses

    # Get phases per station-class
    # If phases is None, compute from proc structure
    if sn.phases is not None:
        phases = np.asarray(sn.phases, dtype=int)
    else:
        phases = np.ones((M, K), dtype=int)
        # Try to extract phases from proc structure
        if hasattr(sn, 'proc') and sn.proc is not None:
            for i in range(min(M, len(sn.proc))):
                if sn.proc[i] is not None:
                    for r in range(min(K, len(sn.proc[i]))):
                        proc_ir = sn.proc[i][r]
                        if isinstance(proc_ir, dict):
                            # Dict-based format: {'k': n, 'mu': rate} for Erlang
                            if 'k' in proc_ir:
                                phases[i, r] = proc_ir['k']
                            elif 'nphases' in proc_ir:
                                phases[i, r] = proc_ir['nphases']
                            else:
                                phases[i, r] = 1  # Exponential
                        elif isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 1:
                            # Matrix-based format: [D0, D1] phase-type
                            D0 = np.asarray(proc_ir[0])
                            if D0.ndim == 2:
                                phases[i, r] = D0.shape[0]
                            else:
                                phases[i, r] = 1

    # Check if station/class has valid phases
    if phases[station_idx, job_class] == 0:
        # No service at this station-class, return trivial CDF
        t = np.array([0.0, 1.0])
        cdf = np.array([1.0, 1.0])
        return t, cdf

    # Get service rates
    rates = sn.rates if sn.rates is not None else np.ones((M, K))
    rates = np.asarray(rates)
    service_rate = rates[station_idx, job_class] if rates[station_idx, job_class] > 0 else 1.0

    # Get process structures from sn - convert dict-based to matrix-based if needed
    pie = sn.pie if hasattr(sn, 'pie') and sn.pie else {}

    # Build proper proc structure with phase-type matrices
    proc = {}
    if hasattr(sn, 'proc') and sn.proc is not None:
        for i in range(min(M, len(sn.proc))):
            proc[i] = {}
            if sn.proc[i] is not None:
                for r in range(min(K, len(sn.proc[i]))):
                    proc_ir = sn.proc[i][r]
                    if isinstance(proc_ir, dict):
                        # Convert dict-based to matrix-based
                        n_phases = phases[i, r]
                        if 'k' in proc_ir and 'mu' in proc_ir:
                            # Erlang(k, mu): D0 = -mu*I + mu*(superdiag), D1 = last row
                            k = proc_ir['k']
                            mu = proc_ir['mu']
                            D0 = np.zeros((k, k))
                            for p in range(k):
                                D0[p, p] = -mu
                                if p < k - 1:
                                    D0[p, p + 1] = mu
                            D1 = np.zeros((k, 1))
                            D1[k - 1, 0] = mu
                            proc[i][r] = [D0, D1]
                        elif 'rate' in proc_ir:
                            # Exponential
                            rate = proc_ir['rate']
                            proc[i][r] = [np.array([[-rate]]), np.array([[rate]])]
                        else:
                            # Default exponential fallback
                            rate = rates[i, r] if rates[i, r] > 0 else 1.0
                            proc[i][r] = [np.array([[-rate]]), np.array([[rate]])]
                    elif isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 2:
                        # Already matrix-based
                        proc[i][r] = proc_ir
                    else:
                        # Default exponential
                        rate = rates[i, r] if rates[i, r] > 0 else 1.0
                        proc[i][r] = [np.array([[-rate]]), np.array([[rate]])]
    rt = sn.rt if sn.rt is not None else None
    nservers = sn.nservers.flatten() if sn.nservers is not None else np.ones(M)

    # Compute slowest rate for time scaling
    slowrate = np.zeros((M, K))
    for i in range(M):
        for r in range(K):
            slowrate[i, r] = np.inf
            if rates[i, r] > 0:
                slowrate[i, r] = rates[i, r]

    nonzero_rates = slowrate.flatten()
    nonzero_rates = nonzero_rates[nonzero_rates > 0]
    nonzero_rates = nonzero_rates[~np.isinf(nonzero_rates)]
    if len(nonzero_rates) > 0:
        min_rate = np.min(nonzero_rates)
    else:
        min_rate = 1.0

    # Time span for integration
    if t_span is None:
        T = 100.0 / min_rate
        t_span = (0.0, T)
    T = t_span[1]

    # Get chain information
    # sn.chains can be in two formats:
    # 1) Index format (1D): chains[class_idx] = chain_index that class belongs to
    # 2) Membership format (2D): chains[chain_idx, class_idx] = 1 if class in chain
    chains_raw = sn.chains if sn.chains is not None else np.arange(K)
    chains_raw = np.asarray(chains_raw)

    if chains_raw.ndim == 1:
        # Index format: chains_raw[class_idx] = chain index
        # Convert to find classes in same chain as job_class
        job_chain = int(chains_raw[job_class]) if job_class < len(chains_raw) else 0
        classes_in_chain = np.where(chains_raw == job_chain)[0]
    else:
        # Membership format: chains_raw[chain_idx, class_idx] = 1 if member
        nchains = chains_raw.shape[0]
        chain_idx = 0
        for k in range(nchains):
            if job_class < chains_raw.shape[1] and chains_raw[k, job_class] == 1:
                chain_idx = k
                break
        classes_in_chain = np.where(chains_raw[chain_idx, :] == 1)[0]

    # Ensure job_class is in classes_in_chain (fallback if chain detection failed)
    if len(classes_in_chain) == 0 or job_class not in classes_in_chain:
        classes_in_chain = np.array([job_class])

    # Build augmented system (K+1 classes)
    Kc = K + 1
    phases_c = np.zeros((M, Kc), dtype=int)
    phases_c[:, :K] = phases
    phases_c[:, K] = phases[:, job_class]  # Transient class has same phases as job_class

    # Total number of phases in augmented system
    total_phases_c = int(np.sum(phases_c))

    # Build augmented proc structure
    # proc[i][r] = [psi_matrix, completion_matrix] for each station-class
    new_proc = {}
    new_pie = {}
    for i in range(M):
        new_proc[i] = {}
        new_pie[i] = {}
        for r in range(K):
            if proc and i in proc and r in proc[i]:
                new_proc[i][r] = proc[i][r]
            else:
                # Default exponential
                rate = rates[i, r] if rates[i, r] > 0 else 1.0
                new_proc[i][r] = [np.array([[-rate]]), np.array([[rate]])]
            if pie and i in pie and r in pie[i]:
                new_pie[i][r] = pie[i][r]
            else:
                new_pie[i][r] = np.array([1.0])
        # Transient class (K) copies from job_class
        new_proc[i][K] = new_proc[i][job_class]
        new_pie[i][K] = new_pie[i][job_class]

    # Build augmented routing matrix
    # new_rt is (M*Kc) x (M*Kc)
    new_rt = np.zeros((M * Kc, M * Kc))

    # Copy original routing among basic classes
    if rt is not None:
        for l in range(K):
            for m in range(K):
                for i in range(M):
                    for j in range(M):
                        src = i * K + l
                        dst = j * K + m
                        if src < rt.shape[0] and dst < rt.shape[1]:
                            new_src = i * Kc + l
                            new_dst = j * Kc + m
                            new_rt[new_src, new_dst] = rt[src, dst]

    # Copy routing for transient class (follows same routing as job_class)
    # Except at target station where it returns to original classes
    if rt is not None:
        for i in range(M):
            for j in range(M):
                if i != station_idx:
                    # Not at target station: transient stays transient
                    src = i * K + job_class
                    dst = j * K + job_class
                    if src < rt.shape[0] and dst < rt.shape[1]:
                        new_src = i * Kc + K
                        new_dst = j * Kc + K
                        new_rt[new_src, new_dst] = rt[src, dst]
                else:
                    # At target station: transient returns to original classes
                    for l in classes_in_chain:
                        src = i * K + job_class
                        dst = j * K + l
                        if src < rt.shape[0] and dst < rt.shape[1]:
                            new_src = i * Kc + K
                            new_dst = j * Kc + l
                            new_rt[new_src, new_dst] = rt[src, dst]
    else:
        # Default cyclic routing
        for i in range(M):
            j = (i + 1) % M
            for r in range(Kc):
                new_rt[i * Kc + r, j * Kc + r] = 1.0

    # Build W matrix using the same formulation as handler.py
    # W = psi + (B * P * A')^T
    W, q_indices = _build_augmented_W(M, Kc, phases_c, new_proc, new_pie, new_rt, rates)

    # Setup initial state
    if steady_state_vec is None:
        # Default: uniform distribution
        total_phases_orig = int(np.sum(phases))
        steady_state_vec = np.ones(total_phases_orig)
        if sn.nclosedjobs > 0:
            steady_state_vec = steady_state_vec * sn.nclosedjobs / total_phases_orig

    y0 = np.asarray(steady_state_vec).flatten()

    # Build y0_c for augmented system
    y0_c = np.zeros(total_phases_c)
    fluid_c = 0.0

    # Map from original phases to augmented phases
    idx_orig = 0
    for i in range(M):
        for r in range(K):
            n_phases_orig = int(phases[i, r])

            # Get augmented index for this station/class (in augmented system)
            idx_aug = int(sum(sum(phases_c[ii, rr] for rr in range(Kc)) for ii in range(i)) + sum(phases_c[i, :r]))

            # Get augmented index for transient class at station i
            idx_trans = int(sum(sum(phases_c[ii, rr] for rr in range(Kc)) for ii in range(i)) + sum(phases_c[i, :K]))

            for k in range(n_phases_orig):
                if idx_orig < len(y0):
                    if i == station_idx and r == job_class:
                        # Move fluid to transient class (all to phase 0)
                        if k == 0:
                            y0_c[idx_trans] = np.sum(y0[idx_orig:idx_orig + n_phases_orig])
                            fluid_c += np.sum(y0[idx_orig:idx_orig + n_phases_orig])
                    else:
                        # Keep in original class
                        if idx_aug + k < len(y0_c):
                            y0_c[idx_aug + k] = y0[idx_orig]
                    idx_orig += 1

    if fluid_c < 1e-14:
        # No fluid at this station-class
        t = np.array([0.0, 1.0])
        cdf = np.array([1.0, 1.0])
        return t, cdf

    # Get indices of transient class states at TARGET station only
    # (MATLAB code comment: "this used to be at all stations" - now only target station)
    transient_indices = []
    for i in [station_idx]:  # Only track transient at target station
        base_idx = int(sum(sum(phases_c[ii, rr] for rr in range(Kc)) for ii in range(i)))
        base_idx += int(sum(phases_c[i, :K]))
        n_trans = int(phases_c[i, K])
        for k in range(n_trans):
            transient_indices.append(base_idx + k)

    # Server capacities
    nservers_aug = np.zeros(M)
    for i in range(M):
        if np.isinf(nservers[i]):
            nservers_aug[i] = sn.nclosedjobs if sn.nclosedjobs > 0 else 1000
        else:
            nservers_aug[i] = nservers[i]

    # Sa (server capacities per phase)
    Sa = np.zeros(total_phases_c)
    for i in range(M):
        for r in range(Kc):
            n_phases_r = int(phases_c[i, r])
            for k in range(n_phases_r):
                idx = q_indices[i, r] + k
                Sa[idx] = nservers_aug[i]

    # Build phase-to-station mapping for computing station queue lengths
    phase_to_station = np.zeros(total_phases_c, dtype=int)
    idx = 0
    for i in range(M):
        for r in range(Kc):
            nphases = int(phases_c[i, r])
            for k in range(nphases):
                phase_to_station[idx] = i
                idx += 1

    # Compute station index ranges for each station
    station_idx_ranges = []
    for i in range(M):
        idx_start = int(sum(sum(phases_c[ii, rr] for rr in range(Kc)) for ii in range(i)))
        idx_end = idx_start + int(sum(phases_c[i, :]))
        station_idx_ranges.append((idx_start, idx_end))

    # Get scheduling strategy for each station
    station_sched = np.full(M, SchedStrategy.PS, dtype=int)  # Default to PS
    if hasattr(sn, 'sched') and sn.sched is not None:
        for i in range(M):
            if isinstance(sn.sched, dict) and i in sn.sched:
                station_sched[i] = sn.sched[i]
            elif isinstance(sn.sched, (list, np.ndarray)) and i < len(sn.sched):
                station_sched[i] = sn.sched[i]

    # Define ODE function using state-dependent method (same as MATLAB's ode_statedep)
    # Different scheduling strategies have different rate scaling:
    # - INF: No scaling (all jobs get full service rate, infinite servers)
    # - PS/FCFS: Scale by min(ni, ci) / ni when ni > ci
    def ode_rhs(t, x):
        x = np.maximum(x, 0)

        # Compute total queue at each station
        station_queues = np.zeros(M)
        for i in range(M):
            idx_start, idx_end = station_idx_ranges[i]
            station_queues[i] = np.sum(x[idx_start:idx_end])

        # Compute scaling factor for each phase based on station queue vs capacity
        # and scheduling strategy (matches MATLAB's ode_statedep behavior)
        ghat = np.zeros_like(x)
        for idx in range(len(x)):
            station_i = phase_to_station[idx]
            sched_i = station_sched[station_i]
            ni = station_queues[station_i]
            ci = nservers_aug[station_i]

            if sched_i == SchedStrategy.INF:
                # INF (delay/infinite server): all jobs served in parallel, no scaling
                ghat[idx] = 1.0
            elif sched_i in (SchedStrategy.PS, SchedStrategy.FCFS, SchedStrategy.DPS):
                # PS/FCFS/DPS: scale by c/n when queue exceeds capacity
                if ni > ci:
                    ghat[idx] = ci / ni
                else:
                    ghat[idx] = 1.0
            else:
                # Default behavior (other strategies): treat like PS
                if ni > ci:
                    ghat[idx] = ci / ni
                else:
                    ghat[idx] = 1.0

        return W @ (x * ghat)

    # Solve ODE
    tol = options.tol if options.tol else 1e-6
    iter_max = options.iter_max if options.iter_max else 100

    full_t = []
    full_y = []
    t_current = 0.0
    y_current = y0_c.copy()
    T_step = T / 10

    for iteration in range(iter_max):
        t_end = min(t_current + T_step, T * (1 + iteration * 0.5))

        try:
            method = 'LSODA' if options.stiff else 'RK45'
            odeopt = {'rtol': tol, 'atol': tol * 1e-3}
            sol = solve_ivp(
                ode_rhs,
                [t_current, t_end],
                y_current,
                method=method,
                **odeopt,
                dense_output=True,
                max_step=(t_end - t_current) / 10,
            )

            if len(full_t) == 0:
                full_t.extend(sol.t.tolist())
                full_y.extend(sol.y.T.tolist())
            else:
                full_t.extend(sol.t[1:].tolist())
                full_y.extend(sol.y.T[1:].tolist())

            y_current = sol.y[:, -1]
            t_current = sol.t[-1]

            # Check if transient fluid is depleted
            transient_fluid = np.sum(np.maximum(y_current[transient_indices], 0))
            if transient_fluid < 1e-10 * fluid_c:
                break

        except Exception as e:
            warnings.warn(f"ODE integration failed: {str(e)}")
            break

    if len(full_t) == 0:
        return _exponential_cdf_fallback(service_rate, t_span)

    t = np.array(full_t)
    y = np.array(full_y)

    # Compute CDF: F(t) = 1 - transient_fluid(t) / initial_fluid
    transient_over_time = np.sum(np.maximum(y[:, transient_indices], 0), axis=1)
    cdf = 1.0 - transient_over_time / fluid_c
    cdf = np.clip(cdf, 0.0, 1.0)

    # Ensure CDF is monotonically increasing
    for i in range(1, len(cdf)):
        cdf[i] = max(cdf[i], cdf[i-1])

    # Adaptive refinement for large CDF jumps
    t, cdf = _refine_cdf_jumps(t, cdf, max_jump=0.0005)

    # Extend time if CDF doesn't reach high enough
    if len(cdf) > 0 and cdf[-1] < 0.995:
        remaining = 1.0 - cdf[-1]
        if remaining > 0.001:
            t_extend = np.linspace(t[-1], t[-1] * 3, 50)
            # Exponential decay of remaining mass
            lambda_tail = service_rate
            cdf_extend = 1.0 - remaining * np.exp(-lambda_tail * (t_extend - t[-1]))
            t = np.concatenate([t, t_extend[1:]])
            cdf = np.concatenate([cdf, cdf_extend[1:]])

    return t, cdf


def _build_augmented_W(
    M: int,
    Kc: int,
    phases_c: np.ndarray,
    proc: Dict,
    pie: Dict,
    rt: np.ndarray,
    rates: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build transition rate matrix W for augmented system with phase-type support.

    For ODE dx/dt = W @ (x * g(x)), W[j,i] represents rate from state i to state j.

    Parameters
    ----------
    M : int
        Number of stations
    Kc : int
        Number of classes (K+1 with transient class)
    phases_c : np.ndarray
        (M x Kc) phases per station-class
    proc : dict
        Process structures: proc[i][r] = [D0, D1] phase-type matrices
    pie : dict
        Initial phase probabilities: pie[i][r] = probability_vector
    rt : np.ndarray
        Routing probability matrix (M*Kc x M*Kc)
    rates : np.ndarray
        Service rates (M x K) for original classes

    Returns
    -------
    W : np.ndarray
        Transition rate matrix (total_phases x total_phases)
    q_indices : np.ndarray
        (M x Kc) phase index for start of each station-class
    """
    # Compute q_indices: phase index for start of each station-class
    q_indices = np.zeros((M, Kc), dtype=int)
    cumsum = 0
    for i in range(M):
        for r in range(Kc):
            q_indices[i, r] = cumsum
            cumsum += int(phases_c[i, r])

    total_phases = cumsum
    W = np.zeros((total_phases, total_phases))

    # Step 1: Add internal phase transitions from D0 matrices
    for i in range(M):
        for r in range(Kc):
            nphases = int(phases_c[i, r])
            if nphases == 0:
                continue

            base_idx = q_indices[i, r]

            # Get D0 matrix
            if proc and i in proc and r in proc[i]:
                proc_ir = proc[i][r]
                if isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 2:
                    D0 = np.asarray(proc_ir[0])
                else:
                    D0 = np.asarray(proc_ir)
            else:
                # Default exponential
                rate = rates[i, min(r, rates.shape[1] - 1)] if rates is not None else 1.0
                D0 = np.array([[-rate]])

            # Add D0 to W with proper transpose for off-diagonals
            # D0[k_from, k_to] = rate FROM phase k_from TO phase k_to
            # W[j, i] = rate from state i to state j, so W[k_to, k_from] = D0[k_from, k_to]
            for k_from in range(min(nphases, D0.shape[0])):
                for k_to in range(min(nphases, D0.shape[1])):
                    if k_from == k_to:
                        # Diagonal: total exit rate (negative)
                        W[base_idx + k_from, base_idx + k_to] = D0[k_from, k_to]
                    else:
                        # Off-diagonal: transpose for W convention
                        W[base_idx + k_to, base_idx + k_from] = D0[k_from, k_to]

    # Step 2: Add routing transitions (completions going to next station-class)
    if rt is not None:
        for src_i in range(M):
            for src_r in range(Kc):
                src_nphases = int(phases_c[src_i, src_r])
                if src_nphases == 0:
                    continue

                src_base = q_indices[src_i, src_r]
                src_rt_idx = src_i * Kc + src_r

                # Get D1 (completion rates) for source
                if proc and src_i in proc and src_r in proc[src_i]:
                    proc_ir = proc[src_i][src_r]
                    if isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 2:
                        D1 = np.asarray(proc_ir[1])
                    else:
                        D0_src = np.asarray(proc_ir)
                        D1 = -np.sum(D0_src, axis=1, keepdims=True)
                else:
                    rate = rates[src_i, min(src_r, rates.shape[1] - 1)] if rates is not None else 1.0
                    D1 = np.array([[rate]])

                # Sum D1 rows to get completion rate per phase
                if D1.ndim == 1:
                    completion_rates = D1.flatten()
                else:
                    completion_rates = np.sum(D1, axis=1).flatten()

                # Route to each destination
                for dst_i in range(M):
                    for dst_r in range(Kc):
                        dst_rt_idx = dst_i * Kc + dst_r

                        if src_rt_idx >= rt.shape[0] or dst_rt_idx >= rt.shape[1]:
                            continue

                        p_route = rt[src_rt_idx, dst_rt_idx]
                        if p_route <= 0:
                            continue

                        dst_nphases = int(phases_c[dst_i, dst_r])
                        if dst_nphases == 0:
                            continue

                        dst_base = q_indices[dst_i, dst_r]

                        # Get initial phase probabilities for destination
                        if pie and dst_i in pie and dst_r in pie[dst_i]:
                            pie_dst = np.asarray(pie[dst_i][dst_r]).flatten()
                        else:
                            pie_dst = np.zeros(dst_nphases)
                            pie_dst[0] = 1.0

                        # Normalize pie if needed
                        if np.sum(pie_dst) > 0:
                            pie_dst = pie_dst / np.sum(pie_dst)
                        else:
                            pie_dst = np.zeros(dst_nphases)
                            pie_dst[0] = 1.0

                        # Add transition: completion from src phase -> dst initial phases
                        for k_src in range(min(src_nphases, len(completion_rates))):
                            rate_out = completion_rates[k_src] * p_route
                            if rate_out <= 0:
                                continue

                            for k_dst in range(dst_nphases):
                                # W[dst, src] = rate from src to dst
                                W[dst_base + k_dst, src_base + k_src] += rate_out * pie_dst[k_dst]

    return W, q_indices


def _exponential_cdf_fallback(
    service_rate: float,
    t_span: Tuple[float, float]
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fallback exponential CDF approximation.

    Used when augmented system construction fails.
    """
    t = np.linspace(t_span[0], t_span[1], 100)
    cdf = 1.0 - np.exp(-service_rate * t)
    return t, cdf


def _refine_cdf_jumps(
    t: np.ndarray,
    cdf: np.ndarray,
    max_jump: float = 0.001
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Adaptively refine CDF at large jumps.

    Detects jumps in CDF exceeding max_jump threshold and interpolates
    additional points to smoothly capture dynamics.
    """
    if len(t) < 2:
        return t, cdf

    # Detect jumps
    jumps = np.abs(np.diff(cdf))
    large_jump_indices = np.where(jumps > max_jump)[0]

    if len(large_jump_indices) == 0:
        return t, cdf

    # Refine around large jumps
    t_refined = [t[0]]
    cdf_refined = [cdf[0]]

    for i in range(len(t) - 1):
        t_refined.append(t[i])
        cdf_refined.append(cdf[i])

        if i in large_jump_indices:
            # Insert refined points (linear interpolation)
            n_refine = 20
            t_fine = np.linspace(t[i], t[i + 1], n_refine + 2)[1:-1]
            cdf_fine = np.interp(t_fine, [t[i], t[i + 1]], [cdf[i], cdf[i + 1]])

            t_refined.extend(t_fine)
            cdf_refined.extend(cdf_fine)

    t_refined.append(t[-1])
    cdf_refined.append(cdf[-1])

    return np.array(t_refined), np.array(cdf_refined)


class PassageTimeMethod:
    """Passage Time solver using network augmentation."""

    def __init__(
        self,
        sn,
        station_idx: int = 0,
        job_class: int = 0,
        options: Optional[SolverFLDOptions] = None,
        steady_state_vec: Optional[np.ndarray] = None
    ):
        """Initialize Passage Time analysis.

        Parameters
        ----------
        sn : NetworkStruct
            Network structure
        station_idx : int, optional
            Station index for response time analysis (default: 0)
        job_class : int, optional
            Job class index (default: 0)
        options : SolverFLDOptions, optional
            Solver configuration
        steady_state_vec : np.ndarray, optional
            ODE state vector from steady-state solution
        """
        self.sn = sn
        self.station_idx = station_idx
        self.job_class = job_class
        self.options = options or SolverFLDOptions()
        self.steady_state_vec = steady_state_vec

    def compute_cdf(
        self,
        t_span: Optional[Tuple[float, float]] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute response time CDF via network augmentation.

        Parameters
        ----------
        t_span : tuple, optional
            Time interval (t_min, t_max) for integration

        Returns
        -------
        t : np.ndarray
            Time points for CDF evaluation
        cdf : np.ndarray
            CDF values F(t) = P(response_time <= t)
        """
        return compute_passage_time_cdf(
            self.sn,
            self.station_idx,
            self.job_class,
            self.options,
            self.steady_state_vec,
            t_span
        )
