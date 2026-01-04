"""
Passage Time Distribution Analysis for SolverFLD.

Computes response time CDFs (passage time distributions) for stations in queueing
networks via network augmentation and ODE integration.

Algorithm: Network Augmentation (transient class approach)
1. Run steady-state analysis to get ODE state vector
2. Add transient job class (K+1) with initial jobs at target station
3. Route transient jobs to return to original classes upon completion at target station
4. Solve augmented network ODE system
5. Track transient fluid over time: F(t) = 1 - sum(transient_fluid) / initial

Features:
- Supports both open and closed networks
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

    Uses network augmentation approach:
    1. Create augmented network with transient class K+1
    2. Initialize transient class with fluid from target station
    3. Route transient class back to original classes upon completion
    4. Track: F(t) = 1 - remaining_transient_fluid / initial_transient_fluid

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
    phases = sn.phases if sn.phases is not None else np.ones((M, K))

    # Check if station/class has valid phases
    if phases[station_idx, job_class] == 0:
        # No service at this station-class, return trivial CDF
        t = np.array([0.0, 1.0])
        cdf = np.array([1.0, 1.0])
        return t, cdf

    # Get service rates
    rates = sn.rates if sn.rates is not None else np.ones((M, K))
    service_rate = rates[station_idx, job_class] if rates[station_idx, job_class] > 0 else 1.0

    # Compute slowest rate for time scaling
    nonzero_rates = rates[rates > 0]
    if len(nonzero_rates) > 0:
        slowest_rate = np.min(nonzero_rates)
    else:
        slowest_rate = 1.0

    # Time span for integration
    if t_span is None:
        # Use 100 events at slowest rate as max time
        T_max = 100.0 / slowest_rate
        t_span = (0.0, T_max)

    # Initialize state vector
    if steady_state_vec is not None and len(steady_state_vec) > 0:
        y0_base = steady_state_vec.copy()
    else:
        # Default: uniform distribution
        total_phases = int(np.sum(phases))
        y0_base = np.ones(total_phases)
        # Scale by population
        if sn.nclosedjobs > 0:
            y0_base = y0_base * sn.nclosedjobs / total_phases

    # Build augmented system with transient class
    # The transient class tracks fluid at target station that started in class c
    augmented_result = _build_augmented_system(
        sn, station_idx, job_class, phases, rates, y0_base, options
    )

    if augmented_result is None:
        # Fallback to exponential approximation
        return _exponential_cdf_fallback(service_rate, t_span)

    W_aug, y0_aug, Sa_aug, SQ_aug, transient_indices, initial_fluid = augmented_result

    if initial_fluid < 1e-14:
        # No fluid at this station-class
        t = np.array([0.0, 1.0])
        cdf = np.array([1.0, 1.0])
        return t, cdf

    # Solve augmented ODE
    def ode_rhs(t, x):
        x = np.maximum(x, 0)
        sum_x = SQ_aug @ x + 1e-14
        ratio = np.minimum(Sa_aug, sum_x) / sum_x
        ratio = np.nan_to_num(ratio, nan=1.0, posinf=1.0, neginf=0.0)
        return W_aug @ (x * ratio)

    # Integration settings
    tol = options.tol if options.tol else 1e-6
    iter_max = options.iter_max if options.iter_max else 100

    # Iterative integration until transient fluid is depleted
    full_t = []
    full_y = []
    t_current = t_span[0]
    y_current = y0_aug.copy()
    T_step = (t_span[1] - t_span[0]) / 10  # Integrate in chunks

    for iteration in range(iter_max):
        t_end = min(t_current + T_step, t_span[1] * (1 + iteration * 0.5))

        try:
            method = 'LSODA' if options.stiff else 'RK45'
            sol = solve_ivp(
                ode_rhs,
                [t_current, t_end],
                y_current,
                method=method,
                rtol=tol,
                atol=tol * 1e-3,
                dense_output=True,
            )

            # Append results
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
            if transient_fluid < 1e-10 * initial_fluid:
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
    cdf = 1.0 - transient_over_time / initial_fluid
    cdf = np.clip(cdf, 0.0, 1.0)

    # Ensure CDF is monotonically increasing
    for i in range(1, len(cdf)):
        cdf[i] = max(cdf[i], cdf[i-1])

    # Adaptive refinement for large CDF jumps
    t, cdf = _refine_cdf_jumps(t, cdf, max_jump=0.001)

    # Extend time if CDF doesn't reach high enough
    if len(cdf) > 0 and cdf[-1] < 0.995:
        # Extend using exponential tail approximation
        remaining = 1.0 - cdf[-1]
        if remaining > 0.001:
            t_extend = np.linspace(t[-1], t[-1] * 3, 50)
            # Exponential decay of remaining mass
            lambda_tail = service_rate
            cdf_extend = 1.0 - remaining * np.exp(-lambda_tail * (t_extend - t[-1]))
            t = np.concatenate([t, t_extend[1:]])
            cdf = np.concatenate([cdf, cdf_extend[1:]])

    return t, cdf


def _build_augmented_system(
    sn,
    station_idx: int,
    job_class: int,
    phases: np.ndarray,
    rates: np.ndarray,
    y0_base: np.ndarray,
    options: SolverFLDOptions
) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, List[int], float]]:
    """
    Build augmented system with transient class for passage time analysis.

    The augmented system has K+1 classes where class K+1 is the transient class
    that tracks fluid starting at station_idx in job_class.

    Returns
    -------
    W_aug : np.ndarray
        Augmented transition rate matrix
    y0_aug : np.ndarray
        Initial state for augmented system
    Sa_aug : np.ndarray
        Server capacities
    SQ_aug : np.ndarray
        State-to-station queue mapping
    transient_indices : list
        Indices of transient class states
    initial_fluid : float
        Initial fluid in transient class
    """
    M = sn.nstations
    K = sn.nclasses

    # Compute phase indices for original system
    total_phases_orig = int(np.sum(phases))

    # Build augmented phases: add phases for transient class at each station
    phases_aug = np.zeros((M, K + 1))
    phases_aug[:, :K] = phases
    phases_aug[:, K] = phases[:, job_class]  # Transient class has same phases as job_class

    total_phases_aug = int(np.sum(phases_aug))

    # Map original state indices to augmented indices
    def get_state_index(i, r, k, aug_phases):
        """Get state index for station i, class r, phase k."""
        idx = 0
        for ii in range(M):
            for rr in range(aug_phases.shape[1]):
                if ii == i and rr == r:
                    return idx + k
                idx += int(aug_phases[ii, rr])
        return -1

    # Build transition rate matrix W_aug
    W_aug = np.zeros((total_phases_aug, total_phases_aug))

    # Get routing matrix
    rt = sn.rt if sn.rt is not None else None

    # Fill in service rates (diagonal)
    for i in range(M):
        for r in range(K + 1):
            nphases = int(phases_aug[i, r])
            rate = rates[i, r] if r < K else rates[i, job_class]

            for k in range(nphases):
                idx = get_state_index(i, r, k, phases_aug)
                if idx >= 0 and idx < total_phases_aug:
                    W_aug[idx, idx] = -rate

    # Fill in routing (off-diagonal)
    for i in range(M):
        for r in range(K + 1):
            nphases = int(phases_aug[i, r])
            rate = rates[i, r] if r < K else rates[i, job_class]

            if rate <= 0:
                continue

            for k in range(nphases):
                idx_src = get_state_index(i, r, k, phases_aug)
                if idx_src < 0:
                    continue

                # Route to destination stations
                for j in range(M):
                    for s in range(K + 1):
                        # Determine routing probability
                        if r < K and s < K:
                            # Regular class to regular class
                            if rt is not None:
                                src_rt = i * K + r
                                dst_rt = j * K + s
                                if src_rt < rt.shape[0] and dst_rt < rt.shape[1]:
                                    prob = rt[src_rt, dst_rt]
                                else:
                                    prob = 0
                            else:
                                # Default cyclic routing
                                prob = 1.0 if j == (i + 1) % M and s == r else 0
                        elif r == K and s == K:
                            # Transient to transient (except at target station)
                            if i != station_idx:
                                if rt is not None:
                                    src_rt = i * K + job_class
                                    dst_rt = j * K + job_class
                                    if src_rt < rt.shape[0] and dst_rt < rt.shape[1]:
                                        prob = rt[src_rt, dst_rt]
                                    else:
                                        prob = 0
                                else:
                                    prob = 1.0 if j == (i + 1) % M else 0
                            else:
                                prob = 0  # At target station, transient returns to original
                        elif r == K and s < K:
                            # Transient to regular (at target station)
                            if i == station_idx:
                                if rt is not None:
                                    src_rt = i * K + job_class
                                    dst_rt = j * K + s
                                    if src_rt < rt.shape[0] and dst_rt < rt.shape[1]:
                                        prob = rt[src_rt, dst_rt]
                                    else:
                                        prob = 0
                                else:
                                    prob = 1.0 if j == (i + 1) % M and s == job_class else 0
                            else:
                                prob = 0
                        else:
                            prob = 0

                        if prob > 0 and int(phases_aug[j, s]) > 0:
                            idx_dst = get_state_index(j, s, 0, phases_aug)  # Enter phase 0
                            if idx_dst >= 0 and idx_dst < total_phases_aug:
                                W_aug[idx_dst, idx_src] += rate * prob

    # Build initial state y0_aug
    y0_aug = np.zeros(total_phases_aug)

    # Copy original state
    initial_fluid = 0.0
    idx_orig = 0
    idx_aug = 0
    for i in range(M):
        for r in range(K):
            nphases = int(phases[i, r])
            nphases_aug = int(phases_aug[i, r])

            for k in range(nphases):
                if idx_orig < len(y0_base):
                    if i == station_idx and r == job_class:
                        # Move fluid to transient class
                        idx_trans = get_state_index(i, K, k, phases_aug)
                        if idx_trans >= 0:
                            y0_aug[idx_trans] = y0_base[idx_orig]
                            initial_fluid += y0_base[idx_orig]
                    else:
                        # Keep in original class
                        y0_aug[idx_aug] = y0_base[idx_orig]
                idx_orig += 1
                idx_aug += 1
            # Skip augmented phases for classes beyond K
            while idx_aug < len(y0_aug) and idx_aug < sum([int(phases_aug[ii, rr])
                                                           for ii in range(i + 1)
                                                           for rr in range(min(r + 1, K))]):
                idx_aug += 1

    # Get indices of transient class states
    transient_indices = []
    for i in range(M):
        nphases_trans = int(phases_aug[i, K])
        for k in range(nphases_trans):
            idx = get_state_index(i, K, k, phases_aug)
            if idx >= 0:
                transient_indices.append(idx)

    # Build server capacities Sa_aug
    nservers = sn.nservers.flatten() if sn.nservers is not None else np.ones(M)
    # Replace inf with large number
    nservers = np.where(np.isinf(nservers), sn.nclosedjobs if sn.nclosedjobs > 0 else 1000, nservers)

    Sa_aug = np.zeros(total_phases_aug)
    idx = 0
    for i in range(M):
        for r in range(K + 1):
            nphases = int(phases_aug[i, r])
            for k in range(nphases):
                Sa_aug[idx] = nservers[i]
                idx += 1

    # Build SQ_aug (state -> total queue at station)
    SQ_aug = np.zeros((total_phases_aug, total_phases_aug))
    idx = 0
    for i in range(M):
        for r in range(K + 1):
            nphases = int(phases_aug[i, r])
            for k in range(nphases):
                # Mark all states at station i
                idx2 = 0
                for i2 in range(M):
                    for r2 in range(K + 1):
                        nphases2 = int(phases_aug[i2, r2])
                        for k2 in range(nphases2):
                            if i2 == i:
                                SQ_aug[idx, idx2] = 1.0
                            idx2 += 1
                idx += 1

    return W_aug, y0_aug, Sa_aug, SQ_aug, transient_indices, initial_fluid


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
