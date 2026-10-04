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
- Automatic horizon extension while the tail misses the law, CDF(end) < 0.99
- Compatible with all SolverFLD methods (matrix, softmin, etc.)

Reference: Solver_Fluid_Passage_Time.m from MATLAB LINE implementation
"""

import numpy as np
from scipy.integrate import solve_ivp
from typing import Optional, Tuple, Dict, Any, List
import warnings

from ..options import SolverFLDOptions
from ..utils.closures import capacity_closure
from ....api.sn import SchedStrategy


def compute_passage_time_cdf(
    sn,
    station_idx: int,
    job_class: int,
    options: SolverFLDOptions,
    steady_state_vec: Optional[np.ndarray] = None,
    t_span: Optional[Tuple[float, float]] = None,
    sigma2: Optional[np.ndarray] = None
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
    sigma2 : np.ndarray, optional
        Per-station population variance the mean solve closed its drift at, i.e.
        the moment-closure `sigma2Drift`. The passage time is a SECOND solve on
        that fixed point, so it has to be driven by the same drift: closing
        min(n_i,c_i) at zero variance drains a station the mean solve holds below
        capacity at full rate, and the distribution then contradicts the mean the
        same solver reports. None or all-zero gives the first-order min().

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

    # proc AND pie must come from the same derivation: the native struct leaves
    # sn.pie None, and a pie shorter than its phase count silently truncates the
    # phase structure (see _kb/06-solver-catalog.md, Fluid passage time)
    from ..utils.phase_type import prepare_phase_type_structures
    proc, pie, _ = prepare_phase_type_structures(sn)
    rt = sn.rt if sn.rt is not None else None
    nservers = sn.nservers.flatten() if sn.nservers is not None else np.ones(M)

    # see _kb/06-solver-catalog.md (Fluid: "Python passage-time augmented
    # system") for why a Disabled service process must be zeroed here
    for i in range(M):
        for r in range(K):
            proc_ir = proc.get(i, {}).get(r) if isinstance(proc, dict) else None
            is_disabled = False
            if proc_ir is not None and isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 1:
                for mat in proc_ir:
                    arr = np.asarray(mat, dtype=float)
                    if arr.size and (np.any(np.isnan(arr)) or np.any(np.isinf(arr))):
                        is_disabled = True
                        break
            if is_disabled:
                phases[i, r] = 0

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

    # sn.chains: index format (1D, chains[class]=chain_idx) or membership
    # format (2D, chains[chain,class]=1 if member); see _kb/04-networkstruct.md
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

    # see _kb/06-solver-catalog.md (Fluid: "Python passage-time augmented
    # system") -- EXT mass conservation excludes the transient tracer class K
    ext_phase_indices = set()  # All phase indices belonging to EXT stations
    ext_class_ranges = {}  # (station, class) -> (first_phase_idx, [other_phase_indices])
    for i in range(M):
        if station_sched[i] == SchedStrategy.EXT:
            idx_start, idx_end = station_idx_ranges[i]
            for idx in range(idx_start, idx_end):
                ext_phase_indices.add(idx)
            for r in range(K):  # Only original classes, not transient class K
                n_ph = int(phases_c[i, r])
                if n_ph > 0:
                    first_idx = int(q_indices[i, r])
                    other_indices = list(range(first_idx + 1, first_idx + n_ph))
                    ext_class_ranges[(i, r)] = (first_idx, other_indices)

    # The closure the mean solve finished at, per station. Only the per-STATION
    # variance transfers to the augmented model: the transient class only
    # relabels a station population, whereas the coordinate covariance is indexed
    # by a state layout that class changes, so the class share stays the plug-in
    # ratio (as in SolverFluid.passageTimeOptions and getCdfRespT.m).
    sigma2_drift = np.zeros(M)
    if sigma2 is not None:
        s2arr = np.asarray(sigma2, dtype=float).ravel()
        sigma2_drift[:min(M, s2arr.size)] = s2arr[:min(M, s2arr.size)]

    # Define ODE function using closing method (matches MATLAB's ode_rates_closing)
    # - EXT: Mass conservation (total mass per class = 1, source doesn't evolve)
    # - INF: No scaling (all jobs get full service rate, infinite servers)
    # - PS/FCFS: Scale by min(ni, ci) / ni when ni > ci
    def ode_rhs(t, x):
        # x is NOT projected here: MATLAB leaves nonnegativity to
        # odeset('NonNegative'), which acts on the accepted step, so projecting the
        # drift's argument would make it discontinuous at x=0 (see closing.py)
        # Build rates vector (same as MATLAB: rates = x, then modify per strategy)
        rates = np.array(x, dtype=float, copy=True)

        # Compute total queue at each station
        station_queues = np.zeros(M)
        for i in range(M):
            idx_start, idx_end = station_idx_ranges[i]
            station_queues[i] = np.sum(x[idx_start:idx_end])

        for i in range(M):
            sched_i = station_sched[i]
            if sched_i == SchedStrategy.EXT:
                # EXT (Source): enforce mass conservation per class
                # First phase rate = 1 - sum(other phases), keeping total mass = 1
                for r in range(Kc):
                    key = (i, r)
                    if key in ext_class_ranges:
                        first_idx, other_indices = ext_class_ranges[key]
                        if other_indices:
                            rates[first_idx] = 1.0 - np.sum(x[other_indices])
                        else:
                            rates[first_idx] = 1.0
            elif sched_i == SchedStrategy.INF:
                # INF: rates = x (no scaling needed, already set)
                pass
            elif sched_i in (SchedStrategy.PS, SchedStrategy.FCFS, SchedStrategy.DPS):
                # PS/FCFS/DPS: scale by psi(n)/n, psi = min(n,c) closed at s2
                ni = station_queues[i]
                ci = nservers_aug[i]
                s2i = sigma2_drift[i]
                if s2i > 0:
                    if ni > 0:
                        h = capacity_closure(ni, ci, s2i)[0]
                        idx_start, idx_end = station_idx_ranges[i]
                        rates[idx_start:idx_end] = x[idx_start:idx_end] / ni * h
                elif ni > ci:
                    idx_start, idx_end = station_idx_ranges[i]
                    rates[idx_start:idx_end] = x[idx_start:idx_end] / ni * ci
            else:
                # Default: treat like PS
                ni = station_queues[i]
                ci = nservers_aug[i]
                s2i = sigma2_drift[i]
                if s2i > 0:
                    if ni > 0:
                        h = capacity_closure(ni, ci, s2i)[0]
                        idx_start, idx_end = station_idx_ranges[i]
                        rates[idx_start:idx_end] = x[idx_start:idx_end] / ni * h
                elif ni > ci:
                    idx_start, idx_end = station_idx_ranges[i]
                    rates[idx_start:idx_end] = x[idx_start:idx_end] / ni * ci

        dx = W @ rates
        # see _kb/06-solver-catalog.md (Fluid: "Python passage-time augmented
        # system") -- freeze EXT phases, mass conservation stays in rates above
        for idx in ext_phase_indices:
            dx[idx] = 0.0
        return dx

    # Solve ODE - match MATLAB solver_fluid_passage_time.m integration strategy
    tol = options.tol if options.tol else 1e-6
    iter_max = options.iter_max if options.iter_max else 100
    ode_method = 'LSODA' if options.stiff else 'RK45'

    # see _kb/06-solver-catalog.md (Fluid: "Python passage-time augmented
    # system") for why min_step is bounded to fail fast on stiff systems
    min_step = max(T * 1e-7, 1e-12)

    fullt = np.array([], dtype=float)
    fully = np.empty((0, total_phases_c), dtype=float)
    iter_count = 1
    finished = False
    tref = 0.0
    y_current = y0_c.copy()

    # Match MATLAB: integrate [0, T] per iteration, offset by tref
    while iter_count <= iter_max and not finished:
        try:
            sol = solve_ivp(
                ode_rhs,
                [0.0, T],
                y_current,
                method=ode_method,
                rtol=tol,
                atol=tol,
                min_step=min_step,
            )
            t_iter = sol.t
            y_iter = sol.y.T
            if not sol.success or not np.all(np.isfinite(y_iter)):
                # Stiff failure: no valid distribution, the caller falls back on the mean
                if len(fullt) == 0:
                    raise RuntimeError('passage-time integration failed (stiff augmented system)')
                break
        except RuntimeError:
            raise
        except Exception:
            if len(fullt) == 0:
                raise RuntimeError('passage-time integration failed (stiff augmented system)')
            break

        iter_count += 1
        if len(fullt) == 0:
            fullt = t_iter + tref
            fully = y_iter
        else:
            fullt = np.concatenate([fullt, t_iter[1:] + tref])
            fully = np.vstack([fully, y_iter[1:]])

        # Check if transient fluid is depleted (MATLAB: sum < 10e-10)
        if np.sum(np.maximum(y_iter[-1][transient_indices], 0)) < 1e-9:
            finished = True

        tref += t_iter[-1]
        y_current = y_iter[-1]

    if len(fullt) == 0:
        raise RuntimeError('passage-time integration produced no trajectory')

    # Compute CDF: F(t) = 1 - transient_fluid(t) / initial_fluid
    transient_over_time = np.sum(np.maximum(fully[:, transient_indices], 0), axis=1)
    cdf = 1.0 - transient_over_time / fluid_c

    # Adaptive CDF Refinement, the rule of MATLAB solver_fluid_passage_time.m,
    # the C++ fluid_passage_time and the JAR SolverFluid.passageTime: while some
    # adjacent pair of CDF values differs by more than max_cdf_jump, split EVERY
    # offending interval and re-integrate the whole curve on the new grid in one
    # call. Refining ONE interval per round instead spends the round cap on five
    # intervals and leaves the jump target unmet, and the curve is read back by
    # quadrature: on cdf_respt_closed_threeclasses a 130-point grid over [0,200]
    # made the right-endpoint mean read 1.1366 for an exactly Exp(1) response
    # time. Both caps bound WORK, not accuracy.
    if fluid_c > 0:
        max_cdf_jump = 0.0005
        max_refinement_rounds = 5
        max_points = 20001
        n_refined = 20

        for _round in range(max_refinement_rounds):
            if len(fullt) >= max_points:
                break
            dcdf = np.diff(cdf)
            dt = np.diff(fullt)
            # a jump at equal times is an ATOM of the law, not a resolution
            # failure: its linspace points would all be the same instant
            offending = np.where((dcdf > max_cdf_jump) & (dt > 0))[0]
            if offending.size == 0:
                break
            pieces = []
            for j in range(len(fullt) - 1):
                pieces.append(np.array([fullt[j]]))
                if j in offending:
                    extra = np.linspace(fullt[j], fullt[j + 1], n_refined)
                    pieces.append(extra[1:-1])
            pieces.append(np.array([fullt[-1]]))
            new_t = np.unique(np.concatenate(pieces))
            if new_t.size > max_points or new_t.size < 2:
                break
            try:
                sol_ref = solve_ivp(
                    ode_rhs, [new_t[0], new_t[-1]], y0_c, t_eval=new_t,
                    method=ode_method, rtol=tol, atol=tol,
                )
                if not sol_ref.success or sol_ref.t.size != new_t.size:
                    break
            except Exception:
                break
            fullt = sol_ref.t
            fully = np.maximum(sol_ref.y.T, 0.0)
            transient_over_time = np.sum(np.maximum(fully[:, transient_indices], 0), axis=1)
            cdf = 1.0 - transient_over_time / fluid_c

        # Horizon extension - extend while the TAIL misses the law.
        # The test is on the LAST grid point. It used to read the FIRST,
        # cdf[0], which is the CDF at the start of the horizon: the marked
        # class holds all of fluid_c at t=0 by construction, so that value is 0
        # whatever the horizon is, and lengthening the horizon cannot move it.
        # The loop it guarded was unreachable, and reachable only into harm --
        # its body REPLACED the refined curve with a fresh solve over a new
        # horizon, discarding the grid the refinement rounds above had just
        # paid for. Same fix as MATLAB solver_fluid_passage_time.m and the JAR
        # SolverFluid.
        max_extend_iterations = 10
        extend_iter = 0
        while cdf[-1] < 0.99 and extend_iter < max_extend_iterations:
            extend_iter += 1
            # CONTINUE the same trajectory from where it stopped and APPEND, as
            # the window loop above does: y_current is the end state and tref
            # the elapsed time, and ode_rhs is autonomous, so [0, extended_T]
            # from it is the next stretch of the SAME passage. Doubling each
            # round reaches a 1024x horizon within the cap instead of 11x.
            extended_T = T * (2 ** extend_iter)
            try:
                sol_ext = solve_ivp(
                    ode_rhs, [0.0, extended_T], y_current,
                    method=ode_method, rtol=tol, atol=tol,
                )
                if not sol_ext.success:
                    break
            except Exception:
                break
            t_ext = sol_ext.t
            y_ext = sol_ext.y.T
            if t_ext.size < 2 or not np.all(np.isfinite(y_ext)):
                break
            # drop the duplicated first row: it repeats the instant the curve
            # already ends on
            fullt = np.concatenate([fullt, t_ext[1:] + tref])
            fully = np.vstack([fully, np.maximum(y_ext[1:], 0.0)])
            tref += t_ext[-1]
            y_current = y_ext[-1]
            transient_over_time = np.sum(np.maximum(fully[:, transient_indices], 0), axis=1)
            cdf = 1.0 - transient_over_time / fluid_c

    return fullt, cdf


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

            # W is transposed relative to D0: D0[from,to] but W[to,from]
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
