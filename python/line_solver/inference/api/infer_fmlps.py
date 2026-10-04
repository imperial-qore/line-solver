import numpy as np
from scipy.optimize import minimize
from scipy.integrate import solve_ivp

from line_solver import SchedStrategy

from line_solver.inference.api.infer_fluid_ps_rt_likelihood import infer_fluid_ps_rt_likelihood


def infer_fmlps(model, node, rt, cls, ql, W):
    """FMLPS demand estimation using sn struct-level operations.

    Estimates service demands at a PS queue using the Fluid Maximum
    Likelihood for Processor Sharing method. Uses sn_set_service for
    fast parameter updates.

    Args:
        model: LINE Network model
        node: PS queue node
        rt: response time samples (n,)
        cls: class of each sample (n,), 1-based
        ql: queue lengths at arrival (n x R)
        W: total population

    Returns:
        demand_est: 1-D array of estimated demands (R,)

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """
    from line_solver.api.sn.transforms import sn_set_service

    rt = np.asarray(rt, dtype=float).flatten()
    cls = np.asarray(cls, dtype=int).flatten()
    ql = np.asarray(ql, dtype=float)

    sn = model.refreshStruct() or model.getStruct()
    R = sn.nclasses
    V = node.getNumberOfServers()
    st_idx = node.getStationIndex() - 1  # 0-based

    x_lb = np.full(R, np.min(rt) / W)
    x_ub = np.full(R, np.max(rt))

    # Initial point estimate
    mean_ql = np.mean(np.sum(ql, axis=1))
    v_tilde = min(mean_ql, V)
    x0 = np.zeros(R)
    for j in range(R):
        mask = cls == (j + 1)
        if np.any(mask):
            x0[j] = v_tilde * np.mean(rt[mask]) / mean_ql
        else:
            x0[j] = x_lb[j]

    # Cache station indices
    delay_idx = -1
    ref_idx = -1
    for ii in range(sn.nstations):
        if sn.sched[ii] == SchedStrategy.INF or sn.sched[ii] == SchedStrategy.INF.value:
            delay_idx = ii
        else:
            ref_idx = ii

    # Cache delay rates
    delay_rates = np.zeros(R)
    for kk in range(R):
        delay_rates[kk] = sn.rates[delay_idx, kk]

    def objfun(x):
        TOL = 1e-6

        # Update service rates
        for r in range(R):
            sn_set_service(sn, st_idx, r, 1.0 / x[r], 1.0)

        unique_tc = np.unique(cls)
        ftemp = np.zeros(len(rt))

        for tc in unique_tc:
            mask = cls == tc
            rt_tc = rt[mask]
            ql_tc = ql[mask, :]

            # Build augmented model and ODE handle
            ode_h, q_idx, aug_phases, _ = infer_fluid_ps_rt_likelihood(sn, int(tc))

            ftemp_tc = np.zeros(len(rt_tc))
            for rr in range(len(rt_tc)):
                like = _solve_fmlps_sample(
                    ode_h, q_idx, aug_phases,
                    delay_rates, delay_idx, ref_idx, R,
                    ql_tc[rr, :], W, int(tc), rt_tc[rr])
                ftemp_tc[rr] = np.log(TOL + like)
            ftemp[mask] = ftemp_tc

        return -np.sum(ftemp)

    bounds = list(zip(x_lb, x_ub))
    result = minimize(objfun, x0, method='L-BFGS-B', bounds=bounds,
                      options={'maxiter': 10000, 'ftol': 1e-10})
    return result.x


def _solve_fmlps_sample(ode_h, q_indices, aug_phases,
                         delay_rates, delay_idx, ref_idx, K,
                         a_queue, W, tagged_class, r_sampled):
    """Solve the fluid ODE for a single sample and extract likelihood."""
    new_k = K + 1
    new_fluid = 1.0

    # Compute initial fluid levels
    y0_levels = np.zeros((aug_phases.shape[0], K))

    # Delay station: distribute remaining fluid
    delay_jobs = (W - np.sum(a_queue)) * delay_rates / np.sum(delay_rates)
    y0_levels[delay_idx, :] = delay_jobs

    # Queue station: observed queue lengths
    y0_levels[ref_idx, :] = a_queue

    # Build state vector
    total_phases = int(np.sum(aug_phases))
    y0 = np.zeros(total_phases)
    M = aug_phases.shape[0]
    for i in range(M):
        for k in range(K):
            if aug_phases[i, k] > 0:
                y0[q_indices[i, k]] = y0_levels[i, k]

    # Move fluid from tagged class to Tagged at ref node
    tc_idx = tagged_class - 1  # 0-based
    y0[q_indices[ref_idx, tc_idx]] -= new_fluid
    y0[q_indices[ref_idx, new_k - 1]] = new_fluid

    # Solve ODE
    ref_tag_idx = q_indices[ref_idx, new_k - 1]

    def event_fn(t, y):
        return y[ref_tag_idx]
    event_fn.terminal = True
    event_fn.direction = 0

    sol = solve_ivp(ode_h, [0, r_sampled], y0, method='BDF',
                    events=event_fn, rtol=1e-5, atol=1e-8,
                    max_step=r_sampled / 10)

    if r_sampled <= sol.t[-1]:
        last_state = sol.y[:, -1]
        last_rates = ode_h(sol.t[-1], last_state)
        like = -last_rates[ref_tag_idx] / new_fluid
        return max(0.0, float(like))
    else:
        return 0.0
