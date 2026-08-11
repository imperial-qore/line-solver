import numpy as np
from scipy.integrate import solve_ivp

from line_solver import SchedStrategy


def infer_fluid_ps_rt_likelihood(sn, tagged_class, y0_levels=None, r_sampled=None):
    """Fluid-based response time likelihood.

    Builds an augmented model with tagged class K+1 by expanding the
    NetworkStruct arrays, then constructs the ODE right-hand side.

    Args:
        sn: NetworkStruct (from model.get_struct())
        tagged_class: class index of the tagged job (1-based)
        y0_levels: M x K matrix of fluid levels (optional)
        r_sampled: observed response time (optional)

    Returns:
        ode_h: ODE right-hand side function
        q_indices: M x Kc index matrix
        aug_phases: M x Kc phase count matrix
        LIKE: likelihood value (only if y0_levels and r_sampled provided)

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """
    M = sn.nstations
    K = sn.nclasses
    N = sn.nclosedjobs if hasattr(sn, 'nclosedjobs') else int(np.sum(sn.njobs[sn.njobs < np.inf]))
    Kc = K + 1

    # Build mu/phi from sn.rates if sn.mu is None (Python LINE compatibility)
    if sn.mu is not None:
        Lambda = sn.mu
        Pi = sn.phi
    else:
        Lambda = [[None] * K for _ in range(M)]
        Pi = [[None] * K for _ in range(M)]
        for j in range(M):
            for k in range(K):
                rate = sn.rates[j, k]
                if rate > 0 and not np.isnan(rate):
                    Lambda[j][k] = np.array([rate])
                    Pi[j][k] = np.array([1.0])
    rt_matrix = sn.rt
    S = sn.nservers.copy()

    # Find reference station and handle inf servers
    ref_idx = -1
    for i in range(M):
        if sn.sched[i] != SchedStrategy.INF and sn.sched[i] != SchedStrategy.INF.value:
            ref_idx = i
        if np.isinf(S[i]):
            S[i] = N

    # Build augmented service rates
    tc = tagged_class - 1  # 0-based
    new_mu = [[None] * Kc for _ in range(M)]
    new_pi = [[None] * Kc for _ in range(M)]
    for j in range(M):
        for k in range(K):
            new_mu[j][k] = Lambda[j][k] if Lambda[j][k] is not None else None
            new_pi[j][k] = Pi[j][k] if Pi[j][k] is not None else None
        new_mu[j][Kc - 1] = Lambda[j][tc] if Lambda[j][tc] is not None else None
        new_pi[j][Kc - 1] = Pi[j][tc] if Pi[j][tc] is not None else None

    # Handle NaN/disabled entries
    for j in range(M):
        for c in range(Kc):
            if new_mu[j][c] is not None and np.isscalar(new_mu[j][c]) and np.isnan(new_mu[j][c]):
                new_mu[j][c] = None
                new_pi[j][c] = None

    # Compute augmented phases
    aug_phases = np.zeros((M, Kc), dtype=int)
    for i in range(M):
        for c in range(Kc):
            if new_mu[i][c] is not None:
                mu_val = np.atleast_1d(new_mu[i][c])
                aug_phases[i, c] = len(mu_val)

    # Compute q_indices: maps (station, class) -> index in state vector
    q_indices = np.zeros((M, Kc), dtype=int)
    idx = 0
    for i in range(M):
        for c in range(Kc):
            q_indices[i, c] = idx
            if aug_phases[i, c] > 0:
                idx += aug_phases[i, c]

    total_phases = int(np.sum(aug_phases))

    # Build expanded routing table
    new_rt = np.zeros((M * Kc, M * Kc))
    for l in range(K):
        for m in range(K):
            for i in range(M):
                for j in range(M):
                    new_rt[i * Kc + l, j * Kc + m] = rt_matrix[i * K + l, j * K + m]

    # Tagged class routes like tagged_class
    for i in range(M):
        for j in range(M):
            new_rt[i * Kc + (Kc - 1), j * Kc + (Kc - 1)] = rt_matrix[i * K + tc, j * K + tc]

    # Absorption at ref_idx: tagged -> original classes
    chain_idx = -1
    if hasattr(sn, 'chains'):
        for ch in range(sn.chains.shape[0]):
            if sn.chains[ch, tc] == 1:
                chain_idx = ch
                break
    if chain_idx >= 0:
        classes_in_chain = np.where(sn.chains[chain_idx, :] == 1)[0]
        for l in classes_in_chain:
            for j in range(M):
                new_rt[ref_idx * Kc + (Kc - 1), j * Kc + l] = rt_matrix[ref_idx * K + tc, j * K + l]
    # Zero out tagged->tagged at ref_idx
    for j in range(M):
        new_rt[ref_idx * Kc + (Kc - 1), j * Kc + (Kc - 1)] = 0

    # Build ODE RHS for fluid model
    def ode_rhs(t, y):
        dydt = np.zeros(total_phases)
        for i in range(M):
            n_at_i = 0
            for c in range(Kc):
                if aug_phases[i, c] > 0:
                    n_at_i += y[q_indices[i, c]]

            capacity = S[i]
            for c in range(Kc):
                if aug_phases[i, c] > 0:
                    qi = q_indices[i, c]
                    mu_val = np.atleast_1d(new_mu[i][c])[0] if new_mu[i][c] is not None else 0
                    phi_val = np.atleast_1d(new_pi[i][c])[0] if new_pi[i][c] is not None else 1.0

                    if n_at_i > 0 and mu_val > 0:
                        # PS rate: mu * min(n, S) / n for the class share
                        share = y[qi] / max(n_at_i, 1e-10)
                        effective_rate = mu_val * min(n_at_i, capacity) * share

                        # Departures
                        dydt[qi] -= effective_rate * phi_val

                        # Route to other stations/classes
                        for j in range(M):
                            for d in range(Kc):
                                p_cd = new_rt[i * Kc + c, j * Kc + d]
                                if p_cd > 0 and aug_phases[j, d] > 0:
                                    qj = q_indices[j, d]
                                    dydt[qj] += effective_rate * phi_val * p_cd
        return dydt

    LIKE = np.nan
    if y0_levels is not None and r_sampled is not None:
        new_fluid = 1.0

        total_p = int(np.sum(aug_phases))
        y0 = np.zeros(total_p)

        for i in range(M):
            for k in range(K):
                if aug_phases[i, k] > 0:
                    y0[q_indices[i, k]] = y0_levels[i, k]

        y0[q_indices[ref_idx, tc]] -= new_fluid
        y0[q_indices[ref_idx, Kc - 1]] = new_fluid

        ref_tag_idx = q_indices[ref_idx, Kc - 1]

        def event_fn(t, y):
            return y[ref_tag_idx]
        event_fn.terminal = True
        event_fn.direction = 0

        sol = solve_ivp(ode_rhs, [0, r_sampled], y0, method='BDF',
                        events=event_fn, rtol=1e-5, atol=1e-8)

        if r_sampled <= sol.t[-1]:
            last_state = sol.y[:, -1]
            last_rates = ode_rhs(sol.t[-1], last_state)
            LIKE = -last_rates[ref_tag_idx] / new_fluid
            LIKE = max(0.0, LIKE)
        else:
            LIKE = 0.0

    return ode_rhs, q_indices, aug_phases, LIKE
