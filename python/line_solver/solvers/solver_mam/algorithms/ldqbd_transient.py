"""
Transient analysis of open queues via standard QBD (matrix exponentiation).

Supports single-class open models (Source -> Queue -> Sink).
For M/M/c: scalar QBD levels (any number of servers).
For M/PH/1: m-phase QBD levels (single server only).

Infinite capacity: truncation to finite buffer + expm.
Finite capacity: direct matrix exponentiation (expm).

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np
from scipy.linalg import expm


def solver_mam_ldqbd_transient(sn, options):
    """Transient QBD analysis for single-class open queues.

    Args:
        sn: NetworkStruct
        options: SolverOptions with timespan, tol

    Returns:
        (Qt, Ut, Tt): Each is M x K list of lists. Non-None entries are
        numpy arrays of shape (nTimePoints, 2) with [metric, time] columns.
    """
    M = sn.nstations
    K = sn.nclasses

    if K != 1:
        raise ValueError("Transient QBD method requires a single-class model.")

    njobs = np.asarray(sn.njobs).flatten()
    if not np.isinf(njobs[0]):
        raise ValueError("Transient QBD method requires an open model.")

    # Identify stations
    source_idx = -1
    queue_idx = -1
    sched = sn.sched if isinstance(sn.sched, dict) else {}
    for i in range(M):
        s = sched.get(i, None)
        if s is None:
            continue
        sname = s.name if hasattr(s, 'name') else str(s)
        sval = s.value if hasattr(s, 'value') else int(s)
        if sname == 'EXT' or sval == 16:
            source_idx = i
        elif sname == 'FCFS' or sval == 0:
            queue_idx = i

    if source_idx < 0 or queue_idx < 0:
        raise ValueError("Transient QBD method requires exactly one Source and one FCFS Queue.")

    # Extract parameters
    rates = np.asarray(sn.rates)
    lam = rates[source_idx, 0]

    # Service process - handle various Python native proc formats
    proc = sn.proc
    if isinstance(proc, list):
        ph_queue = proc[queue_idx][0] if isinstance(proc[queue_idx], list) else proc[queue_idx]
    elif isinstance(proc, dict):
        ph_queue = proc.get((queue_idx, 0), proc.get(queue_idx, None))
    else:
        ph_queue = proc[queue_idx]

    # Convert proc entry to D0/D1 matrices
    if isinstance(ph_queue, dict) and 'rate' in ph_queue:
        # Simple exponential: {'rate': mu}
        rate_val = ph_queue['rate']
        D0 = np.array([[-rate_val]])
        D1 = np.array([[rate_val]])
    elif isinstance(ph_queue, (list, tuple)):
        D0 = np.atleast_2d(np.asarray(ph_queue[0], dtype=float))
        D1 = np.atleast_2d(np.asarray(ph_queue[1], dtype=float)) if len(ph_queue) > 1 else None
    elif isinstance(ph_queue, dict) and 0 in ph_queue:
        D0 = np.atleast_2d(np.asarray(ph_queue[0], dtype=float))
        D1 = np.atleast_2d(np.asarray(ph_queue[1], dtype=float)) if 1 in ph_queue else None
    else:
        D0 = np.atleast_2d(np.asarray(ph_queue, dtype=float))
        D1 = None

    nservers_arr = np.asarray(sn.nservers).flatten()
    n_servers = int(nservers_arr[queue_idx])

    cap_arr = np.asarray(sn.cap).flatten()
    buf_cap = cap_arr[queue_idx]

    # Determine if exponential or PH
    if D0.shape[0] == 1 and D0.shape[1] == 1:
        mu = -D0[0, 0]
        n_phases = 1
        is_ph = False
    else:
        n_phases = D0.shape[0]
        is_ph = True
        mu = 0.0

    if is_ph and n_servers > 1:
        raise ValueError("Transient QBD with PH service supports single-server only.")

    # PH-specific
    if is_ph:
        alpha = _map_pie(D0, D1)
        t_exit = -D0 @ np.ones((n_phases, 1))
    else:
        alpha = None
        t_exit = None

    # Time parameters
    T_start = options.timespan[0]
    T_end = options.timespan[1]
    T_duration = T_end - T_start

    # Determine capacity
    if np.isinf(buf_cap) or buf_cap > 1000000:
        if is_ph:
            mean_svc = _map_mean(D0, D1)
            mu_eff = 1.0 / mean_svc if mean_svc > 0 else 1.0
        else:
            mu_eff = mu
        rho = lam / (n_servers * mu_eff)
        tol = options.tol if hasattr(options, 'tol') and options.tol > 0 else 1e-6
        if rho < 1.0:
            Cap = min(10000, max(200, int(np.ceil(-np.log(tol) / (-np.log(rho))))))
        else:
            Cap = 10000
    else:
        Cap = int(buf_cap)

    # Build generator Q
    if not is_ph:
        dim = Cap + 1
        Q = np.zeros((dim, dim))
        for n in range(Cap + 1):
            dep = min(n, n_servers) * mu
            arr = lam if n < Cap else 0.0
            if n > 0:
                Q[n, n - 1] = dep
            if n < Cap:
                Q[n, n + 1] = arr
            Q[n, n] = -(dep + arr)
    else:
        dim = 1 + Cap * n_phases
        Q = np.zeros((dim, dim))

        # Level 0
        Q[0, 0] = -lam
        Q[0, 1:1 + n_phases] = lam * alpha

        for n in range(1, Cap + 1):
            rs = 1 + (n - 1) * n_phases
            re = rs + n_phases
            rows = slice(rs, re)

            # Internal transitions (D0)
            if n < Cap:
                Q[rs:re, rs:re] = D0 - lam * np.eye(n_phases)
            else:
                Q[rs:re, rs:re] = D0

            # Arrivals: level n -> n+1
            if n < Cap:
                ns = rs + n_phases
                Q[rs:re, ns:ns + n_phases] = lam * np.eye(n_phases)

            # Departures: level n -> n-1
            if n == 1:
                Q[rs:re, 0:1] = t_exit
            else:
                ps = rs - n_phases
                Q[rs:re, ps:ps + n_phases] = D1

    # Time grid
    n_time_points = min(101, max(11, round(T_duration * 10)))
    dt = T_duration / (n_time_points - 1)
    times = np.linspace(T_start, T_end, n_time_points)

    # Initial distribution: empty queue
    pi_t = np.zeros((1, dim))
    pi_t[0, 0] = 1.0

    # Compute eQdt = expm(Q * dt)
    eQdt = expm(Q * dt)

    queue_lengths = np.zeros(n_time_points)
    util_values = np.zeros(n_time_points)
    tput_values = np.zeros(n_time_points)

    for t_idx in range(n_time_points):
        q = 0.0
        u = 0.0
        tput = 0.0

        for n in range(Cap + 1):
            if not is_ph:
                p_n = pi_t[0, n]
            else:
                if n == 0:
                    p_n = pi_t[0, 0]
                else:
                    idx_s = 1 + (n - 1) * n_phases
                    p_n = np.sum(pi_t[0, idx_s:idx_s + n_phases])

            q += n * p_n
            if n >= 1:
                u += (min(n, n_servers) / n_servers) * p_n
                if not is_ph:
                    tput += min(n, n_servers) * mu * p_n
                else:
                    idx_s = 1 + (n - 1) * n_phases
                    tput += pi_t[0, idx_s:idx_s + n_phases] @ t_exit

        queue_lengths[t_idx] = q
        util_values[t_idx] = u
        tput_values[t_idx] = float(tput)

        if t_idx < n_time_points - 1:
            pi_t = pi_t @ eQdt

    # Package results as M x K lists
    Qt = [[None] * K for _ in range(M)]
    Ut = [[None] * K for _ in range(M)]
    Tt = [[None] * K for _ in range(M)]

    Qt[queue_idx][0] = np.column_stack([queue_lengths, times])
    Ut[queue_idx][0] = np.column_stack([util_values, times])
    Tt[queue_idx][0] = np.column_stack([tput_values, times])

    return Qt, Ut, Tt


def _map_pie(D0, D1):
    """Compute stationary probability vector of the phase process."""
    S = D0 + D1
    n = S.shape[0]
    # Solve pi * S = 0, sum(pi) = 1
    A = S.T.copy()
    A[-1, :] = 1.0
    b = np.zeros(n)
    b[-1] = 1.0
    try:
        pi = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        pi = np.ones(n) / n
    # MAP initial (entry) probability vector: alpha = (pi @ D1) / sum(pi @ D1)
    piD1 = pi @ D1
    s = np.sum(piD1)
    if s > 0:
        return piD1 / s
    return np.ones(n) / n


def _map_mean(D0, D1):
    """Compute mean inter-event time of a MAP."""
    S = D0 + D1
    n = S.shape[0]
    A = S.T.copy()
    A[-1, :] = 1.0
    b = np.zeros(n)
    b[-1] = 1.0
    try:
        pi = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        pi = np.ones(n) / n
    # Mean = pi * (-D0)^{-1} * ones / (pi * D1 * ones)
    try:
        neg_D0_inv = np.linalg.inv(-D0)
    except np.linalg.LinAlgError:
        return 1.0
    ones = np.ones(n)
    num = pi @ neg_D0_inv @ ones
    den = pi @ D1 @ ones
    if den > 0:
        return num / den
    return 1.0
