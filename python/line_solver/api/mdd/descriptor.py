"""
Kronecker rate descriptor of a single-class closed queueing network.

For the Miner-Ciardo-Donatelli approximate-aggregation solver (mdd_mcd), after
A.S. Miner, G. Ciardo, S. Donatelli, "Using the exact state space of a Markov
model to compute approximate stationary measures", SIGMETRICS 2000.

The transition rate matrix is expressed compositionally as
    R = sum_e ( kron_k W_k^e )   restricted to the reachable set,
with W_k^e[i_k,j_k] = lambda_k^e[i_k] * Prob_k^e(i_k,j_k) (Eq. 1). Each level k
is a station and each event e is a completion at station a routed to b.

Exponential stations
--------------------
The local state is the population alone:
    W_a^e[i,i-1] = mu[a]*min(i,servers[a])*P[a,b]   (i >= 1)   -- departure
    W_b^e[i,i+1] = 1                                (i <= N-1) -- arrival
    W_l^e        = I                                (l not in {a,b})

Phase-type stations
-------------------
The local state is the PAIR (population, phase of the job in service), encoded
in one level rather than two. Splitting them does not work: on a completion
routed into station b the phase at b restarts only when b was empty, a joint
condition on b's two components, which is not a product of per-level terms.
Merging them keeps every event local. The encoding is
    index 0                   : station empty
    index 1 + (n-1)*h + (a-1) : n jobs present, job in service in phase a
so the domain is 1 + N*h and h = 1 reproduces the exponential encoding
index = n exactly. With exit vector t = D1*1 and entry law pie,
    departure  (n,a) -> (n-1,b)  at t[a]*P[a,b]*pie[b]   for n >= 2
               (1,a) -> 0        at t[a]*P[a,b]
    arrival    0     -> (1,b)    at pie[b]
               (m,c) -> (m+1,c)  at 1                    for m >= 1
    internal   (n,a) -> (n,b)    at D0[a,b]              for n >= 1, a != b

Restrictions
------------
A phase-type station must be single-server. With c > 1 or an infinite server the
local state would have to count jobs per phase rather than name one phase, which
is a different and much larger encoding.

A phase-type station must also be NON-preemptive. The composite level names the
phase of the one job in service and restarts it at pie when the next job starts;
under preemptive resume an arrival suspends that job and its phase has to be
remembered, so the local state would need a stack of phases. This matters for
LCFSPR, which is BCMP type 2 and stays product-form under general service: that
insensitivity is real but is NOT reachable through this encoding. Exponential
service is unaffected, preemption being immaterial by memorylessness. Pass
`sched` to have the case rejected rather than silently modelled as
non-preemptive.

See also: mdd_mcd, mdd_reachset, mdd_ps.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

from typing import Dict, List, Optional

import numpy as np
import scipy.sparse as sp

from ..io.logging import line_error

# Preemptive-resume and shared-server disciplines the composite encoding cannot
# represent. Compared by NAME because the SchedStrategy numeric values are not
# shared across the codebases.
_PREEMPTIVE_RESUME = ('LCFSPR', 'FCFSPR', 'LCFSPRPRIO', 'FCFSPRPRIO')
_SHARED_SERVER = ('PS', 'DPS', 'GPS')


def _sched_name(s) -> str:
    """Normalise a discipline to its upper-case name."""
    if s is None:
        return ''
    if hasattr(s, 'name'):
        return str(s.name).upper()
    return str(s).upper().split('.')[-1]


def _service_law(pr, mu_i, i):
    """(D0, D1, entry law, phases) of station i from a proc entry."""
    if pr is None or (isinstance(pr, (list, tuple)) and len(pr) == 0):
        return np.array([[-float(mu_i)]]), np.array([[float(mu_i)]]), np.array([1.0]), 1
    if isinstance(pr, (list, tuple)):
        D0 = np.asarray(pr[0], dtype=float)
        D1 = np.asarray(pr[1], dtype=float)
    else:
        getter = getattr(pr, 'getProcess', None) or getattr(pr, 'get_process', None)
        if getter is None:
            line_error('mdd_descriptor',
                       'station %d carries a service law that exposes neither a (D0,D1) pair nor '
                       'a getProcess method' % (i + 1))
        pc = getter()
        D0 = np.asarray(pc[0], dtype=float)
        D1 = np.asarray(pc[1], dtype=float)
    return D0, D1, None, D0.shape[0]


def _entry_law(pr, pie_opt, h, i, D1=None, caller='mdd_descriptor'):
    """Entry law of a phase-type station, from options, the object, or D1.

    A (D0,D1) pair carries its own restart law: for a renewal process
    D1 = t0 * pie, so every row with a positive exit rate is proportional to
    pie. Deriving it is not optional -- defaulting to e_1 instead silently
    replaces a hyperexponential (whose D0 is diagonal, so a job entering phase 1
    can never leave it) by an exponential at the phase-1 rate.
    """
    if pie_opt is not None and len(np.ravel(np.asarray(pie_opt))) > 0:
        pie = np.asarray(pie_opt, dtype=float).ravel()
        return pie / pie.sum()
    if not isinstance(pr, (list, tuple)):
        getter = getattr(pr, 'getInitProb', None) or getattr(pr, 'get_init_prob', None)
        if getter is not None:
            pie = np.asarray(getter(), dtype=float).ravel()
            return pie / pie.sum()
    if D1 is not None:
        t0 = np.asarray(D1, dtype=float).sum(axis=1)
        live = np.nonzero(t0 > 0)[0]
        if live.size:
            rows = np.asarray(D1, dtype=float)[live, :] / t0[live, None]
            # A non-renewal MAP restarts in a law that depends on the phase it
            # left from, which one entry vector cannot express; the composite
            # level would silently model the renewal process instead.
            if live.size > 1 and np.max(np.abs(rows - rows[0, :])) > 1e-9:
                line_error(caller,
                           'station %d carries a service law whose restart distribution depends on '
                           'the completing phase (a non-renewal MAP); the local state names one '
                           'entry law, so this encoding cannot represent it' % (i + 1))
            return rows[0, :]
    pie = np.zeros(h)
    pie[0] = 1.0
    return pie


def mdd_descriptor(mu, P, servers, N, options: Optional[Dict] = None) -> Dict:
    """Kronecker rate descriptor of a single-class closed queueing network.

    Parameters
    ----------
    mu : length-K station service rates; entry i is ignored when station i is
        given a phase-type law through options['proc']
    P : K x K routing matrix (row-stochastic)
    servers : length-K servers per station (inf for delay/IS)
    N : closed population
    options : dict with optional keys
        proc  - length-K list; entry i None for an exponential station, or a
                Markovian distribution object, or a (D0,D1) pair
        pie   - length-K list of entry laws; taken from the object, or e_1, when
                omitted
        sched - length-K list of SchedStrategy values. Optional, and only
                consulted to REJECT a phase-type law at a preemptive-resume or
                shared-server station. Supply it whenever the stations are not
                all non-preemptive, because the descriptor otherwise has no way
                to detect that case.

    Returns
    -------
    dict with keys K, N, domain, mu, servers, P, nphases, valuemap, init,
    nextfun, events.
    """
    if options is None:
        options = {}
    proc = list(options.get('proc') or [])
    pie_opt = list(options.get('pie') or [])
    sched = list(options.get('sched') or [])

    mu = np.asarray(mu, dtype=float).ravel()
    servers = np.asarray(servers, dtype=float).ravel()
    P = np.asarray(P, dtype=float)
    N = int(N)
    K = mu.size

    # ---- per-station service law
    D0: List[np.ndarray] = [None] * K
    D1: List[np.ndarray] = [None] * K
    pie: List[np.ndarray] = [None] * K
    h = np.ones(K, dtype=int)
    for i in range(K):
        pr = proc[i] if i < len(proc) else None
        D0[i], D1[i], pie_i, h[i] = _service_law(pr, mu[i], i)
        if h[i] == 1 and (pr is None or (isinstance(pr, (list, tuple)) and len(pr) == 0)):
            pie[i] = np.array([1.0])
            continue
        pie[i] = _entry_law(pr, pie_opt[i] if i < len(pie_opt) else None, int(h[i]), i,
                            D1=D1[i], caller='mdd_descriptor')
        if h[i] > 1 and servers[i] != 1:
            line_error('mdd_descriptor',
                       'station %d has a phase-type service law and %g servers; a multi-server or '
                       'delay station would have to count jobs per phase rather than name the '
                       'phase of one job in service, which this encoding does not carry'
                       % (i + 1, servers[i]))
        # The composite level names ONE in-service phase and restarts it at pie
        # when the next job starts, i.e. NON-preemptive service. Under
        # preemptive-resume an arrival suspends the job in service and its phase
        # must be remembered, so the local state would need a STACK of phases.
        if h[i] > 1 and i < len(sched) and sched[i] is not None:
            nm = _sched_name(sched[i])
            if nm in _PREEMPTIVE_RESUME:
                line_error('mdd_descriptor',
                           'station %d combines a phase-type service law with a preemptive-resume '
                           'discipline; the suspended jobs\' phases would have to be stacked in the '
                           'local state, which this encoding does not carry, and the descriptor '
                           'would silently model the non-preemptive chain instead' % (i + 1))
            # Under a shared server every job present is in service and holds its
            # own phase, so naming one in-service phase is the wrong state.
            if nm in _SHARED_SERVER:
                line_error('mdd_descriptor',
                           'station %d combines a phase-type service law with a shared-server '
                           'discipline; every job present is in service and holds its own phase, '
                           'which this encoding does not carry. Use mdd_ps, whose local state is '
                           'the per-phase count vector.' % (i + 1))

    d = 1 + N * h                              # local domain per station

    # ---- index maps
    valuemap = []
    for i in range(K):
        vm = np.zeros(int(d[i]))
        for n in range(1, N + 1):
            vm[(1 + (n - 1) * h[i]):(1 + n * h[i])] = n
        valuemap.append(vm)

    # ---- initial state: all jobs at station 1, in its entry phase
    init = np.zeros(K, dtype=int)
    init[0] = _idx(N, int(np.nonzero(pie[0] > 0)[0][0]) + 1, int(h[0]))

    desc = {
        'K': K,
        'N': N,
        'domain': d.astype(int),
        'mu': mu,
        'servers': servers,
        'P': P,
        'nphases': h,
        'valuemap': valuemap,
        'init': init,
    }
    desc['nextfun'] = lambda s: _next(s, D0, D1, pie, h, N, P)

    # ---- events
    events = []
    # internal phase changes, one event per phase-type station
    for i in range(K):
        if h[i] == 1:
            continue
        Wi = _internal(D0[i], int(h[i]), N, int(d[i]))
        if Wi.nnz == 0:
            continue
        events.append({'a': i, 'b': i, 'lev': [i], 'W': [Wi]})

    # completions routed a -> b
    aa, bb = np.nonzero(P)
    for a, b in zip(aa.tolist(), bb.tolist()):
        if a == b:
            continue
        pr_ab = float(P[a, b])
        Wa = _departure(D1[a], pie[a], int(h[a]), N, int(d[a]), float(mu[a]),
                        float(servers[a]), pr_ab)
        Wb = _arrival(pie[b], int(h[b]), N, int(d[b]))
        events.append({'a': a, 'b': b, 'lev': [a, b], 'W': [Wa, Wb]})
    desc['events'] = events
    return desc


# ---------------------------------------------------------------------------
def _idx(n, a, h):
    """Local index of (population n, service phase a); 0 when empty."""
    if n == 0:
        return 0
    return 1 + (n - 1) * h + (a - 1)


def _decode(idx, h):
    """(population, phase) of a local index; (0,0) when empty."""
    if idx == 0:
        return 0, 0
    return (idx - 1) // h + 1, (idx - 1) % h + 1


# ---------------------------------------------------------------------------
def _internal(D0i, h, N, d):
    """Phase changes that do not complete a service, at any population n >= 1."""
    ri, ci, vv = [], [], []
    for n in range(1, N + 1):
        for a in range(1, h + 1):
            for b in range(1, h + 1):
                if a == b or D0i[a - 1, b - 1] == 0:
                    continue
                ri.append(_idx(n, a, h))
                ci.append(_idx(n, b, h))
                vv.append(D0i[a - 1, b - 1])
    return sp.coo_matrix((vv, (ri, ci)), shape=(d, d)).tocsr()


# ---------------------------------------------------------------------------
def _departure(D1i, piei, h, N, d, mui, srv, pr):
    """Completion at this station, routed out with probability pr."""
    ri, ci, vv = [], [], []
    if h == 1:
        # exponential: the multi-server and delay rate laws live here
        for n in range(1, N + 1):
            ri.append(n)
            ci.append(n - 1)
            vv.append(mui * min(n, srv) * pr)
    else:
        t = D1i.sum(axis=1)                      # exit rate per phase
        for n in range(1, N + 1):
            for a in range(1, h + 1):
                if t[a - 1] == 0:
                    continue
                if n == 1:
                    ri.append(_idx(1, a, h))
                    ci.append(0)
                    vv.append(t[a - 1] * pr)
                else:
                    for b in range(1, h + 1):
                        if piei[b - 1] == 0:
                            continue
                        ri.append(_idx(n, a, h))
                        ci.append(_idx(n - 1, b, h))
                        vv.append(t[a - 1] * pr * piei[b - 1])
    return sp.coo_matrix((vv, (ri, ci)), shape=(d, d)).tocsr()


# ---------------------------------------------------------------------------
def _arrival(pieb, h, N, d):
    """An arrival starts service only when the station was empty."""
    ri, ci, vv = [], [], []
    if h == 1:
        for m in range(N):
            ri.append(m)
            ci.append(m + 1)
            vv.append(1.0)
    else:
        for b in range(1, h + 1):
            if pieb[b - 1] == 0:
                continue
            ri.append(0)
            ci.append(_idx(1, b, h))
            vv.append(pieb[b - 1])
        for m in range(1, N):
            for c in range(1, h + 1):
                ri.append(_idx(m, c, h))
                ci.append(_idx(m + 1, c, h))
                vv.append(1.0)
    return sp.coo_matrix((vv, (ri, ci)), shape=(d, d)).tocsr()


# ---------------------------------------------------------------------------
def _next(s, D0, D1, pie, h, N, P):
    """Successor local-index vectors of s, used to generate the reachable set."""
    K = len(s)
    T = []
    for i in range(K):
        ni, ai = _decode(int(s[i]), int(h[i]))
        if ni == 0:
            continue
        # internal phase change
        if h[i] > 1:
            for b in range(1, int(h[i]) + 1):
                if b == ai or D0[i][ai - 1, b - 1] == 0:
                    continue
                t = list(s)
                t[i] = _idx(ni, b, int(h[i]))
                T.append(t)
        # completion routed to j
        if h[i] == 1:
            exits = 1.0
        else:
            exits = float(D1[i].sum(axis=1)[ai - 1])
        if exits == 0:
            continue
        for j in range(K):
            if j == i or P[i, j] <= 0:
                continue
            nj, aj = _decode(int(s[j]), int(h[j]))
            newi = _idx(ni - 1, 1, 1) if h[i] == 1 else None
            for bi in range(1, int(h[i]) + 1):
                if h[i] > 1:
                    if ni == 1:
                        newi = 0
                    elif pie[i][bi - 1] == 0:
                        continue
                    else:
                        newi = _idx(ni - 1, bi, int(h[i]))
                elif bi > 1:
                    continue
                for bj in range(1, int(h[j]) + 1):
                    if nj == 0:
                        if pie[j][bj - 1] == 0:
                            continue
                        newj = _idx(1, bj, int(h[j]))
                    elif bj > 1:
                        continue
                    else:
                        newj = _idx(nj + 1, aj, int(h[j]))
                    t = list(s)
                    t[i] = newi
                    t[j] = newj
                    T.append(t)
                if ni == 1 and h[i] > 1:
                    break
    return T
