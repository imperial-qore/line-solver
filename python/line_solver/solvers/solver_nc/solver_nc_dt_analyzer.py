"""
Discrete-time (slotted) normalizing-constant analyzer.

Native Python port of matlab/src/solvers/NC/nc_is_dt_model.m and
solver_nc_dt_analyzer.m.

Two families of discrete-time product-form models are covered, both from
H. Daduna, Queueing Networks with Discrete Time Scale, LNCS 2046, Springer,
2001:

    chapter 2  a Bernoulli server fed by a Bernoulli arrival stream, with an
               unbounded buffer (theorem 2.3, corollary 2.7), a finite buffer
               (corollary 2.8) or a load-dependent service probability
               (example 2.10);
    chapter 3  a closed cycle of Bernoulli servers (theorem 3.2,
               corollary 3.4).

Every metric is on the slot lattice: a rate is a per-slot probability and a
time is a number of slots.
"""

import math
import time
from typing import Any, Dict, Optional

import numpy as np

from ...api.dpfqn import dpfqn_nc, dpfqn_ncld
from ...api.dqsys import dqsys_bernoulli1
from ...constants import ProcessType, SchedStrategy

__all__ = ['nc_is_dt_model', 'solver_nc_dt_analyzer', 'is_slotted']

_TOL = 1e-8
_COARSE = 1e-3


def is_slotted(options) -> bool:
    """True when the caller asked for the discrete-time route.

    Never inferred from the model: a Geometric service time is an ordinary
    continuous-time model unless the caller declares the slot lattice with
    ``options.config['slotted']``.
    """
    cfg = getattr(options, 'config', None)
    if isinstance(cfg, dict):
        return bool(cfg.get('slotted', False))
    if cfg is not None and hasattr(cfg, 'slotted'):
        return bool(getattr(cfg, 'slotted') or False)
    return False


def _slot_length(options) -> float:
    cfg = getattr(options, 'config', None)
    d = 1.0
    if isinstance(cfg, dict):
        d = cfg.get('slotlength', cfg.get('slot_length', 1.0))
    elif cfg is not None and hasattr(cfg, 'slotlength'):
        d = getattr(cfg, 'slotlength')
    if d is None:
        d = 1.0
    d = float(d)
    if not math.isfinite(d) or d <= 0:
        raise ValueError('options.config slotlength must be a positive finite scalar')
    return d


def nc_is_dt_model(sn, options=None) -> Dict[str, Any]:
    """Classify sn against the discrete-time product-form families.

    Returns a dict with key ``kind`` in ``{'bernoulli1', 'cycle', 'none'}``.
    For ``'none'`` the key ``reason`` says why, and the analyzer turns that
    into an error rather than falling back to a continuous-time approximation.

    The admissible feature set is narrow because the discrete-time product
    form is narrow. Beyond the geometric service requirement:

    - a cycle is the only topology; Daduna, section 4.1, records that general
      discrete-time topologies of FCFS Bernoulli servers have no product form;
    - every station must be a single server. Pestien and Ramakrishnan, quoted
      before example 2.10, proved that a multiserver node inside a cycle of
      geometrical queues destroys the product form for any finite server
      count;
    - class switching is rejected, and a multichain cycle is admitted only
      through the aggregate population.
    """
    dt: Dict[str, Any] = {'kind': 'none', 'reason': ''}
    if sn is None:
        dt['reason'] = 'no network structure available'
        return dt

    M = int(sn.nstations)
    R = int(sn.nclasses)
    rates = np.asarray(sn.rates, dtype=float).reshape(M, R)
    procid = np.asarray(sn.procid, dtype=object).reshape(M, R)

    for i in range(M):
        for r in range(R):
            if np.isfinite(rates[i, r]) and rates[i, r] > 0:
                if _name(procid[i, r]) != _name(ProcessType.GEOMETRIC):
                    dt['reason'] = ('station %d serves class %d with a %s process; a '
                                    'discrete-time model needs Geometric service and '
                                    'interarrival times' % (i, r, _name(procid[i, r])))
                    return dt

    csmask = getattr(sn, 'csmask', None)
    if csmask is not None and np.size(csmask) > 0:
        cs = np.array(csmask, dtype=bool).reshape(R, R).copy()
        np.fill_diagonal(cs, False)
        if cs.any():
            dt['reason'] = 'class switching is not covered by the discrete-time product form'
            return dt
    for name, label in (('cdscaling', 'class-dependent'), ('jdscaling', 'joint-dependent')):
        sc = getattr(sn, name, None)
        if sc is not None and any(s is not None for s in np.atleast_1d(np.asarray(sc, dtype=object)).ravel()):
            dt['reason'] = '%s scaling is not covered by the discrete-time product form' % label
            return dt

    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    if np.any(np.isinf(njobs)):
        return _dt_open_single(sn, dt, rates, njobs)
    return _dt_closed_cycle(sn, dt, rates, njobs)


def _dt_open_single(sn, dt, rates, njobs):
    """Chapter 2: Source -> Bernoulli server -> Sink, one open class."""
    M = int(sn.nstations)
    R = int(sn.nclasses)
    if R != 1:
        dt['reason'] = 'the discrete-time single-node route handles one open class'
        return dt
    srcs = [i for i in range(M) if _sched(sn, i) == _name(SchedStrategy.EXT)]
    if len(srcs) != 1:
        dt['reason'] = 'an open discrete-time model needs exactly one Source'
        return dt
    queues = [i for i in range(M) if i not in srcs]
    if len(queues) != 1:
        dt['reason'] = ('the discrete-time single-node route handles one queueing station, '
                        'found %d' % len(queues))
        return dt
    ist = queues[0]
    if _sched(sn, ist) != _name(SchedStrategy.FCFS):
        dt['reason'] = 'a Bernoulli server is a FCFS station'
        return dt
    c = float(np.asarray(sn.nservers, dtype=float).ravel()[ist])
    if math.isfinite(c) and c != 1:
        dt['reason'] = ('a Bernoulli server is a single-server station; use load dependence '
                        'for the multiserver approximation of example 2.10')
        return dt

    b = float(rates[srcs[0], 0])
    p = float(rates[ist, 0])
    if not math.isfinite(b) or b <= 0 or b > 1:
        dt['reason'] = 'the source arrival probability must lie in (0,1]'
        return dt
    if not math.isfinite(p) or p <= 0 or p > 1:
        dt['reason'] = 'the service probability must lie in (0,1]'
        return dt

    cap = _dt_capacity(sn, ist, 0)
    lld = _dt_lld(sn, ist)
    if math.isinf(cap) and lld is not None:
        dt['reason'] = ('a load-dependent Bernoulli server needs a finite capacity to bound '
                        'the state space')
        return dt
    if math.isinf(cap) and b >= p:
        dt['reason'] = ('an unbounded discrete-time queue needs an arrival probability below '
                        'the service probability')
        return dt

    if math.isinf(cap):
        svc = p
    else:
        svc = p * _dt_lld_vector(lld, int(cap))
        if np.any(svc <= 0) or np.any(svc > 1):
            dt['reason'] = 'load dependence must keep the service probability inside (0,1]'
            return dt

    dt.update({'kind': 'bernoulli1', 'station': ist, 'source': srcs[0],
               'arrivalProb': b, 'serviceProb': svc, 'capacity': cap})
    return dt


def _dt_closed_cycle(sn, dt, rates, njobs):
    """Chapter 3: a closed cycle of state dependent Bernoulli servers."""
    M = int(sn.nstations)
    R = int(sn.nclasses)
    N = float(np.sum(njobs[np.isfinite(njobs)]))
    if N <= 0 or N != math.floor(N):
        dt['reason'] = 'the closed population must be a positive integer'
        return dt
    N = int(N)

    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    for i in range(M):
        if _sched(sn, i) != _name(SchedStrategy.FCFS):
            dt['reason'] = ('station %d is not FCFS; a cycle of Bernoulli servers is FCFS '
                            'throughout' % i)
            return dt
        if math.isfinite(nservers[i]) and nservers[i] != 1:
            dt['reason'] = ('station %d has %d servers; a multiserver node inside a cycle of '
                            'geometrical queues has no product form' % (i, int(nservers[i])))
            return dt

    p = np.zeros(M)
    for i in range(M):
        rr = rates[i, :]
        rr = rr[np.isfinite(rr) & (rr > 0)]
        if rr.size == 0:
            dt['reason'] = 'station %d serves no class' % i
            return dt
        if rr.max() - rr.min() > _TOL * rr.max():
            dt['reason'] = ('station %d has a class-dependent service probability; the '
                            'discrete-time cycle needs one Bernoulli server per node' % i)
            return dt
        p[i] = rr[0]
        if p[i] <= 0 or p[i] >= 1:
            dt['reason'] = ('station %d has service probability %g; the product form of '
                            'theorem 3.2 needs p in (0,1)' % (i, p[i]))
            return dt

    order = _dt_cycle_order(sn)
    if order is None:
        dt['reason'] = ('the stations do not form a single deterministic cycle; discrete-time '
                        'FCFS networks of other topologies have no product form')
        return dt

    P = np.zeros((M, max(N, 1)))
    for k, i in enumerate(order):
        P[k, :] = p[i] * _dt_lld_vector(_dt_lld(sn, i), N)
    if np.any(P <= 0) or np.any(P > 1):
        dt['reason'] = 'load dependence must keep every service probability inside (0,1]'
        return dt

    dt.update({'kind': 'cycle', 'order': order, 'serviceProb': P, 'population': N})
    return dt


def _dt_cycle_order(sn) -> Optional[list]:
    """Station order along the cycle, or None when the routing is not a single
    deterministic cycle visiting every station exactly once."""
    from ...api.sn.transforms import sn_rt_stations
    M = int(sn.nstations)
    R = int(sn.nclasses)
    out = sn_rt_stations(sn)
    rtst = np.asarray(out[0] if isinstance(out, tuple) else out, dtype=float)

    succ = [0] * M
    for i in range(M):
        tgt = -1
        for j in range(M):
            w = float(rtst[i * R:(i + 1) * R, j * R:(j + 1) * R].sum())
            if w > _TOL:
                if abs(w - R) > _COARSE and abs(w - 1) > _COARSE:
                    return None       # fractional routing out of station i
                if tgt != -1:
                    return None       # more than one successor
                tgt = j
        if tgt == -1 or tgt == i:
            return None
        succ[i] = tgt

    visited = [False] * M
    order = []
    cur = 0
    for _ in range(M):
        if visited[cur]:
            return None
        visited[cur] = True
        order.append(cur)
        cur = succ[cur]
    if cur != 0 or not all(visited):
        return None
    return order


def _name(x) -> str:
    """Enum member name, or the plain string form.

    Two SchedStrategy/ProcessType enums coexist in the native package
    (line_solver.constants and line_solver.lang.base), and sn carries members
    of the lang.base one. Their integer values are aligned but the classes are
    distinct, so ``member == Other.MEMBER`` is False; the numeric ids also
    differ from MATLAB's. Comparing names is the only form that holds in both
    directions, and it is the rule CLAUDE.md states for cross-codebase enum
    comparisons.
    """
    return getattr(x, 'name', None) or str(x)


def _sched(sn, i: int) -> str:
    """Scheduling strategy name of station i. Python stores sn.sched as a dict
    keyed by station index, MATLAB and the JAR as a vector."""
    sc = getattr(sn, 'sched', None)
    if isinstance(sc, dict):
        return _name(sc.get(i, None))
    arr = np.atleast_1d(np.asarray(sc, dtype=object)).ravel()
    return _name(arr[i]) if i < arr.size else ''


def _dt_capacity(sn, ist: int, r: int) -> float:
    """Effective buffer capacity of station ist for class r, inf when unbounded."""
    cap = math.inf
    scap = getattr(sn, 'cap', None)
    if scap is not None and np.size(scap) > ist:
        cap = min(cap, float(np.asarray(scap, dtype=float).ravel()[ist]))
    ccap = getattr(sn, 'classcap', None)
    if ccap is not None and np.size(ccap) > 0:
        arr = np.asarray(ccap, dtype=float)
        if arr.ndim == 2 and arr.shape[0] > ist and arr.shape[1] > r:
            cap = min(cap, float(arr[ist, r]))
    return cap


def _dt_lld(sn, ist: int):
    """Load-dependent scaling vector of station ist, None when undeclared."""
    lld = getattr(sn, 'lldscaling', None)
    if lld is None or np.size(lld) == 0:
        return None
    arr = np.asarray(lld, dtype=float)
    if arr.ndim != 2 or arr.shape[0] <= ist:
        return None
    row = arr[ist, :]
    if not np.any(np.abs(row - 1.0) > _TOL):
        return None
    return row


def _dt_lld_vector(lld, N: int) -> np.ndarray:
    """Expand a load-dependent scaling to alpha(1..N), holding the last declared
    value beyond the tabulated range, as the load-dependent NC solvers do."""
    alpha = np.ones(N)
    if lld is None:
        return alpha
    lld = np.asarray(lld, dtype=float).ravel()
    n = min(N, lld.size)
    alpha[:n] = lld[:n]
    if N > lld.size:
        alpha[lld.size:] = lld[-1]
    return alpha


def solver_nc_dt_analyzer(sn, options):
    """Exact discrete-time analysis of sn.

    Returns ``(Q, U, R, T, C, X, lG, runtime, iter, method)`` in the contract
    of the other NC analyzers.

    On the cycle route the per-class split is proportional to the per-class
    population. Service in the cycle is type independent and FCFS forbids
    overtaking, so the cyclic order of the jobs is frozen; the marginal law of
    the queue lengths carries no class information, and the long-run share of
    station j held by chain g is its population share N_g/N. That is the sense
    in which section 3.2 of the reference calls the multichain case a direct
    adaptation of the unichain one.
    """
    tstart = time.time()
    dt = nc_is_dt_model(sn, options)
    if dt['kind'] == 'bernoulli1':
        method = 'dt.bernoulli1'
        Q, U, R, T, C, X, lG = _dt_single(sn, dt)
    elif dt['kind'] == 'cycle':
        P = dt['serviceProb']
        if P.shape[1] >= 1 and np.all(np.abs(P - P[:, [0]]) <= _TOL):
            method = 'dt.cycle'
        else:
            method = 'dt.cycleld'
        Q, U, R, T, C, X, lG = _dt_cycle(sn, dt, method)
    else:
        raise RuntimeError('options.config slotted was requested but the model is not a '
                           'discrete-time product-form model: %s' % dt['reason'])

    d = _slot_length(options)
    if d != 1.0:
        T = T / d
        X = X / d
        R = R * d
        C = C * d
    return Q, U, R, T, C, X, lG, time.time() - tstart, 1, method


def _dt_single(sn, dt):
    """Chapter 2 route: one Bernoulli server, one open class."""
    M = int(sn.nstations)
    K = int(sn.nclasses)
    Q = np.zeros((M, K)); U = np.zeros((M, K))
    R = np.zeros((M, K)); T = np.zeros((M, K))
    C = np.zeros(K); X = np.zeros(K)

    res = dqsys_bernoulli1(dt['arrivalProb'], dt['serviceProb'], dt['capacity'])
    ist = dt['station']
    Q[ist, 0] = res['meanQueueLength']
    U[ist, 0] = res['utilization']
    T[ist, 0] = res['throughput']
    R[ist, 0] = res['meanSojournTime']
    X[0] = res['throughput']
    C[0] = res['meanSojournTime']
    # The Source row carries the offered stream, as on the continuous-time route.
    T[dt['source'], 0] = dt['arrivalProb']
    return Q, U, R, T, C, X, math.log(res['normConst'])


def _dt_cycle(sn, dt, method):
    """Chapter 3 route: closed cycle of Bernoulli servers."""
    M = int(sn.nstations)
    K = int(sn.nclasses)
    N = dt['population']
    order = dt['order']
    P = dt['serviceProb']

    Qs = np.zeros(M); Us = np.zeros(M); Ts = np.zeros(M)
    if method == 'dt.cycle':
        # State independent: propositions 3.18 and 3.19 end to end, with the
        # index of corollary 3.20(a) corrected (see dpfqn_nc).
        p = P[:, 0]
        q = 1.0 - p
        lG, G, G1 = dpfqn_nc(p, N)
        xtput = G1[N] / G
        for k, ist in enumerate(order):
            tail = 0.0
            for n in range(1, N + 1):
                tail += (q[k] / p[k]) ** n / q[k] * G1[N - n + 1] / G
            Qs[ist] = tail               # E[X_j] = sum_{n>=1} P(X_j>=n)
            Us[ist] = xtput / p[k]       # P(X_j >= 1)
            Ts[ist] = xtput
    else:
        # State dependent: theorem 3.2 through the convolution of dpfqn_ncld.
        lG, G, W, Gc, _ = dpfqn_ncld(P, N)
        for k, ist in enumerate(order):
            marg = W[k] * Gc[k][::-1] / G[N]
            Qs[ist] = float((marg * np.arange(N + 1)).sum())
            Us[ist] = float(1.0 - marg[0])
            Ts[ist] = float((marg[1:] * P[k, :N]).sum())

    # Per-class split by population share; see the docstring for why it is exact.
    Nr = np.asarray(sn.njobs, dtype=float).ravel().copy()
    Nr[~np.isfinite(Nr)] = 0.0
    share = Nr / N if N > 0 else np.zeros(K)

    Q = np.zeros((M, K)); U = np.zeros((M, K))
    R = np.zeros((M, K)); T = np.zeros((M, K))
    C = np.zeros(K); X = np.zeros(K)
    for ist in range(M):
        for r in range(K):
            Q[ist, r] = Qs[ist] * share[r]
            U[ist, r] = Us[ist] * share[r]
            T[ist, r] = Ts[ist] * share[r]
            if T[ist, r] > 0:
                R[ist, r] = Q[ist, r] / T[ist, r]
    refstat = np.asarray(sn.refstat, dtype=int).ravel()
    for r in range(K):
        ref = int(refstat[r]) if r < refstat.size else -1
        if 0 <= ref < M:
            X[r] = T[ref, r]
        if X[r] > 0:
            C[r] = Nr[r] / X[r]
    return Q, U, R, T, C, X, lG
