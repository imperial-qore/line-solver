"""Extract delayed-hit retrieval algorithm inputs from a NetworkStruct (native Python).

Mirror of matlab/src/api/retrieval/cache_retrieval_inputs.m and
java/.../jline/api/retrieval/Cache_retrieval_inputs.java. Single read class (IRM).
Supported station types: IS, PS, SIRO, FCFS, LCFSPR; SIRO/FCFS require exponential
service with identical per-class rates (LCFSPR is exempt).
"""
import numpy as np

from ..sn import NodeType, SchedStrategy
from ..cache.rrm import cache_gamma_lp


def _proc_to_ph(proc_entry, rate):
    """Return (alpha, T) PH representation of a station/class service process.

    Handles the native-Python ``sn.proc`` encodings: an exponential/Erlang/
    HyperExp parameter dict, or a (D0, D1) matrix pair. ``rate`` is the fallback
    mean rate from ``sn.rates`` (used for the exponential case).
    """
    if isinstance(proc_entry, dict):
        if 'probs' in proc_entry and 'rates' in proc_entry:           # HyperExp
            probs = np.asarray(proc_entry['probs'], dtype=float).ravel()
            rates = np.asarray(proc_entry['rates'], dtype=float).ravel()
            return probs.copy(), np.diag(-rates)
        if 'k' in proc_entry:                                          # Erlang
            k = int(proc_entry['k'])
            r = float(proc_entry.get('rate', proc_entry.get('mu', rate)))
            T = np.zeros((k, k))
            for a in range(k):
                T[a, a] = -r
                if a + 1 < k:
                    T[a, a + 1] = r
            al = np.zeros(k)
            al[0] = 1.0
            return al, T
        mu = float(proc_entry.get('rate', proc_entry.get('mu', proc_entry.get('lambda', rate))))
        return np.array([1.0]), np.array([[-mu]])
    if isinstance(proc_entry, (list, tuple)) and len(proc_entry) >= 1:  # (D0, D1) matrices
        D0 = np.asarray(proc_entry[0], dtype=float)
        ph = D0.shape[0]
        if ph == 1:
            return np.array([1.0]), D0.copy()
        # PH initial distribution alpha = -(alpha solves) ; for a PH/MAP without an
        # explicit entry vector, use the embedded one from D1 row sums normalised.
        if len(proc_entry) >= 2:
            D1 = np.asarray(proc_entry[1], dtype=float)
            w = D1.sum(axis=1)
            s = w.sum()
            al = (w / s) if s > 0 else np.concatenate(([1.0], np.zeros(ph - 1)))
        else:
            al = np.concatenate(([1.0], np.zeros(ph - 1)))
        return al, D0.copy()
    # scalar rate fallback
    return np.array([1.0]), np.array([[-float(rate)]])


def _ph_mean(al, T):
    al = np.asarray(al, dtype=float)
    T = np.asarray(T, dtype=float)
    z = np.linalg.solve(-T, np.ones(T.shape[0]))
    return float(al @ z)


def cache_retrieval_inputs(sn):
    """Return a dict with keys m, lambda, gamma, eta, alpha, T, R, station_type,
    jobin_class, queue_node, queue_station, source_rate."""
    K = sn.nclasses

    cache_idx = -1
    for i, nt in enumerate(sn.nodetype):
        if nt == NodeType.CACHE or int(nt) == int(NodeType.CACHE):
            cache_idx = i
            break
    if cache_idx < 0:
        raise RuntimeError("Retrieval analysis requires a Cache node.")
    ch = sn.nodeparam[cache_idx]
    if int(getattr(ch, 'retrieval_system_capacity', 0)) <= 0:
        raise RuntimeError("The Cache node has no retrieval system.")

    m = np.atleast_1d(np.asarray(ch.itemcap, dtype=float)).ravel()
    h = len(m)
    n = int(ch.nitems)

    qidx = ch.retrieval_system_queue_indices
    if len(qidx) != 1:
        raise RuntimeError("Retrieval analysis supports a single read class.")
    jobin_class = next(iter(qidx.keys()))
    queue_node = list(qidx[jobin_class])
    S = len(queue_node)
    queue_station = [int(sn.nodeToStation[q]) for q in queue_node]

    source_node = -1
    for i, nt in enumerate(sn.nodetype):
        if nt == NodeType.SOURCE or int(nt) == int(NodeType.SOURCE):
            source_node = i
            break
    source_ist = int(sn.nodeToStation[source_node])
    source_rate = float(sn.rates[source_ist, jobin_class])
    if np.isnan(source_rate):
        source_rate = 0.0
    pread = np.asarray(ch.pread[jobin_class], dtype=float).ravel()
    lambda_ = source_rate * pread

    # gamma via cache_gamma_lp (single read class)
    lambd3d = np.zeros((1, n, h + 1))
    for i in range(n):
        lambd3d[0, i, :] = lambda_[i]
    accost = getattr(ch, 'accost', None)
    if accost is None:
        R_cost = []
        for _ in range(n):
            Rmat = np.zeros((h + 1, h + 1))
            for j in range(h):
                Rmat[j, j + 1] = 1.0
            Rmat[h, h] = 1.0
            R_cost.append(Rmat)
        R_cost = [R_cost]
    else:
        R_cost = accost
    gamma_res = cache_gamma_lp(lambd3d, R_cost)
    gamma = np.asarray(gamma_res[0], dtype=float).reshape(n, h)

    # station types
    station_type = []
    for s in range(S):
        sd = sn.sched[queue_station[s]]
        sdv = sd.value if hasattr(sd, 'value') else int(sd)
        name = {int(SchedStrategy.INF): "IS", int(SchedStrategy.PS): "PS",
                int(SchedStrategy.SIRO): "SIRO", int(SchedStrategy.FCFS): "FCFS",
                int(SchedStrategy.LCFSPR): "LCFSPR"}.get(int(sdv))
        if name is None:
            raise RuntimeError("Retrieval analysis supports only IS, PS, SIRO, FCFS, LCFSPR stations; node %d" % queue_node[s])
        station_type.append(name)

    # per-item PH service (alpha, T) and routing R
    alpha = [[None] * n for _ in range(S)]
    T = [[None] * n for _ in range(S)]
    R = [np.zeros((S + 1, S + 1)) for _ in range(n)]
    for i in range(n):
        rc = int(ch.retrieval_classes[i, jobin_class])
        for s in range(S):
            ist = queue_station[s]
            al, Tm = _proc_to_ph(sn.proc[ist][rc], sn.rates[ist, rc])
            alpha[s][i] = al
            T[s][i] = Tm
        for s in range(S):
            R[i][0, s + 1] = sn.rtnodes[cache_idx * K + rc, queue_node[s] * K + rc]
            R[i][s + 1, 0] = sn.rtnodes[queue_node[s] * K + rc, cache_idx * K + rc]
            for sp in range(S):
                R[i][s + 1, sp + 1] = sn.rtnodes[queue_node[s] * K + rc, queue_node[sp] * K + rc]

    fsz = [np.asarray(T[s][0]).shape[0] for s in range(S)]
    for s in range(S):
        if station_type[s] in ("SIRO", "FCFS"):
            if fsz[s] > 1:
                raise RuntimeError("SIRO/FCFS retrieval stations require exponential (single-phase) service; node %d" % queue_node[s])
            tau0 = _ph_mean(alpha[s][0], T[s][0])
            for i in range(1, n):
                if abs(_ph_mean(alpha[s][i], T[s][i]) - tau0) > 1e-9 * max(tau0, 1e-300):
                    raise RuntimeError("SIRO/FCFS retrieval stations require identical per-class rates; node %d" % queue_node[s])

    is_is = [station_type[s] == "IS" for s in range(S)]
    ps_idx = [s for s in range(S) if not is_is[s]]
    r = len(ps_idx)
    eta = np.zeros((n, r + 1))
    for i in range(n):
        a = R[i][0, 1:S + 1]
        P = R[i][1:S + 1, 1:S + 1]
        visits = np.linalg.solve(np.eye(S) - P.T, a)
        tau = np.array([_ph_mean(alpha[s][i], T[s][i]) for s in range(S)])
        for s in range(S):
            if is_is[s]:
                eta[i, 0] += visits[s] * tau[s]
        for p in range(r):
            eta[i, 1 + p] = visits[ps_idx[p]] * tau[ps_idx[p]]

    return {
        'm': m, 'lambda': lambda_, 'gamma': gamma, 'eta': eta,
        'alpha': alpha, 'T': T, 'R': R, 'station_type': station_type,
        'jobin_class': jobin_class, 'queue_node': queue_node,
        'queue_station': queue_station, 'source_rate': source_rate,
        'cache_idx': cache_idx,
    }
