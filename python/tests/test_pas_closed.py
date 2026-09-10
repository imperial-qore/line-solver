"""
Native-Python CTMC tests for CLOSED pass-and-swap (PAS) networks.

These cover the REDUCIBLE case, which ``test_pas.py`` does not: its models are
all open (Source -> PASQueue -> Sink) and hence irreducible. A closed PAS
network with a non-empty swapping graph conserves the placement order (Comte &
Dorsman, "Pass-and-Swap Queues", 2021, arXiv:2009.12299, Prop. 2), so its
generator has one bottom SCC per reachable placement order and the stationary
distribution is component-dependent. The declared initial placement selects the
component, which is what the analyzer's seeded block decomposition must honour.

Both tests check against an INDEPENDENT analytical oracle rather than a recorded
golden: the pass-and-swap Markov chain built directly from the algorithm, and
the brute-force product-form normalizing constant. A recorded golden cannot
catch a consistent error; these can.
"""

import numpy as np
import pytest

from line_solver import Network, Queue, ClosedClass, SchedStrategy, CTMC, circul

TOL = 1e-6

# Figure 5 swapping graph of Comte & Dorsman: edges 1-3, 1-4, 2-4, 2-5, 3-6, 4-6, 5-6
_EDGES = [(1, 3), (1, 4), (2, 4), (2, 5), (3, 6), (4, 6), (5, 6)]
_G6 = np.zeros((6, 6))
for _a, _b in _EDGES:
    _G6[_a - 1, _b - 1] = 1
    _G6[_b - 1, _a - 1] = 1

MU1, MU2 = 1.0, 1.3


def _qlen(model, cutoff, nstations, nclasses):
    """Mean queue lengths as a flat station-major vector, read off the array API.

    getAvgTable would work too now that IndexedTable no longer swallows boolean
    masks (see test_indexed_table.py), but the array API needs no name matching
    at all, so it cannot be perturbed by table formatting.
    """
    q = np.asarray(CTMC(model, cutoff=cutoff).getAvgQLen(), dtype=float)
    q = q.reshape(nstations, nclasses)
    return q.ravel()


def _tput(model, cutoff, nstations, nclasses):
    t = np.asarray(CTMC(model, cutoff=cutoff).getAvgTput(), dtype=float)
    return t.reshape(nstations, nclasses)


# ---------------------------------------------------------------------------
# 1. Closed tandem, Fig. 6a placement, against the pass-and-swap chain itself
# ---------------------------------------------------------------------------

def _ps_algorithm(lst, p, G):
    """One pass-and-swap completion at position p: classes shift one step along
    the swap chain, the served slot is removed. Returns (new list, departing class)."""
    chain = [p]
    cur = p
    while True:
        nxt = -1
        for j in range(cur + 1, len(lst)):
            if G[lst[cur] - 1, lst[j] - 1]:
                nxt = j
                break
        if nxt < 0:
            break
        chain.append(nxt)
        cur = nxt
    dep = lst[chain[-1]]
    new = list(lst)
    for i in range(len(chain) - 1):
        new[chain[i + 1]] = lst[chain[i]]
    del new[chain[0]]
    return new, dep


def _ps_reference(G, mu1, mu2, nclasses=6):
    """Independent reference: enumerate the pass-and-swap chain over ordered
    placements from the Fig. 6a start and solve it directly."""
    init = (tuple(range(1, nclasses + 1)), ())
    states = [init]
    index = {init: 0}
    edges = []
    fr = 0
    while fr < len(states):
        l1, l2 = states[fr]
        if l1:
            nl, dep = _ps_algorithm(list(l1), 0, G)
            t = (tuple(nl), l2 + (dep,))
            j = index.setdefault(t, len(states))
            if j == len(states):
                states.append(t)
            edges.append((fr, j, mu1))
        if l2:
            nl, dep = _ps_algorithm(list(l2), 0, G)
            t = (l1 + (dep,), tuple(nl))
            j = index.setdefault(t, len(states))
            if j == len(states):
                states.append(t)
            edges.append((fr, j, mu2))
        fr += 1

    ns = len(states)
    Qm = np.zeros((ns, ns))
    for i, j, r in edges:
        Qm[i, j] += r
    for i in range(ns):
        Qm[i, i] = -Qm[i].sum()
    A = np.vstack([Qm.T, np.ones(ns)])
    b = np.append(np.zeros(ns), 1.0)
    pi = np.linalg.lstsq(A, b, rcond=None)[0]

    Q1 = np.zeros(nclasses)
    Q2 = np.zeros(nclasses)
    for s, (l1, l2) in enumerate(states):
        for r in range(nclasses):
            Q1[r] += pi[s] * l1.count(r + 1)
            Q2[r] += pi[s] * l2.count(r + 1)
    return np.concatenate([Q1, Q2])


def _build_closed_tandem(placement):
    model = Network('PASclosedTandem')
    q1 = Queue(model, 'PASQueue1', SchedStrategy.PAS)
    q2 = Queue(model, 'PASQueue2', SchedStrategy.PAS)
    jobclass = [ClosedClass(model, 'Class%d' % (r + 1), 1, q1) for r in range(6)]
    q1.setService(lambda c: MU1)        # head-only service: mu(prefix) = mu1
    q2.setService(lambda c: MU2)
    q1.setSwapGraph(_G6); q1.setNumberOfServers(1); q1.setCap(6)
    q2.setSwapGraph(_G6); q2.setNumberOfServers(1); q2.setCap(6)
    model.addLink(q1, q2)
    model.addLink(q2, q1)
    for r in range(6):
        q1.setProbRouting(jobclass[r], q2, 1.0)
        q2.setProbRouting(jobclass[r], q1, 1.0)
    model.link(model.initRoutingMatrix())
    if placement is not None:
        # PAS state list is 1-based (jobClass+1, 0 = empty)
        q1.setState(np.array(placement))
    return model


def test_pas_closed_tandem_matches_ps_chain():
    """The seeded component must be the one the declared placement selects."""
    model = _build_closed_tandem([1, 2, 3, 4, 5, 6])
    got = _qlen(model, 6, 2, 6)
    ref = _ps_reference(_G6, MU1, MU2)
    assert np.allclose(got, ref, atol=TOL), \
        f'closed PAS QLen does not match the pass-and-swap chain:\n got={got}\n ref={ref}'


def test_pas_closed_tandem_conserves_population():
    """One customer per class, so every class sums to 1 across the two stations."""
    model = _build_closed_tandem([1, 2, 3, 4, 5, 6])
    q = _qlen(model, 6, 2, 6)
    per_class = q[:6] + q[6:]
    assert np.allclose(per_class, np.ones(6), atol=TOL), \
        f'closed PAS population not conserved per class: {per_class}'


def test_pas_closed_tandem_requires_placement():
    """Several recurrent components and no declared placement is ill-posed: the
    solver must refuse rather than average the mirror components."""
    model = _build_closed_tandem(None)
    with pytest.raises(RuntimeError, match='initial job placement'):
        CTMC(model, cutoff=6).getAvgQLen()


# ---------------------------------------------------------------------------
# 2. Closed cyclic OI network against the brute-force product form
# ---------------------------------------------------------------------------

def _oi_sum(counts, mu):
    """sum over ordered placements of prod 1/mu(prefix), by depth-first descent."""
    def dfs(cnt, prefix, w):
        if not any(cnt):
            return w
        s = 0.0
        for r in range(len(cnt)):
            if cnt[r] > 0:
                c2 = list(cnt)
                c2[r] -= 1
                p2 = prefix + [r + 1]
                s += dfs(c2, p2, w / mu(p2))
        return s
    return dfs(list(counts), [], 1.0)


def _enum_splits(K):
    out = [[]]
    for k in K:
        out = [row + [m] for row in out for m in range(k + 1)]
    return out


def _bruteforce(K, mu1, mu2):
    """Exact product-form G, per-class throughputs X_r = G(K-1_r)/G(K), and
    mean queue lengths, by full ordered-placement enumeration."""
    def norm_const(Kv):
        return sum(_oi_sum(m, mu1) * _oi_sum([Kv[r] - m[r] for r in range(len(Kv))], mu2)
                   for m in _enum_splits(Kv))

    G = norm_const(K)
    R = len(K)
    X = np.zeros(R)
    for r in range(R):
        Km = list(K)
        Km[r] -= 1
        X[r] = norm_const(Km) / G

    Q1 = np.zeros(R)
    for m in _enum_splits(K):
        w = _oi_sum(m, mu1) * _oi_sum([K[r] - m[r] for r in range(R)], mu2)
        for r in range(R):
            Q1[r] += m[r] * w
    Q1 /= G
    Q2 = np.array(K, dtype=float) - Q1
    return X, np.concatenate([Q1, Q2])


def test_pas_closed_cyclic_matches_bruteforce_product_form():
    """Two order-independent PAS stations in a cycle are product-form; compare
    against the brute-force normalizing constant."""
    K = [2, 2]
    s1, k1 = 1.0, 2
    beta2 = [1.5, 1.0]
    R = len(K)

    # Total service rate as a function of the ordered class list (1-based ids)
    def mu1(c):
        return min(len(c), k1) * s1

    def mu2(c):
        return float(sum(beta2[int(j) - 1] for j in c))

    Xref, Qref = _bruteforce(K, mu1, mu2)

    model = Network('PAScyclic')
    q1 = Queue(model, 'PASQueue1', SchedStrategy.PAS)
    q2 = Queue(model, 'PASQueue2', SchedStrategy.PAS)
    jobclass = [ClosedClass(model, 'Class%d' % (r + 1), K[r], q1) for r in range(R)]
    # NOTE the index-base asymmetry: setState takes a 1-BASED ordered class list
    # (0 = empty) but the service function receives the list 0-BASED, so the ids
    # are shifted back up before indexing beta2.
    q1.setService(lambda c: mu1([int(j) + 1 for j in np.atleast_1d(c)]))
    q2.setService(lambda c: mu2([int(j) + 1 for j in np.atleast_1d(c)]))
    # Empty swap graph: a plain order-independent queue. The result is
    # swap-graph invariant, so the product form applies.
    q1.setSwapGraph(np.zeros((R, R))); q1.setNumberOfServers(k1); q1.setCap(sum(K))
    q2.setSwapGraph(np.zeros((R, R))); q2.setNumberOfServers(sum(K)); q2.setCap(sum(K))
    model.addLink(q1, q2)
    model.addLink(q2, q1)
    for r in range(R):
        q1.setProbRouting(jobclass[r], q2, 1.0)
        q2.setProbRouting(jobclass[r], q1, 1.0)
    model.link(model.initRoutingMatrix())

    Qgot = _qlen(model, sum(K), 2, R)
    Xgot = _tput(model, sum(K), 2, R)[0]   # one visit per station in a cycle

    assert np.allclose(Qgot, Qref, atol=1e-9), \
        f'closed cyclic PAS QLen vs brute force:\n got={Qgot}\n ref={Qref}'
    assert np.allclose(Xgot, Xref, atol=1e-9), \
        f'closed cyclic PAS Tput vs brute force:\n got={Xgot}\n ref={Xref}'
