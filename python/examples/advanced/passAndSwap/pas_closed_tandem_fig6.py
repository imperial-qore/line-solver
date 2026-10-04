"""
Closed tandem of two pass-and-swap (PAS) queues, reproducing Figures 5 and 6
of Comte & Dorsman, "Pass-and-Swap Queues" (2021, arXiv:2009.12299).

    Topology (Fig. 6):   --> PASQueue1 --> PASQueue2 -->   (closed tandem)

Six classes, one customer each; both queues share the Figure 5 swapping graph
(edges 1-3, 1-4, 2-4, 2-5, 3-6, 4-6, 5-6). Head-only single-server service;
the departing customer is chosen by the pass-and-swap scan and routed to the
back of the other queue.

A closed PAS network with a non-empty swapping graph is REDUCIBLE (the
pass-and-swap mechanism conserves the placement order), so the initial job
placement is a REQUIRED model input that selects the recurrent component on
which the stationary distribution is supported. We use the placement of
Fig. 6a: queue 1 = (1,2,3,4,5,6), queue 2 empty. The PAS state list is 1-based
(jobClass+1, 0 = empty), so setState uses 1..6.

The CTMC result is validated against an independent product-form reference that
builds the Markov chain directly from the pass-and-swap algorithm.
"""

import numpy as np
from line_solver import *

MU1, MU2 = 1.0, 1.3
EDGES = [(1, 3), (1, 4), (2, 4), (2, 5), (3, 6), (4, 6), (5, 6)]
G = np.zeros((6, 6))
for a, b in EDGES:
    G[a - 1, b - 1] = 1
    G[b - 1, a - 1] = 1


# ---- LINE model -------------------------------------------------------
model = Network('PASclosedTandem')
q1 = Queue(model, 'PASQueue1', SchedStrategy.PAS)
q2 = Queue(model, 'PASQueue2', SchedStrategy.PAS)
jobclass = [ClosedClass(model, 'Class%d' % (r + 1), 1, q1) for r in range(6)]
q1.setService(lambda c: MU1)        # head-only: mu(prefix)=mu1
q2.setService(lambda c: MU2)
q1.setSwapGraph(G); q1.setNumberOfServers(1); q1.setCap(6)
q2.setSwapGraph(G); q2.setNumberOfServers(1); q2.setCap(6)
model.addLink(q1, q2)
model.addLink(q2, q1)
for r in range(6):
    q1.setProbRouting(jobclass[r], q2, 1.0)
    q2.setProbRouting(jobclass[r], q1, 1.0)
model.link(model.initRoutingMatrix())

# Required initial placement (Fig. 6a), 1-based ordered class list.
q1.setState(np.array([1, 2, 3, 4, 5, 6]))

T = CTMC(model, cutoff=6).getAvgTable()
print(T)
col = 'Node' if 'Node' in T.columns else 'Station'


# ---- independent product-form reference -------------------------------
def psalgorithm(lst, p):
    """Pass-and-swap algorithm: customer at position p completes; classes shift one
    step along the swap chain and the served slot is removed; returns the new
    list and the departing class."""
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


def reference():
    init = (tuple(range(1, 7)), ())
    states = [init]
    index = {init: 0}
    edges = []
    fr = 0
    while fr < len(states):
        l1, l2 = states[fr]
        if l1:                                   # head of queue 1 completes
            nl, dep = psalgorithm(list(l1), 0)
            t = (tuple(nl), l2 + (dep,))
            j = index.setdefault(t, len(states))
            if j == len(states):
                states.append(t)
            edges.append((fr, j, MU1))
        if l2:                                   # head of queue 2 completes
            nl, dep = psalgorithm(list(l2), 0)
            t = (l1 + (dep,), tuple(nl))
            j = index.setdefault(t, len(states))
            if j == len(states):
                states.append(t)
            edges.append((fr, j, MU2))
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
    Q1 = np.zeros(6)
    Q2 = np.zeros(6)
    for s, (l1, l2) in enumerate(states):
        for r in range(6):
            Q1[r] += pi[s] * l1.count(r + 1)
            Q2[r] += pi[s] * l2.count(r + 1)
    return Q1, Q2


Qref1, Qref2 = reference()
Qref = np.concatenate([Qref1, Qref2])
qc = np.array([float(T[(T[col] == nm) & (T['JobClass'] == 'Class%d' % (r + 1))]['QLen'].iloc[0])
               for nm in ['PASQueue1', 'PASQueue2'] for r in range(6)])

print('\n  station    class   reference     CTMC')
names = ['PASQueue1'] * 6 + ['PASQueue2'] * 6
for i in range(12):
    print('  %-9s  Class%-2d %9.5f %9.5f' % (names[i], (i % 6) + 1, Qref[i], qc[i]))
err = np.max(np.abs(qc - Qref))
print('\nCTMC vs reference: max|dQ| = %.3e' % err)
assert err <= 1e-6, 'CTMC does not match the pass-and-swap reference'
print('PASS: CTMC matches the pass-and-swap reference.')
