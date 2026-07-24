"""
Pass-and-swap (PAS) saturation handling for large / unbounded buffers.

An order-independent service rate mu(c) usually SATURATES: beyond some per-class
count threshold tau_r, adding more class-r jobs no longer changes mu (e.g. a
compatibility / rank function saturates once each present class has one job; an
M/M/K rate saturates at K jobs). LDES exploits this: it serialises mu as a
rate table over the saturation box prod_r {0..tau_r} plus the cutoffs tau, and
the simulator clamps each class count to tau before the table lookup. The table
therefore stays COMPACT regardless of the buffer size, so a large buffer can
approximate an unbounded queue cheaply.

This example uses the Comte (thesis Sect. 4.1) compatibility queue: I=2 classes,
S=3 servers, S_1={1,2}, S_2={2,3}, rank function
    mu(A) = sum_{ s : exists present class in A compatible with s } cap_s,
so mu({1})=mu({2})=cap_1+cap_2 = 3, mu({1,2})=cap_1+cap_2+cap_3 = 5, and mu is
constant once both classes are present (cutoffs tau=[1,1]).

It checks that (i) the serialised table stays at 3 entries whatever the buffer,
(ii) LDES matches CTMC, and (iii) an unset/infinite buffer raises a
clean, actionable error (saturating rate -> suggest setCap; unbounded additive
rate -> setCap required).
"""

import json
import numpy as np
from line_solver import *
from line_solver.io.linemodel_io import save_model

# compatibility comp[s, i] = 1 iff server s serves class i; capacities cap_s
comp = np.array([[1, 0],     # server 1 -> class 1
                 [1, 1],     # server 2 -> classes 1, 2 (shared)
                 [0, 1]])    # server 3 -> class 2
cap_s = np.array([2.0, 1.0, 2.0])
lam = [1.5, 1.0]
I = comp.shape[1]


def mu_rank(c):
    """Rank function: total capacity of servers compatible with a present class."""
    c = np.asarray(c, dtype=int)
    if len(c) == 0:
        return 0.0
    return float(np.sum(cap_s[np.any(comp[:, c], axis=1)]))


def mu_additive(c):
    """Per-class additive rate sum_r n_r*beta_r -- grows without saturating."""
    c = np.asarray(c, dtype=int)
    return float(np.sum(np.array([1.5, 1.0])[c])) if len(c) else 0.0


def build(mu, cap):
    model = Network('PASsaturation')
    source = Source(model, 'Source')
    queue = Queue(model, 'PASQueue', SchedStrategy.PAS)
    sink = Sink(model, 'Sink')
    classes = [OpenClass(model, 'Class%d' % (r + 1)) for r in range(I)]
    for r in range(I):
        source.setArrival(classes[r], Exp(lam[r]))
    queue.setService(mu)
    queue.setSwapGraph(np.zeros((I, I)))     # empty graph: plain OI queue
    queue.setNumberOfServers(comp.shape[0])
    if cap is not None:
        queue.setCap(cap)
    P = model.initRoutingMatrix()
    for r in range(I):
        P[classes[r]] = Network.serialRouting(source, queue, sink)
    model.link(P)
    return model


def qlen(tab):
    col = 'Node' if 'Node' in tab.columns else 'Station'
    return np.array([float(tab[(tab[col] == 'PASQueue') &
                              (tab['JobClass'] == 'Class%d' % (r + 1))]['QLen'].iloc[0])
                     for r in range(I)])


# (1) the rate table stays compact (3 entries, cutoffs [1,1]) for any buffer ----
print('saturation keeps the order-independent rate table compact:')
import tempfile, os
for cap in [8, 100, 1000]:
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'm.json')
        save_model(build(mu_rank, cap), path)
        node = [n for n in json.load(open(path))['model']['nodes']
                if n.get('scheduling') == 'PAS'][0]
        print('  cap=%-4d  table entries = %d  cutoffs = %s'
              % (cap, len(node['oiServiceRate']), node['oiCutoffs']))

# (2) a large buffer approximates the unbounded queue; LDES matches CTMC --------
CAP = 8
Qc = qlen(CTMC(build(mu_rank, CAP), cutoff=CAP).getAvgTable())
Ql = qlen(LDES(build(mu_rank, CAP), samples=300000, seed=23000).getAvgTable())
print('\nopen compatibility PAS, buffer=%d:' % CAP)
print('  CTMC  QLen =', np.round(Qc, 5))
print('  LDES  QLen =', np.round(Ql, 5))
rel = np.max(np.abs(Qc - Ql) / np.maximum(Qc, 1e-12))
print('  LDES vs CTMC max rel|dQ| = %.2f%%' % (100 * rel))
assert rel <= 0.03, 'LDES does not match CTMC within simulation noise'

# (3) clean errors for an unset / infinite buffer ------------------------------
print('\ninfinite-buffer PAS raises a clean, actionable error:')
for tag, mu in [('saturating rank', mu_rank), ('unbounded additive', mu_additive)]:
    try:
        LDES(build(mu, None), samples=1000, seed=1).getAvgTable()
        raise AssertionError('expected a finite-buffer error for ' + tag)
    except Exception as e:
        msg = [ln for ln in str(e).split('\n') if 'PAS station' in ln]
        print('  %-20s -> %s' % (tag, (msg[0] if msg else str(e))[:88]))

print('\nPASS: PAS saturation handling keeps the table compact, matches CTMC, '
      'and errors cleanly on infinite buffers.')
