"""Open Whittle network: a bandwidth-sharing model.

One route holds SEVERAL links at once, which no per-station rate scaling can
express. This is the 2-link linear network: route 1 crosses both links, routes 2
and 3 use one link each.

    A = [1 1 0;      link 1 carries routes 1 and 2
         1 0 1]      link 2 carries routes 1 and 3

Capacity is shared by BALANCED FAIRNESS, whose rates come from the recursion
    Phi(n) = max_l (1/C_l) sum_{s in l} Phi(n-e_s),  x_s(n) = Phi(n-e_s)/Phi(n).
That construction satisfies the Whittle balance property by design, so the
stationary law is pi(n) ~ Phi(n) prod rho_s**n_s and is insensitive.

Each route is modelled as its own PS queue fed by its own open class, so the
queue populations ARE the coordinates of the Whittle state n.

MATLAB twin: matlab/examples/advanced/loadDependent/ld_whittle_bandwidth.m
"""
import itertools

import numpy as np

from line_solver import (Disabled, Exp, Network, OpenClass, Queue,
                         SchedStrategy, Sink, SolverCTMC, Source)

A = np.array([[1, 1, 0], [1, 0, 1]])
C = np.array([1.0, 1.0])
nu = np.array([0.20, 0.30, 0.30])  # flow arrival rates
mu = np.array([1.00, 1.00, 1.00])  # 1/mean file size
cutoff = 6                         # per-route truncation of the open state space


def balance_function(A, C, cut):
    """Phi(n) = max_l (1/C_l) sum_{s in l} Phi(n-e_s), over the bounded lattice."""
    S = A.shape[1]
    base = np.asarray(cut, dtype=int) + 1
    states = np.array(list(itertools.product(*[range(b) for b in base])), dtype=int)
    # itertools.product varies the LAST coordinate fastest; index accordingly
    index = {tuple(n): k for k, n in enumerate(states)}
    Phi = np.zeros(len(states))
    Phi[index[tuple(np.zeros(S, dtype=int))]] = 1.0
    for k in np.argsort(states.sum(axis=1), kind='stable'):
        n = states[k]
        if not n.any():
            continue
        best = 0.0
        for l in range(A.shape[0]):
            acc = 0.0
            for s in range(S):
                if A[l, s] > 0 and n[s] > 0:
                    m = n.copy()
                    m[s] -= 1
                    acc += Phi[index[tuple(m)]]
            best = max(best, acc / C[l])
        Phi[k] = best
    return Phi, index


def bf_rates(n, Phi, index):
    S = len(n)
    x = np.zeros(S)
    key = tuple(int(v) for v in n)
    if not any(key):
        return x
    den = Phi[index[key]]
    for s in range(S):
        if n[s] > 0:
            m = np.asarray(n, dtype=int).copy()
            m[s] -= 1
            x[s] = Phi[index[tuple(m)]] / den
    return x


def bf_means(Phi, index, rho):
    states = np.array(list(index.keys()), dtype=int)
    w = np.array([Phi[index[tuple(n)]] * np.prod(rho ** n) for n in states])
    w = w / w.sum()
    return w @ states


Phi, index = balance_function(A, C, [cutoff] * 3)

model = Network('model')
source = Source(model, 'Source')
routes = [Queue(model, 'Route%d' % (s + 1), SchedStrategy.PS) for s in range(3)]
sink = Sink(model, 'Sink')
classes = [OpenClass(model, 'Route%dFlows' % (s + 1)) for s in range(3)]
for s in range(3):
    source.setArrival(classes[s], Exp(nu[s]))
    for t in range(3):
        if s == t:
            routes[t].setService(classes[s], Exp(mu[s]))
        else:
            routes[t].setService(classes[s], Disabled())
P = model.initRoutingMatrix()
for s in range(3):
    P.set(classes[s], classes[s], source, routes[s], 1.0)
    P.set(classes[s], classes[s], routes[s], sink, 1.0)
model.link(P)

sn = model.getStruct()
# getNodeIndex is 1-based while nodeToStation is a 0-based array
idx = [int(sn.nodeToStation[model.getNodeIndex('Route%d' % (s + 1)) - 1]) for s in range(3)]


def phi(n):
    npop = np.array([n[idx[s], s] for s in range(3)])
    v = np.ones(n.shape)
    x = bf_rates(npop, Phi, index)
    for s in range(3):
        v[idx[s], s] = x[s]
    return v


# The third argument is the per-slot open-class truncation used when phi is
# materialized onto the JSON wire; it matches the cutoff below, which is also the
# range over which the balance function Phi was built.
model.set_global_dependence(phi, 1.0, cutoff)

print(SolverCTMC(model, cutoff=cutoff).getAvgTable())

# Cross-check against the closed-form product form. Truncating a reversible chain
# preserves the conditional law, so this agrees to solver precision rather than to
# a truncation-limited tolerance.
En = bf_means(Phi, index, nu / mu)
print('\nproduct form E[n] = [%.6f %.6f %.6f]' % tuple(En))
