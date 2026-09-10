"""Globally state-dependent (Whittle) model.

set_global_dependence declares a rate scaling phi(n) over the FULL
(nstations, nclasses) population matrix, not just the population local to one
station. Here two PS stations share one unit of capacity, phi_s(n) = n_s/|n|,
which is the single-link allocation every alpha-fair rule collapses to.

phi satisfies the Whittle balance property

    phi_s(n) phi_t(n-e_s) = phi_t(n) phi_s(n-e_t),

so the chain is reversible, has the product form pi(n) ~ Phi(n) prod rho**n and
is INSENSITIVE: the means below do not change when the exponential service is
replaced by an Erlang or a hyperexponential of the same mean.

MATLAB twin: matlab/examples/advanced/loadDependent/ld_global_dependence.m
"""
import numpy as np

from line_solver import (ClosedClass, Exp, Network, Queue, SchedStrategy,
                         SolverCTMC, SolverSSA)
from line_solver.api.sn import sn_gd_balance

N = 3  # number of jobs

model = Network('model')
queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
class1 = ClosedClass(model, 'Class1', N, queue1, 0)
queue1.setService(class1, Exp.fitMean(1.0))
queue2.setService(class1, Exp.fitMean(0.5))
model.link(Network.serialRouting(queue1, queue2))

sn = model.getStruct()
# getNodeIndex is 1-based while nodeToStation is a 0-based array
i1 = int(sn.nodeToStation[model.getNodeIndex('Queue1') - 1])
i2 = int(sn.nodeToStation[model.getNodeIndex('Queue2') - 1])


def share(n):
    """phi as SolverCTMC sees it: an (nstations, nclasses) matrix of scalings."""
    v = np.ones(n.shape)
    tot = n[i1, :].sum() + n[i2, :].sum()
    if tot > 0:
        v[i1, :] = n[i1, :].sum() / tot
        v[i2, :] = n[i2, :].sum() / tot
    return v


# peak scaling is 1: no station ever receives more than the whole link
model.set_global_dependence(share, 1.0)

# Only SolverCTMC and SolverSSA plumb the handle; every other solver rejects the
# model rather than silently solving it unscaled (feature 'GlobalDependence').
print(SolverCTMC(model, 'exact').getAvgTable())

# SolverSSA carries the SAME factorization on the sample path: phi(n) is a
# constant within a state, so it is evaluated once per state and multiplies every
# station service rate there. Its NRM engine cannot (its propensity closures see
# one station's population slice), so the model is routed to the serial engine.
print(SolverSSA(model, seed=23000, samples=200000).getAvgTable())

# The balance property is checkable, and is what separates a Whittle network from
# an arbitrary state-dependent rate.
viol, _ = sn_gd_balance(lambda n: np.asarray(n) / max(np.sum(n), 1e-300), [4, 4])
print('\nworst relative balance violation: %.2e (balanced when ~0)' % viol)
