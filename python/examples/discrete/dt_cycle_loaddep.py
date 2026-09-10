"""
Closed cycle of state dependent Bernoulli servers.

Theorem 3.2 of Daduna (2001) keeps the product form when the service
probability of station j depends on its own queue length, with node weight

    w_j(n) = prod_{h=1}^{n-1} q_j(h) / prod_{h=1}^{n} p_j(h).

Note where the state dependence sits: the missing q_j in the numerator is tied
to the actual queue length, not to the node being non-empty, so the tidy
(1/q_j)^{1{n>0}} of the state independent case cannot be factored out.

Station 2 below runs at p_2(n) = p_2 min(n,2), the discrete-time analogue of
adding a second server. That is admissible only because the dependence is
expressed as a state dependent single server: a genuine multiserver node inside
a cycle of geometrical queues destroys the product form for every finite server
count (Pestien and Ramakrishnan, cited before example 2.10), and SolverNC
rejects one rather than approximating it.
"""
import numpy as np

from line_solver import *
from line_solver.api.dpfqn import dpfqn_ncld

if __name__ == "__main__":
    p = [0.5, 0.25, 0.7]
    N = 5

    model = Network('BernoulliCycleLD')

    station = [Queue(model, 'Queue%d' % (j + 1), SchedStrategy.FCFS)
               for j in range(len(p))]

    job_class = ClosedClass(model, 'Jobs', N, station[0])
    for j in range(len(p)):
        station[j].setService(job_class, Geometric(p[j]))
    station[1].setLoadDependence(np.minimum(np.arange(1, N + 1), 2).astype(float))

    model.link(Network.serialRouting(*station))

    solver = SolverNC(model)
    solver.options.config = {'slotted': True}
    print(solver.getAvgTable())

    # Marginal queue length law of station 2 from the convolution constants.
    P = np.tile(np.array(p).reshape(-1, 1), (1, N))
    P[1, :] = p[1] * np.minimum(np.arange(1, N + 1), 2)
    lG, G, W, Gc, Wa = dpfqn_ncld(P, N)
    marg = W[1] * Gc[1][::-1] / G[N]
    print('P(X_2 = 0..%d) = %s' % (N, np.round(marg, 6)))
