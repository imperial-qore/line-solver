"""Multi-center branches under Krzesinski (1987) state-dependent routing.

Checks that the product form of eq. (16) extends to a branch holding several
centers, and that the coefficients xi are the branch traffic equations rather
than the Section 3.2 shorthand xi = xi_e, which is exact only when the branch
departure center is visited once. See _kb/16-state-dependent-routing.md
"""

import numpy as np

from line_solver import (ClosedClass, Exp, Network, Queue, SchedStrategy,
                         SolverCTMC)
from line_solver.api.pfqn import pfqn_sdr, pfqn_sdrvisits


def multibranch(pback):
    """Branch 2 = B2a -> B2b, with B2b feeding back to B2a w.p. `pback`."""
    mu = [1.0, 0.9, 0.7, 0.5]   # centers 1, 2a, 2b, 3
    N = 3
    model = Network('sdr_multi')
    node = [Queue(model, 'CPU', SchedStrategy.FCFS),
            Queue(model, 'B2a', SchedStrategy.FCFS),
            Queue(model, 'B2b', SchedStrategy.FCFS),
            Queue(model, 'B3', SchedStrategy.FCFS)]
    jobclass = ClosedClass(model, 'Class1', N, node[0], 0)
    for i in range(4):
        node[i].set_service(jobclass, Exp(mu[i]))
    model.add_link(node[0], node[0])   # denied entry: the busy form of waiting
    model.add_link(node[0], node[1])
    model.add_link(node[0], node[3])
    model.add_link(node[1], node[2])
    model.add_link(node[2], node[0])
    model.add_link(node[3], node[0])
    node[1].set_prob_routing(jobclass, node[2], 1.0)
    if pback > 0:
        model.add_link(node[2], node[1])
        node[2].set_prob_routing(jobclass, node[1], pback)
        node[2].set_prob_routing(jobclass, node[0], 1.0 - pback)
    else:
        node[2].set_prob_routing(jobclass, node[0], 1.0)
    node[3].set_prob_routing(jobclass, node[0], 1.0)

    d = np.zeros((2, 3))
    d[0, 1] = 2.0
    d[0, 2] = 2.0
    d[1, 2] = 2.0
    node[0].set_state_dep_routing(jobclass, node[0],
                                  [[], [node[1], node[2]], [node[3]]],
                                  [0, 1, 2], [-1, -1], d)
    return model, mu, N


def run_case(label, pback):
    model, mu, N = multibranch(pback)
    sn = model.get_struct()

    # The state-INDEPENDENT part of the routing: the complement M-V and the arcs
    # inside a branch. The state-dependent arcs out of the entry center are not
    # part of it.
    P = np.zeros((4, 4))
    P[0, 0] = 1.0                      # complement M-V = {1}
    P[1, 2] = 1.0                      # inside branch 2
    if pback > 0:
        P[2, 1] = pback
    xi = pfqn_sdrvisits(sn.sdr, P)
    S = np.array([[1.0 / m] for m in mu])
    Q, X = pfqn_sdr(S, xi, [N], sn.sdr)[:2]

    ctmc = SolverCTMC(model)
    Qc = np.asarray(ctmc.get_avg_qlen(), dtype=float).reshape(-1, 1)
    Xc = np.asarray(ctmc.get_avg_tput(), dtype=float).reshape(-1, 1)

    print('\n  %s' % label)
    print('    xi          = %s' % np.array2string(xi.ravel(), precision=5))
    print('    pfqn_sdr Q  = %s' % np.array2string(Q.ravel(), precision=6))
    print('    CTMC     Q  = %s' % np.array2string(Qc.ravel(), precision=6))
    print('    pfqn_sdr X  = %s' % np.array2string(X.ravel(), precision=6))
    print('    CTMC     X  = %s' % np.array2string(Xc.ravel(), precision=6))
    print('    max|dQ| = %.3e   max|dX| = %.3e'
          % (np.max(np.abs(Qc - Q)), np.max(np.abs(Xc - X))))

    # the reading asserted verbatim by the paper for E+D: xi = xi_e everywhere
    xiflat = np.ones((4, 1))
    Q2, X2 = pfqn_sdr(S, xiflat, [N], sn.sdr)[:2]
    print('    xi==1 everywhere: max|dQ| = %.3e   max|dX| = %.3e'
          % (np.max(np.abs(Qc - Q2)), np.max(np.abs(Xc - X2))))


if __name__ == '__main__':
    print('\n==== SDR with multi-center branches ====')
    run_case('A: branch 2 = 2a -> 2b (series)', 0.0)
    run_case('B: branch 2 = 2a -> 2b, 2b -> 2a w.p. 0.5 (feedback onto the branch departure)', 0.5)
