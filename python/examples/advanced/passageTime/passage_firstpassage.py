"""
First passage times in a Markov chain, and the exact cycle time along an
overtake-free path of a closed tree-like product-form network.

Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions in
Large Markov Chains", 2002.

Twin of matlab/examples/advanced/passageTime/passage_firstpassage.m.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np

from line_solver.api.mc import (ctmc_hitting_time, ctmc_passage_moments,
                                ctmc_passage_time, smp_passage_moments)
from line_solver.api.mc.ctmc import ctmc_makeinfgen
from line_solver.api.pfqn import pfqn_cyclet_ofree


def main():
    # 1. The time for an M/M/1/K queue to fill from empty.
    # This is a first passage into a STATE SET, which no response-time getter
    # can express: it is the chain reaching a marking, not a job finishing.
    K, lam, mu = 6, 1.0, 1.5
    n = K + 1
    Q = np.zeros((n, n))
    for i in range(n):
        if i + 1 < n:
            Q[i, i + 1] = lam
        if i > 0:
            Q[i, i - 1] = mu
    Q = ctmc_makeinfgen(Q)

    pi0 = np.zeros(n)
    pi0[0] = 1.0            # start empty
    target = [n - 1]        # the full buffer

    mall, m = ctmc_passage_moments(Q, pi0, target, 3)
    var = m[1] - m[0] ** 2
    print("Time to fill an M/M/1/%d from empty (lambda=%g, mu=%g)" % (K, lam, mu))
    print("  mean            = %.6f" % m[0])
    print("  variance        = %.6f" % var)
    print("  coeff. of var.  = %.6f" % (np.sqrt(var) / m[0]))
    print("  from state K-1  = %.6f (one arrival away)" % mall[n - 2, 0])

    tset = np.linspace(0.0, 4 * m[0], 400)
    F, f, _ = ctmc_passage_time(Q, pi0, target, tset)
    print("  P(fill <= mean) = %.6f" % np.interp(m[0], tset, F))

    h = ctmc_hitting_time(Q, target)
    print("  hitting times   = %s" % np.round(h, 3))

    # 2. The same passage with a non-exponential sojourn (semi-Markov).
    # The embedded chain is unchanged; only the holding-time law moves. Nothing
    # in a generator can express this, which is why the semi-Markov route exists.
    rate = -np.diag(Q)
    P = np.zeros((n, n))
    for i in range(n):
        if rate[i] > 0:
            P[i, :] = Q[i, :] / rate[i]
            P[i, i] = 0.0
        else:
            P[i, i] = 1.0
    hmom = np.zeros((n, 3))
    for i in range(n):
        d = 1.0 / rate[i]
        hmom[i, :] = [d, d ** 2, d ** 3]
    _, mD = smp_passage_moments(P, hmom, pi0, target, 3)
    varD = mD[1] - mD[0] ** 2
    print("\nSame embedded chain, DETERMINISTIC sojourns of equal mean:")
    print("  mean            = %.6f (unchanged, as it must be)" % mD[0])
    print("  coeff. of var.  = %.6f (was %.6f)"
          % (np.sqrt(varD) / mD[0], np.sqrt(var) / m[0]))

    # 3. The cycle time of the tree-like network of Fig. 6 of the paper.
    mu6 = np.array([3.0, 5.0, 4.0, 6.0, 2.0, 1.0])
    p12, p13, p14 = 0.2, 0.5, 0.3
    v = np.array([1.0, p12, p13, p14, p12, p14])
    N = 18
    paths = [[0, 2], [0, 1, 4], [0, 3, 5]]
    tt = np.arange(0.0, 40.001, 0.25)
    fc, Fc, mom, out = pfqn_cyclet_ofree(v, mu6, N, paths, tt,
                                         pathprob=[p13, p12, p14])
    print("\nTree network of Fig. 6, N = %d customers" % N)
    print("  moments  = %.5f  %.4f  %.3f" % (mom[0], mom[1], mom[2]))
    print("  paper    = 6.12717  53.3067  612.887")
    print("  routes   = %s" % ", ".join(o['method'] for o in out))
    print("  P(cycle <= mean) = %.6f" % np.interp(mom[0], tt, Fc))


if __name__ == '__main__':
    main()
