"""
PAS_NC_SAMPLING  Importance-sampling normalizing-constant analysis of a closed
pass-and-swap (P&S) tandem via SolverNC method='sampling' (Casale, Comte &
Dorsman, 2026). The model is the Figure 6 closed tandem of Comte & Dorsman
(2021, arXiv:2009.12299) also used in PAS_CLOSED_TANDEM_FIG6.

A non-empty swap graph makes the ordered-state chain reducible; the recurrent
communicating class carries the per-class product form pi(c)=Phi_1 Phi_2/G_C.
SolverNC importance sampling estimates G_C and the mean queue lengths by
auto-normalized importance sampling (API PFQN_PAS_IS), scaling to populations
where the exact ordered-state CTMC is intractable. Here we validate it against
the exact CTMC.
"""

import numpy as np

from line_solver import (ClosedClass, CTMC, GlobalConstants, NC, Network,
                         Queue, SchedStrategy, VerboseLevel)


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

    MU1, MU2 = 1.0, 1.3
    EDGES = [(1, 3), (1, 4), (2, 4), (2, 5), (3, 6), (4, 6), (5, 6)]
    G = np.zeros((6, 6))
    for a, b in EDGES:
        G[a - 1, b - 1] = 1
        G[b - 1, a - 1] = 1

    model = Network('PASsampling')
    q1 = Queue(model, 'PASQueue1', SchedStrategy.PAS)
    q2 = Queue(model, 'PASQueue2', SchedStrategy.PAS)
    jobclass = [ClosedClass(model, 'Class%d' % (r + 1), 1, q1) for r in range(6)]
    q1.setService(lambda c: MU1)       # head-only single server
    q2.setService(lambda c: MU2)
    q1.setSwapGraph(G); q1.setNumberOfServers(1); q1.setCap(6)
    q2.setSwapGraph(G); q2.setNumberOfServers(1); q2.setCap(6)
    # MATLAB writes P{class}(q1,q2)=1 into the routing matrix and links THAT;
    # setProbRouting alone leaves the matrix empty here, and the chain then has
    # no arcs at all.
    P = model.initRoutingMatrix()
    for r in range(6):
        P.set(jobclass[r], jobclass[r], q1, q2, 1.0)
        P.set(jobclass[r], jobclass[r], q2, q1, 1.0)
    model.link(P)
    q1.setState(np.array([1, 2, 3, 4, 5, 6]))   # Fig. 6a initial placement

    print('=== CTMC (exact) ===')
    Tc = CTMC(model, 'exact', cutoff=6).getAvgTable()
    print("=== SolverNC method='sampling' (importance sampling, 5e5 samples) ===")
    Tn = NC(model, method='sampling', samples=int(5e5), seed=777,
            verbose=False).getAvgTable()

    qc = np.asarray(Tc['QLen'], dtype=float)
    qn = np.asarray(Tn['QLen'], dtype=float)
    print('\n  station     class    CTMC      NC-samp')
    for i in range(len(qc)):
        print('  %-9s  %-6s %9.5f %9.5f'
              % (Tc['Station'][i], Tc['JobClass'][i], qc[i], qn[i]))

    err = float(np.max(np.abs(qn - qc)))
    print('\nNC-samp vs CTMC:  max|dQ| = %.3e (importance-sampling noise)' % err)
    assert err <= 5e-2, 'NC sampling does not match the exact CTMC within IS tolerance'
    print("PASS: SolverNC 'sampling' matches the exact CTMC within importance-sampling noise.")
