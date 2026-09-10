"""
Large-scale transient fluid analysis with trajectory-based iteration (TBI)

This example builds a vehicle-sharing-style closed network with M queueing
stations and one Erlang transit delay per ordered station pair, giving O(M^2)
stations and, with 16 Erlang phases per delay, a fluid ODE with about 1500 state
variables. At this size the monolithic stiff fluid solution (method 'closing')
takes several minutes, dominated by the Jacobian factorizations, while
trajectory-based iteration (TBI) solves each station cell separately against
frozen inbound trajectories and completes in seconds. See:

  M. Sheldon, D. Tuncer, G. Casale, "TBI: Transient Hierarchical Modeling of
  Large-Scale Vehicle Sharing Systems", IEEE Transactions on Intelligent
  Transportation Systems.

Each cell holds one queueing station together with its outbound transit delays,
mirroring the spatial submodels of the paper.
"""

import time

import numpy as np

from line_solver import (ClosedClass, Delay, Erlang, FLD, GlobalConstants,
                         Network, Queue, SchedStrategy, VerboseLevel)

M = 10     # queueing stations; the model has M + M*(M-1) stations in total
N = 250    # closed population (vehicles)
KPH = 16   # Erlang phases per transit delay
TEND = 20  # transient horizon


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

    model = Network('tbi_largescale')
    rng = np.random.default_rng(1)

    Q = [Queue(model, 'Q%d' % (i + 1), SchedStrategy.PS) for i in range(M)]
    D = {}
    for i in range(M):
        for j in range(M):
            if i != j:
                D[(i, j)] = Delay(model, 'D%d_%d' % (i + 1, j + 1))

    job = ClosedClass(model, 'C1', N, Q[0])
    lam = np.zeros((M, M))
    for i in range(M):
        Q[i].setService(job, Erlang.fitMeanAndOrder(1.0 / (1 + 3 * rng.random()), 4))
        for j in range(M):
            if i != j:
                D[(i, j)].setService(
                    job, Erlang.fitMeanAndOrder(0.2 + 2 * rng.random(), KPH))
                lam[i, j] = rng.random()

    P = model.initRoutingMatrix()
    for i in range(M):
        pr = lam[i, :] / lam[i, :].sum()
        for j in range(M):
            if i != j:
                P.set(job, job, Q[i], D[(i, j)], float(pr[j]))
                P.set(job, job, D[(i, j)], Q[j], 1.0)
    model.link(P)
    model.initDefault()

    # one cell per queueing station plus its outbound transit delays
    names = [s.getName() for s in model.getStations()]
    cells = []
    for i in range(M):
        idx = [names.index('Q%d' % (i + 1))]
        for j in range(M):
            if i != j:
                idx.append(names.index('D%d_%d' % (i + 1, j + 1)))
        cells.append(idx)

    solver = FLD(model, method='tbi', timespan=[0, TEND], stiff=True)
    solver.options.config['tbi_cells'] = cells
    t0 = time.time()
    avg_table = solver.getAvgTable()
    print('TBI solved %d stations (%d ODE variables) in %.1f seconds.'
          % (model.getNumberOfStations(), M * (M - 1) * KPH + M * 4, time.time() - t0))
    print(avg_table.head(M))   # queue-length summary at the queueing stations

    # For comparison, the undecomposed solution of the same model:
    #   solver = FLD(model, method='closing', timespan=[0, TEND], stiff=True)
    # takes several minutes on the same machine (about 80 seconds already at
    # M=8, and beyond 10 minutes at M=12).
