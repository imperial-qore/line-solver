#!/usr/bin/env python3
"""Transient analysis of a MAP/MAP/1 queue via the Laplace-domain transient QBD
solver in MAM (native Python).

Computes the time-dependent mean queue length, utilization, and throughput of a
single-server queue with correlated (MAP) arrivals and correlated (MAP) service,
starting empty. The Laplace transient QBD method is auto-selected by getTranAvg
because arrival and service are correlated MAPs (the libQBD/expm fast path only
handles Poisson arrivals).
"""

import numpy as np
from line_solver import *
from line_solver.api.mam.map_analysis import map_lambda


def mam_transient_mapmap1():
    D0 = np.array([[-8., 1, 3], [0, -6, 4], [2, 0, -3]])
    D1 = np.array([[3., 1, 0], [0, 2, 0], [0, 0, 1]])
    S0 = np.array([[-3., 1], [6, -7]])
    S1 = np.array([[0., 2], [1, 0]])
    # Scale service to rho = 0.6 (native map_scale multiplies the generator).
    fac = (map_lambda(D0, D1) / 0.6) / map_lambda(S0, S1)
    Ds0, Ds1 = fac * S0, fac * S1

    model = Network('MAP/MAP/1 transient')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1')
    source.setArrival(oclass, MAP(D0, D1))
    queue.setService(oclass, MAP(Ds0, Ds1))
    model.link(Network.serial_routing([source, queue, sink]))

    solver = MAM(model, timespan=[0, 40])
    QNt, UNt, TNt = solver.getTranAvg()

    q = QNt[1][0]
    print('MAP/MAP/1 transient (rho=0.6), start empty:')
    print('  t=%5.1f  E[N]=%.5f  U=%.5f  Tput=%.5f'
          % (q.t[-1], q.metric[-1], UNt[1][0].metric[-1], TNt[1][0].metric[-1]))
    print('  steady-state E[N] approaches 1.5458 as t -> inf.')
    return model


if __name__ == '__main__':
    mam_transient_mapmap1()
