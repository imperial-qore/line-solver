"""
Open queueing network with an NHPP (cyclic) arrival process.

NHPP is a non-homogeneous Poisson process with a piecewise-constant intensity:
segment i covers [breakpoints[i], breakpoints[i+1]) and carries rate rates[i].
With cyclic=True the schedule repeats with period breakpoints[-1]-breakpoints[0],
giving a cyclic Poisson process. The LDES simulation engine honours the exact
schedule; SolverFLD honours it in getTranAvg, where the intensity enters the
closing fluid ODE as a time-varying rate multiplier. Steady state is the
time-average rate.
"""

import numpy as np

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    jobclass = OpenClass(model, 'OpenClass', 0)

    # Rates 2,8,4 held for 3,1,2 time units, repeating cyclically.
    source.set_arrival(jobclass, NHPP([0, 3, 4, 6], [2, 8, 4], True))
    queue.set_service(jobclass, Exp(10))

    model.link(Network.serial_routing([source, queue, sink]))

    avg_table_1 = LDES(model, seed=1234, samples=100000).avg_table()
    print('DES Result:')
    print(avg_table_1)

    # Fluid transient over two periods: the queue throughput tracks lambda(t).
    nhpp = source.getArrivalProcess(jobclass)
    _, _, TNt = SolverFLD(model, timespan=[0, 12], verbose=0).getTranAvg()
    t = np.asarray(TNt[1][0].t).ravel()
    y = np.asarray(TNt[1][0].metric).ravel()
    print('FLD transient:')
    for tt in [1.5, 3.5, 5.0, 7.5, 9.5, 11.0]:
        print('  t=%5.2f lambda(t)=%6.3f queue Tput=%6.3f'
              % (tt, nhpp.getRateAt(tt), float(np.interp(tt, t, y))))
