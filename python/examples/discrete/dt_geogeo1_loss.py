"""
Geo/Geo/1/L: the discrete-time loss system of Daduna (2001), corollary 2.8.

The buffer holds at most L jobs. An arrival that lands in a slot which already
holds L jobs is lost, i.e. b(n) = 0 for n >= L, which is exactly the assumption
corollary 2.8 places on the arrival probabilities. The queue length law stays
the birth-death form of theorem 2.3, now normalized over the finite state
space, and the loss probability is the stationary probability of a full system.
"""
import numpy as np

from line_solver import *
from line_solver.api.dqsys import dqsys_bernoulli1

if __name__ == "__main__":
    a = 0.2   # per-slot arrival probability
    s = 0.5   # per-slot service completion probability
    L = 4     # buffer capacity in jobs

    model = Network('GeoGeo1L')

    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    job_class = OpenClass(model, 'Class1')
    source.setArrival(job_class, Geometric(a))
    queue.setService(job_class, Geometric(s))
    queue.setCapacity(L)

    model.link(Network.serialRouting(source, queue, sink))

    solver = SolverNC(model)
    solver.options.config = {'slotted': True}
    print(solver.getAvgTable())

    r = dqsys_bernoulli1(a, s, L)
    print('queue length law   : %s' % np.round(r['pmf'], 6))
    print('loss probability   : %g' % r['lossProb'])
    print('carried throughput : %g of the %g offered per slot' % (r['throughput'], a))
