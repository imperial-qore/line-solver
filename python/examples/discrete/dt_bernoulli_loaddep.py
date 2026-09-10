"""
Load-dependent Bernoulli server, and the discrete-time arrival theorem.

Example 2.10 of Daduna (2001) notes that a discrete-time M/M/c queue has no
exactly equivalent state dependent single server, but that p(n) = p min(n,c)
reproduces its conditional service intensity. That is a load dependence in LINE
terms, so the model is a single Bernoulli server whose service probability
rises with the queue length up to c servers' worth of capacity.

The example also prints the arrival distribution of theorem 2.11, the law an
arriving job sees with itself not counted. In continuous time with Poisson
arrivals that law would coincide with the time-stationary one (PASTA); discrete
time has no such analogue, and the two rows below differ.
"""
import numpy as np

from line_solver import *
from line_solver.api.dqsys import dqsys_bernoulli1

if __name__ == "__main__":
    a = 0.6   # per-slot arrival probability
    s = 0.3   # per-slot service completion probability of one server
    c = 3     # servers' worth of capacity
    L = 20    # buffer capacity, which bounds the state space

    model = Network('LoadDepBernoulli')

    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    job_class = OpenClass(model, 'Class1')
    source.setArrival(job_class, Geometric(a))
    queue.setService(job_class, Geometric(s))
    queue.setCapacity(L)
    alpha = np.minimum(np.arange(1, L + 1), c).astype(float)
    queue.setLoadDependence(alpha)          # p(n) = s * min(n,c), example 2.10

    model.link(Network.serialRouting(source, queue, sink))

    solver = SolverNC(model)
    solver.options.config = {'slotted': True}
    print(solver.getAvgTable())

    r = dqsys_bernoulli1(a, s * alpha, L)
    print('time-stationary law pi(0..6)  : %s' % np.round(r['pmf'][:7], 6))
    print('arrival law         pi_1(0..6): %s' % np.round(r['arrivalPmf'][:7], 6))
    print('no PASTA in discrete time: the two rows above are different laws')
