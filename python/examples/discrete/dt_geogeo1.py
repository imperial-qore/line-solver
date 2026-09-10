"""
Geo/Geo/1: the discrete-time single-server queue with an unbounded buffer.

Arrivals occur with probability a in each slot, a service completes with
probability s. SolverNC recognizes the slotted model and returns the exact
closed form of Daduna (2001), theorem 2.3, which for constant a and s is the
geometric law of corollary 2.7. Mean queue length a(1-a)/(s-a) and mean
sojourn time (1-a)/(s-a) slots.
"""
from line_solver import *
from line_solver.api.dqsys import dqsys_bernoulli1

if __name__ == "__main__":
    a = 0.2   # per-slot arrival probability
    s = 0.5   # per-slot service completion probability

    model = Network('GeoGeo1')

    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    job_class = OpenClass(model, 'Class1')
    source.setArrival(job_class, Geometric(a))
    queue.setService(job_class, Geometric(s))

    model.link(Network.serialRouting(source, queue, sink))

    solver = SolverNC(model)
    solver.options.config = {'slotted': True}   # run on the slot lattice
    print(solver.getAvgTable())

    print('closed form: E[N] = %g, E[T] = %g slots'
          % (a * (1 - a) / (s - a), (1 - a) / (s - a)))

    # The same numbers straight from the single-queue formula.
    r = dqsys_bernoulli1(a, s)
    print('E[N] = %g, U = %g, X = %g per slot'
          % (r['meanQueueLength'], r['utilization'], r['throughput']))
