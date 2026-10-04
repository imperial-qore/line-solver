"""
Multichain closed cycle of Bernoulli servers.

Section 3.2 of Daduna (2001) splits the circulating jobs into chains that never
mix. Service in the cycle is type independent and FCFS forbids overtaking, so
the cyclic order of the jobs is frozen for all time and the joint queue length
law is the unichain one at the aggregate population. Each chain then holds a
share of every station equal to its share of the population, which is what
SolverNC reports below: the per-station totals reproduce dt_cycle at
N = N_1 + N_2, and the per-class rows split them in the ratio N_1 : N_2.
"""
import numpy as np

from line_solver import *
from line_solver.api.dpfqn import dpfqn_nc

if __name__ == "__main__":
    p = [0.5, 0.25, 0.7]
    N1 = 3   # jobs in chain 1
    N2 = 2   # jobs in chain 2
    J = len(p)

    model = Network('BernoulliCycleMC')

    station = [Queue(model, 'Queue%d' % (j + 1), SchedStrategy.FCFS) for j in range(J)]

    class1 = ClosedClass(model, 'Chain1', N1, station[0])
    class2 = ClosedClass(model, 'Chain2', N2, station[0])
    for j in range(J):
        station[j].setService(class1, Geometric(p[j]))
        station[j].setService(class2, Geometric(p[j]))

    # Both chains follow the same cycle and never switch class.
    cycle = np.roll(np.eye(J), 1, axis=1)
    zero = np.zeros((J, J))
    P = model.initRoutingMatrix()
    P[0][0] = cycle
    P[0][1] = zero
    P[1][0] = zero
    P[1][1] = cycle
    model.link(P)

    solver = SolverNC(model)
    solver.options.config = {'slotted': True}
    print(solver.getAvgTable())

    lG, G, G1 = dpfqn_nc(p, N1 + N2)
    X = G1[N1 + N2] / G
    print('aggregate throughput  = %g jobs per slot' % X)
    print('per-chain throughput  = %g and %g' % (X * N1 / (N1 + N2), X * N2 / (N1 + N2)))
