"""
Pass-and-swap (PAS) station as a multiserver order-independent (M/M/K) queue.

The M/M/K queue satisfies the order-independence conditions (Dorsman & Gardner
2024, Sect. 2.2): with K unit-rate servers the total service rate in state c is
mu(c) = min(n, K), where n is the number of jobs. With an empty swapping graph
the PAS station reduces to a plain order-independent queue.
"""

import numpy as np
from line_solver import *

K = 2   # number of servers


def mu_fun(c):
    c = np.asarray(c)
    return float(min(len(c), K))


model = Network('PASmmk')
source = Source(model, 'Source')
queue = Queue(model, 'PASQueue', SchedStrategy.PAS)
sink = Sink(model, 'Sink')

lam = [0.7, 0.5]
classes = [OpenClass(model, f'Class{r+1}') for r in range(2)]
for r in range(2):
    source.setArrival(classes[r], Exp(lam[r]))

queue.setService(mu_fun)
queue.setSwapGraph(np.zeros((2, 2)))   # empty graph: plain order-independent queue
queue.setNumberOfServers(K)
queue.setCap(4)

P = model.initRoutingMatrix()
for r in range(2):
    P[classes[r]] = Network.serialRouting(source, queue, sink)
model.link(P)

AvgTable = CTMC(model, cutoff=4).getAvgTable()
print(AvgTable)
