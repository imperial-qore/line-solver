"""
Pass-and-swap (PAS) queue with a swapping graph that includes a self-loop.

Per Dorsman & Gardner (2024), Sect. 2.3, self-loops in the swapping graph are
permissible: a completing class-i job may take the place of another class-i job
further back in the queue. The order-independent product form (Theorem 2) still
holds and remains invariant to the swapping graph.
"""

import numpy as np
from line_solver import *

# Three unit-rate servers, one per class: mu(c) = number of distinct classes in c.
comp = np.eye(3)


def mu_fun(c):
    c = np.asarray(c, dtype=int)
    return float(np.sum(np.any(comp[:, c], axis=1))) if len(c) else 0.0


model = Network('PASselfloop')
source = Source(model, 'Source')
queue = Queue(model, 'PASQueue', SchedStrategy.PAS)
sink = Sink(model, 'Sink')

lam = [0.6, 0.4, 0.3]
classes = [OpenClass(model, f'Class{r+1}') for r in range(3)]
for r in range(3):
    source.setArrival(classes[r], Exp(lam[r]))

queue.setService(mu_fun)
# Swapping graph with a self-loop on class 0 and an edge (1,2) (0-based)
G = np.zeros((3, 3))
G[0, 0] = 1            # self-loop: class 0 swaps with class 0
G[1, 2] = 1
G[2, 1] = 1
queue.setSwapGraph(G)
queue.setNumberOfServers(3)
queue.setCap(3)

P = model.initRoutingMatrix()
for r in range(3):
    P[classes[r]] = Network.serialRouting(source, queue, sink)
model.link(P)

AvgTable = CTMC(model, cutoff=3).getAvgTable()
print(AvgTable)
