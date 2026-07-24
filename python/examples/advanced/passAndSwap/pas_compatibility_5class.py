"""
Pass-and-swap (PAS) queue: an order-independent queue with a swapping graph.

Reproduces the five-class, three-server compatibility example of
Dorsman & Gardner (2024), "New directions in pass-and-swap queues",
Queueing Systems 107:205-256 (Figs 1-2). A PAS station is parameterized by the
total service-rate function mu(c) of the ordered state vector c (0-based class
indices, c[0] the oldest job) and by a swapping graph G; both are properties of
the Queue object. The stationary distribution is the order-independent product
form and is invariant to G.
"""

import numpy as np
from line_solver import *

# Server compatibility (Fig 1): three unit-rate servers over five classes.
# mu(c) = number of servers compatible with at least one class present in c.
comp = np.array([[1, 0, 0, 1, 0],    # server 1: classes {0,3}
                 [0, 1, 0, 1, 0],    # server 2: classes {1,3}
                 [0, 0, 1, 0, 1]])   # server 3: classes {2,4}


def mu_fun(c):
    c = np.asarray(c, dtype=int)
    return float(np.sum(np.any(comp[:, c], axis=1))) if len(c) else 0.0


model = Network('PAScompatibility')
source = Source(model, 'Source')
queue = Queue(model, 'PASQueue', SchedStrategy.PAS)
sink = Sink(model, 'Sink')

lam = [0.5, 0.4, 0.3, 0.2, 0.1]
classes = [OpenClass(model, f'Class{r+1}') for r in range(5)]
for r in range(5):
    source.setArrival(classes[r], Exp(lam[r]))

# PAS service: specify the rate function mu(c) as a whole (no per-class rates)
queue.setService(mu_fun)
# Swapping graph (Fig 2a): E = {(0,2),(0,4),(1,3),(2,3),(3,4)} (0-based)
G = np.zeros((5, 5))
for i, j in [(0, 2), (0, 4), (1, 3), (2, 3), (3, 4)]:
    G[i, j] = 1
    G[j, i] = 1
queue.setSwapGraph(G)
queue.setNumberOfServers(3)   # three servers (used for utilization reporting)
queue.setCap(3)               # finite buffer (order-independent loss model)

P = model.initRoutingMatrix()
for r in range(5):
    P[classes[r]] = Network.serialRouting(source, queue, sink)
model.link(P)

AvgTable = CTMC(model, cutoff=3).getAvgTable()
print(AvgTable)
