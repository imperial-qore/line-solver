"""
State-Dependent Routing - Join the Shortest Queue

This example demonstrates JSQ routing in an open queueing network.
The router sends each arriving job to the queue with the fewest jobs.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from line_solver import *

model = Network('myModel')

N = 3
rho = 0.7
mu = 1
lam = N * rho * mu

# Block 1: nodes
source = Source(model, 'Source')
router = Router(model, 'Router')
queue = []
for i in range(N):
    queue.append(Queue(model, f'Queue{i+1}', SchedStrategy.FCFS))
sink = Sink(model, 'Sink')

# Block 2: classes
oclass = OpenClass(model, 'Class1')
source.setArrival(oclass, Exp(lam))
for i in range(N):
    queue[i].setService(oclass, Exp(mu))

# Block 3: topology
model.addLink(source, router)
for i in range(N):
    model.addLink(router, queue[i])
    model.addLink(queue[i], sink)

router.setRouting(oclass, RoutingStrategy.JSQ)

solvers = [
    JMT(model, seed=23000),
    LDES(model, seed=23000),
]

AvgTable = []
for solver in solvers:
    print(f'SOLVER: {solver.getName()}')
    avg_table = solver.getAvgNodeTable()
    print(avg_table)
    AvgTable.append(avg_table)
    print()
