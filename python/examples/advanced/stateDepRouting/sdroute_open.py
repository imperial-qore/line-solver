"""
State-Dependent Routing - Open Network

This example demonstrates state-dependent routing in an open queueing network.
The router uses round-robin routing to alternate between two queues.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create open network with state-dependent routing
model = Network('myModel')

# Block 1: nodes
source = Source(model, 'Source')
router = Router(model, 'Router')
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')

# Block 2: classes
oclass = OpenClass(model, 'Class1')
source.setArrival(oclass, Exp(1))
queue1.setService(oclass, Exp(2))
queue2.setService(oclass, Exp(2))

# Block 3: topology
model.addLink(source, router)
model.addLink(router, queue1)
model.addLink(router, queue2)
model.addLink(queue1, sink)
model.addLink(queue2, sink)

# Set state-dependent routing: round-robin between Queue1 and Queue2
router.setRouting(oclass, RoutingStrategy.RROBIN)

print('=== State-Dependent Routing - Open Network ===\n')

# Solve with JMT
print('JMT Simulation:')
solver_jmt = JMT(model, seed=23000)
avg_table_jmt = solver_jmt.getAvgNodeTable()
print(avg_table_jmt)

# Solve with CTMC (with cutoff to limit state space)
print('\nCTMC (with cutoff=5):')
solver_ctmc = CTMC(model, cutoff=5)
avg_table_ctmc = solver_ctmc.getAvgNodeTable()
print(avg_table_ctmc)

print('\nNote: Round-robin routing alternates jobs between Queue1 and Queue2.')
print('      This is state-dependent because routing decisions depend on the')
print('      number of jobs previously routed to each destination.')
