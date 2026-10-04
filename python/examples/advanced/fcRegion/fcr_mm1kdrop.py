"""
Finite Capacity Region with Dropping vs M/M/1/K

This example shows that an FCR with dropping around a single queue
behaves like an M/M/1/K queue, where K is the FCR capacity.
Jobs arriving when the region is full are dropped (lost).
"""

from line_solver import *

# Parameters
arrival_rate = 0.8
service_rate = 1.0
K = 3

# Model 1: Queue with FCR (dropping)
model1 = Network('FCR Dropping')

source1 = Source(model1, 'Source')
queue1 = Queue(model1, 'Queue', SchedStrategy.FCFS)
sink1 = Sink(model1, 'Sink')

jobclass1 = OpenClass(model1, 'Class1', 0)

source1.set_arrival(jobclass1, Exp.fit_rate(arrival_rate))
queue1.set_service(jobclass1, Exp.fit_rate(service_rate))

P1 = model1.init_routing_matrix()
P1.set(jobclass1, jobclass1, source1, queue1, 1.0)
P1.set(jobclass1, jobclass1, queue1, sink1, 1.0)
model1.link(P1)

# Add FCR with dropping - jobs are lost when region is full
fcr = model1.add_region(queue1)
fcr.set_global_max_jobs(K)
fcr.set_drop_rule(jobclass1, True)  # true = drop jobs

# Model 2: M/M/1/K using queue capacity
model2 = Network('M/M/1/K')

source2 = Source(model2, 'Source')
queue2 = Queue(model2, 'Queue', SchedStrategy.FCFS)
queue2.set_number_of_servers(1)
queue2.set_capacity(K)  # Set queue capacity to K
sink2 = Sink(model2, 'Sink')

jobclass2 = OpenClass(model2, 'Class1', 0)

source2.set_arrival(jobclass2, Exp.fit_rate(arrival_rate))
queue2.set_service(jobclass2, Exp.fit_rate(service_rate))

P2 = model2.init_routing_matrix()
P2.set(jobclass2, jobclass2, source2, queue2, 1.0)
P2.set(jobclass2, jobclass2, queue2, sink2, 1.0)
model2.link(P2)

# Solve both models and compare results
solver1 = JMT(model1, seed=23000, samples=100000)
avgTable1 = solver1.get_avg_table()
print(avgTable1)

solver2 = JMT(model2, seed=23000, samples=100000)
avgTable2 = solver2.get_avg_table()
print(avgTable2)
