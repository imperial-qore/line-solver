"""
Finite Capacity Region with Dropping (Multiclass Open Network)

This example creates a multiclass open network with 2 queues and a
finite capacity region with global and per-class constraints.
When the region is full, arriving jobs are dropped (lost).
"""

from line_solver import *

# Create network
model = Network('FCR Dropping Example')

# Add nodes
source = Source(model, 'Source')
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')

# Add job classes
class1 = OpenClass(model, 'Class1', 0)
class2 = OpenClass(model, 'Class2', 1)

# Set arrival and service rates
source.set_arrival(class1, Exp.fit_rate(0.4))
source.set_arrival(class2, Exp.fit_rate(0.3))
queue1.set_service(class1, Exp.fit_rate(1.0))
queue1.set_service(class2, Exp.fit_rate(0.9))
queue2.set_service(class1, Exp.fit_rate(1.1))
queue2.set_service(class2, Exp.fit_rate(1.0))

# Create routing matrix
P = model.init_routing_matrix()
# Class1 routing
P.set(class1, class1, source, queue1, 0.5)
P.set(class1, class1, source, queue2, 0.5)
P.set(class1, class1, queue1, queue2, 0.3)
P.set(class1, class1, queue1, sink, 0.7)
P.set(class1, class1, queue2, sink, 1.0)

# Class2 routing
P.set(class2, class2, source, queue1, 0.6)
P.set(class2, class2, source, queue2, 0.4)
P.set(class2, class2, queue1, queue2, 0.5)
P.set(class2, class2, queue1, sink, 0.5)
P.set(class2, class2, queue2, sink, 1.0)

model.link(P)

# Add finite capacity region with constraints
# When region is full, jobs are dropped (lost)
fcr = model.add_region('FCR', queue1, queue2)
fcr.set_global_max_jobs(8)          # Global: max 8 jobs in region
fcr.set_class_max_jobs(class1, 5)   # Class1: max 5 jobs
fcr.set_class_max_jobs(class2, 4)   # Class2: max 4 jobs
fcr.set_drop_rule(class1, True)     # True = drop jobs
fcr.set_drop_rule(class2, True)     # True = drop jobs

# Run JMT
solver = JMT(model, seed=23000, samples=50000, verbose=VerboseLevel.SILENT)

# Get results
avg_table = solver.get_avg_table()
print(avg_table)
