"""
Finite Capacity Region - Constraint Types

This example demonstrates different FCR constraint types:
- Global max jobs: limits total jobs across all classes
- Per-class max jobs: limits jobs of a specific class
- Drop rule: determines blocking vs dropping behavior
"""

from line_solver import *

# Create a multiclass open network
model = Network('FCR Constraints Demo')

# Nodes
source = Source(model, 'Source')
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')

# Job classes
high_priority = OpenClass(model, 'HighPriority', 0)
low_priority = OpenClass(model, 'LowPriority', 1)

# Arrival and service rates
source.set_arrival(high_priority, Exp.fit_rate(0.3))  # High priority: lower arrival rate
source.set_arrival(low_priority, Exp.fit_rate(0.5))   # Low priority: higher arrival rate
queue1.set_service(high_priority, Exp.fit_rate(1.0))
queue1.set_service(low_priority, Exp.fit_rate(0.8))
queue2.set_service(high_priority, Exp.fit_rate(1.2))
queue2.set_service(low_priority, Exp.fit_rate(1.0))

# Routing: both classes go through both queues
P = model.init_routing_matrix()
P.set(high_priority, high_priority, source, queue1, 1.0)
P.set(high_priority, high_priority, queue1, queue2, 1.0)
P.set(high_priority, high_priority, queue2, sink, 1.0)
P.set(low_priority, low_priority, source, queue1, 1.0)
P.set(low_priority, low_priority, queue1, queue2, 1.0)
P.set(low_priority, low_priority, queue2, sink, 1.0)
model.link(P)

# Add Finite Capacity Region with multiple constraints
fcr = model.add_region('FCR', queue1, queue2)

# Global constraint: max 6 jobs total in the region
fcr.set_global_max_jobs(6)

# Per-class constraints: high priority gets more space
fcr.set_class_max_jobs(high_priority, 4)  # HighPriority: max 4 jobs
fcr.set_class_max_jobs(low_priority, 3)   # LowPriority: max 3 jobs
# Note: per-class limits must sum to >= global limit for consistent behavior

# Drop rules: high priority jobs wait, low priority jobs are dropped
fcr.set_drop_rule(high_priority, False)  # False = block (wait)
fcr.set_drop_rule(low_priority, True)    # True = drop

# Display constraint summary
print('\n=== FCR Constraint Configuration ===')
print('Region covers: Queue1, Queue2')
print('Global max jobs: 6')
print('HighPriority max: 4 (blocking)')
print('LowPriority max: 3 (dropping)\n')

# Solve with JMT
solver = JMT(model, seed=23000, samples=100000, verbose=VerboseLevel.SILENT)
avg_table = solver.get_avg_table()
print(avg_table)

print('\nNote: HighPriority jobs experience delays when region is full.')
print('LowPriority jobs are dropped when their class limit or global limit is reached.')
