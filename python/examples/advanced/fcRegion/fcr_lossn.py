"""
Loss Network with Finite Capacity Region (NC Solver)

This example demonstrates the NC solver's ability to analyze open
loss networks using the Erlang fixed-point approximation.
The model has a multiclass Delay node inside an FCR with DROP policy.
"""

from line_solver import *

# Create network
model = Network('FCR Loss Network')

# Add nodes
source = Source(model, 'Source')
delay = Delay(model, 'Delay')
sink = Sink(model, 'Sink')

# Add job classes
class1 = OpenClass(model, 'Class1', 0)
class2 = OpenClass(model, 'Class2', 1)

# Set arrival and service rates
lambda1 = 0.3
lambda2 = 0.2
mu1 = 1.0
mu2 = 0.8

source.set_arrival(class1, Exp.fit_rate(lambda1))
source.set_arrival(class2, Exp.fit_rate(lambda2))
delay.set_service(class1, Exp.fit_rate(mu1))
delay.set_service(class2, Exp.fit_rate(mu2))

# Create routing matrix
P = model.init_routing_matrix()
P.set(class1, class1, source, delay, 1.0)
P.set(class1, class1, delay, sink, 1.0)
P.set(class2, class2, source, delay, 1.0)
P.set(class2, class2, delay, sink, 1.0)
model.link(P)

# Add finite capacity region with constraints
# When region is full, arriving jobs are dropped (lost)
fcr = model.add_region(delay)
fcr.set_global_max_jobs(5)           # Global: max 5 jobs in region
fcr.set_class_max_jobs(class1, 3)    # Class1: max 3 jobs
fcr.set_class_max_jobs(class2, 3)    # Class2: max 3 jobs
fcr.set_drop_rule(class1, True)      # True = drop jobs
fcr.set_drop_rule(class2, True)      # True = drop jobs

# Run NC solver (uses Erlang fixed-point for loss networks)
print('Running NC solver (lossn method)...')
solver_nc = NC(model)
avg_table_nc = solver_nc.get_avg_table()

print('\nNC Results:')
print(avg_table_nc)

# Run JMT for comparison
print('\nRunning JMT for comparison...')
solver_jmt = JMT(model, seed=23000, samples=500000)
avg_table_jmt = solver_jmt.get_avg_table()

print('\nJMT Results:')
print(avg_table_jmt)
