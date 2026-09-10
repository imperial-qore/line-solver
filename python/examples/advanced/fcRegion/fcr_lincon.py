"""
Linear Admission Constraints for Finite Capacity Regions

Open network with a Delay node inside an FCR with DROP policy and
general linear admission constraints: A*n <= b.

Constraint: 2*n_1 + 1*n_2 <= 5, 1*n_1 + 3*n_2 <= 7
(general weighted cross-class coupling)
"""

import numpy as np
from line_solver import *

model = Network('FCR LinCon')

source = Source(model, 'Source')
delay = Delay(model, 'Delay')
sink = Sink(model, 'Sink')

class1 = OpenClass(model, 'Class1', 0)
class2 = OpenClass(model, 'Class2', 0)

source.set_arrival(class1, Exp(0.3))
source.set_arrival(class2, Exp(0.2))
delay.set_service(class1, Exp(1.0))
delay.set_service(class2, Exp(0.8))

P = model.init_routing_matrix()
P.set(class1, class1, source, delay, 1.0)
P.set(class1, class1, delay, sink, 1.0)
P.set(class2, class2, source, delay, 1.0)
P.set(class2, class2, delay, sink, 1.0)
model.link(P)

# Add FCR with DROP policy and linear constraints
fcr = model.add_region(delay)
fcr.set_class_max_jobs(class1, 1000)
fcr.set_class_max_jobs(class2, 1000)
fcr.set_drop_rule(class1, DropStrategy.DROP)
fcr.set_drop_rule(class2, DropStrategy.DROP)
fcr.set_global_max_jobs(1000)
fcr.set_linear_constraints(
    np.array([[2.0, 1.0], [1.0, 3.0]]),
    np.array([5.0, 7.0])
)

solver = LDES(model, seed=23000, samples=500000)
avg_table = solver.avg_table()
print(avg_table)
