# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Create network
model = Network('model')

# Block 1: nodes
node = np.empty(3, dtype=object)
node[0] = Delay(model, 'Delay')
node[1] = Queue(model, 'Queue1', SchedStrategy.PS)
node[2] = Queue(model, 'Queue2', SchedStrategy.PS)

# Block 2: classes (1 job each)
jobclass = np.empty(2, dtype=object)
jobclass[0] = ClosedClass(model, 'Class1', 1, node[0], 0)
jobclass[1] = ClosedClass(model, 'Class2', 1, node[0], 0)

# %%
# Set service distributions
node[0].set_service(jobclass[0], Erlang(3, 2))
node[0].set_service(jobclass[1], HyperExp(0.5, 3.0, 10.0))

node[1].set_service(jobclass[0], HyperExp(0.1, 1.0, 10.0))
node[1].set_service(jobclass[1], MMPP2(1, 2, 3, 4))

node[2].set_service(jobclass[0], HyperExp(0.1, 1.0, 10.0))
node[2].set_service(jobclass[1], Erlang(1, 2))

# %%
# Set up connectivity
model.add_link(node[0], node[0])
model.add_link(node[0], node[1])
model.add_link(node[0], node[2])
model.add_link(node[1], node[0])
model.add_link(node[2], node[0])

# Class1: Probabilistic routing from Delay
node[0].set_prob_routing(jobclass[0], node[0], 0.0)
node[0].set_prob_routing(jobclass[0], node[1], 0.3)
node[0].set_prob_routing(jobclass[0], node[2], 0.7)
node[1].set_prob_routing(jobclass[0], node[0], 1.0)
node[2].set_prob_routing(jobclass[0], node[0], 1.0)

# Class2: Random routing strategy
node[0].set_routing(jobclass[1], RoutingStrategy.RAND)
node[1].set_routing(jobclass[1], RoutingStrategy.RAND)
node[2].set_routing(jobclass[1], RoutingStrategy.RAND)

# %%
# Configure solvers
solver = np.array([], dtype=object)
solver = np.append(solver, JMT(model, seed=23000))
solver = np.append(solver, LDES(model, seed=23000))

# %%
# Solve with each solver
for s in range(len(solver)):
    print(f'SOLVER: {solver[s].get_name()}')
    avg_table = solver[s].avg_table()
    print(avg_table)