# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
model = Network('model')

node = np.empty(2, dtype=object)
node[0] = Delay(model, 'Delay')
node[1] = Queue(model, 'Queue1', SchedStrategy.PS)

jobclass = np.empty(2, dtype=object)
jobclass[0] = ClosedClass(model, 'Class1', 2, node[0], 0)
jobclass[1] = ClosedClass(model, 'Class2', 2, node[0], 0)

node[0].set_service(jobclass[0], Erlang(3, 2))
node[0].set_service(jobclass[1], HyperExp(0.5, 3.0, 10.0))

node[1].set_service(jobclass[0], HyperExp(0.1, 1.0, 10.0))
node[1].set_service(jobclass[1], Exp(1))

P = model.init_routing_matrix()
P.set(jobclass[0], jobclass[0], node[0], node[0], 0.3)
P.set(jobclass[0], jobclass[0], node[0], node[1], 0.1)
P.set(jobclass[0], jobclass[0], node[1], node[0], 0.2)

P.set(jobclass[0], jobclass[1], node[0], node[0], 0.6)
P.set(jobclass[0], jobclass[1], node[1], node[0], 0.8)

P.set(jobclass[1], jobclass[0], node[1], node[0], 1.0)

P.set(jobclass[1], jobclass[1], node[0], node[1], 1.0)

model.link(P)
# %%
# Aligned with JAR test scenarios for cqn_twoclass_hyperl
# JAR tests: CTMC(), JMT(seed=23000, samples=5000), SSA(seed=23000, samples=5000),
#           Fluid(), MVA(method="exact"), NC(method="exact"), MAM(), LDES(seed=23000, samples=5000)

solver = np.array([], dtype=object)

# CTMC with default settings (matches JAR)
solver = np.append(solver, CTMC(model))

# JMT with seed=23000, samples=5000 (matches JAR)
solver = np.append(solver, JMT(model, seed=23000, samples=5000))

# SSA with seed=23000, samples=5000 (matches JAR)
solver = np.append(solver, SSA(model, seed=23000, samples=5000))

# Fluid with default settings (matches JAR)
solver = np.append(solver, FLD(model))

# MVA with method="exact" (matches JAR)
solver = np.append(solver, MVA(model, method='exact'))

# NC with method="exact" (matches JAR)
solver = np.append(solver, NC(model, method='exact'))

# MAM with default settings (matches JAR)
solver = np.append(solver, MAM(model))

# LDES with seed=23000, samples=5000 (matches MATLAB)
solver = np.append(solver, LDES(model, seed=23000, samples=5000))

AvgTable = np.empty(len(solver), dtype=object)
for s in range(len(solver)):
    print(f'\nSOLVER: {solver[s].get_name()}')
    AvgTable[s] = solver[s].avg_table()
    print(AvgTable[s])