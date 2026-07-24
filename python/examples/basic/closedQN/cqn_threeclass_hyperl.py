# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.SILENT)
# %%
model = Network('model')

node = np.empty(2, dtype=object)
node[0] = Delay(model, 'Delay')
node[1] = Queue(model, 'Queue1', SchedStrategy.PS)
node[1].set_number_of_servers(2)

jobclass = np.empty(3, dtype=object)
jobclass[0] = ClosedClass(model, 'Class1', 2, node[0], 0)
jobclass[1] = ClosedClass(model, 'Class2', 0, node[0], 0)
jobclass[2] = ClosedClass(model, 'Class3', 1, node[0], 0)

node[0].set_service(jobclass[0], Erlang(3, 2))
node[0].set_service(jobclass[1], HyperExp(0.5, 3.0, 10.0))
node[0].set_service(jobclass[2], Exp(1))

node[1].set_service(jobclass[0], HyperExp(0.1, 1.0, 10.0))
node[1].set_service(jobclass[1], Exp(2))
node[1].set_service(jobclass[2], Exp(3))

P = model.init_routing_matrix()

# P{1,1} = [0.3,0.1; 0.2,0]
P.set(jobclass[0], jobclass[0], node[0], node[0], 0.3)
P.set(jobclass[0], jobclass[0], node[0], node[1], 0.1)
P.set(jobclass[0], jobclass[0], node[1], node[0], 0.2)

# P{1,2} = [0.6,0; 0.8,0]
P.set(jobclass[0], jobclass[1], node[0], node[0], 0.6)
P.set(jobclass[0], jobclass[1], node[1], node[0], 0.8)

# P{2,1} = [0,0; 1,0]
P.set(jobclass[1], jobclass[0], node[1], node[0], 1.0)

# P{2,2} = [0,1; 0,0]
P.set(jobclass[1], jobclass[1], node[0], node[1], 1.0)

# P{3,3} = circul(M)
P.set(jobclass[2], jobclass[2], node[0], node[1], 1.0)
P.set(jobclass[2], jobclass[2], node[1], node[0], 1.0)

model.link(P)
# %%
solver = np.array([], dtype=object)
solver = np.append(solver, CTMC(model))
solver = np.append(solver, JMT(model, seed=23000, samples=5000, verbose=False))
solver = np.append(solver, SSA(model, seed=23000, samples=5000, verbose=False))
solver = np.append(solver, FLD(model))
solver = np.append(solver, MVA(model))
solver = np.append(solver, NC(model, method='exact'))
solver = np.append(solver, MAM(model))
solver = np.append(solver, LINE(model, seed=23000))
solver = np.append(solver, LDES(model, seed=23000, samples=5000, verbose=False))

AvgTable = np.empty(len(solver), dtype=object)
AvgChainTable = np.empty(len(solver), dtype=object)
AvgSysTable = np.empty(len(solver), dtype=object)

for s in range(len(solver)):
    print(f'\nSOLVER: {solver[s].get_name()}')
    try:
        AvgTable[s] = solver[s].avg_table()
        print(AvgTable[s])
        AvgChainTable[s] = solver[s].avg_chain_table()
        AvgSysTable[s] = solver[s].avg_sys_table()
    except Exception as e:
        print(f'Error with {solver[s].get_name()}: {str(e)[:100]}')
        # Set to None for failed solvers
        AvgTable[s] = None
        AvgChainTable[s] = None
        AvgSysTable[s] = None