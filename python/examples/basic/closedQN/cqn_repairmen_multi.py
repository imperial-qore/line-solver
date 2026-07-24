# %%
from line_solver import *

import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
model = Network('model')

node = np.empty(2, dtype=object)
node[0] = Delay(model, 'Delay')
node[1] = Queue(model, 'Queue1', SchedStrategy.FCFS)
node[1].set_number_of_servers(3)

# Default: scheduling is set as FCFS everywhere, routing as Random
jobclass = np.empty(2, dtype=object)
jobclass[0] = ClosedClass(model, 'Class1', 4, node[0], 0)
jobclass[1] = ClosedClass(model, 'Class2', 2, node[0], 0)

node[0].set_service(jobclass[0], Exp(1))
node[0].set_service(jobclass[1], Exp(1))

node[1].set_service(jobclass[0], Exp(1))
node[1].set_service(jobclass[1], Exp(10))

myP = model.init_routing_matrix()
myP.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)
myP.set(jobclass[0], jobclass[0], node[1], node[0], 1.0)
myP.set(jobclass[1], jobclass[1], node[0], node[1], 1.0)
myP.set(jobclass[1], jobclass[1], node[1], node[0], 1.0)
model.link(myP)
# %%
# Solve with multiple solvers
solver = np.array([], dtype=object)
solver = np.append(solver, CTMC(model))

# Skip QNS solvers if external qnsolver is not available
if QNS.isAvailable():
    solver = np.append(solver, QNS(model, method='conway'))
    solver = np.append(solver, QNS(model, method='reiser'))
    solver = np.append(solver, QNS(model, method='rolia'))
    solver = np.append(solver, QNS(model, method='zhou'))
else:
    print("QNS solvers skipped: external qnsolver not available")

# MVA solvers - pass multiserver config via keyword args
solver = np.append(solver, MVA(model, config={'multiserver': 'softmin'}))
solver = np.append(solver, MVA(model, config={'multiserver': 'seidmann'}))
solver = np.append(solver, NC(model))

AvgTable = np.empty(len(solver), dtype=object)
for s in range(len(solver)):
    print(f'\nSOLVER: {solver[s].get_name()}')
    try:
        AvgTable[s] = solver[s].avg_table()
        print(AvgTable[s])
    except Exception as e:
        print(f'  Error: {e}')
        AvgTable[s] = None