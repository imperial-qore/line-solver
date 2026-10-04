# %%
from line_solver import *
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
model = Network('MRP')
# Block 1: nodes
delay = Delay(model, 'WorkingState')
queue = Queue(model, 'RepairQueue', SchedStrategy.FCFS)
queue.set_number_of_servers(2)
# Block 2: classes
cclass = ClosedClass(model, 'Machines', 3, delay)
delay.set_service(cclass, Exp(0.5))
queue.set_service(cclass, Exp(4.0))
# Block 3: topology
model.link(Network.serial_routing(delay, queue))
# %%
# Block 4: solution
solver = CTMC(model)
ctmcAvgTable = solver.avg_table()
print(ctmcAvgTable)
# %%
stateSpace, nodeStateSpace = solver.state_space()
print("CTMC state space:")
print(stateSpace)
print(nodeStateSpace)
# %%
infGen, eventFilt = solver.generator()
print("CTMC infinitesimal generator:")
print(infGen)
# %%
CTMC.print_inf_gen(infGen, stateSpace)