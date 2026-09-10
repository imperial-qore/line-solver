# %%
from line_solver import *
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
model = Network('RL')
queue = Queue(model, 'Queue', SchedStrategy.FCFS)
K = 3
N = (1, 0, 0)
jobclass = []
for k in range(K):
    jobclass.append(ClosedClass(model, 'Class' + str(k+1), N[k], queue))
    queue.set_service(jobclass[k], Erlang.fit_mean_and_order(1+k, 2))
P = model.init_routing_matrix()
P.set(jobclass[0], jobclass[1], queue, queue, 1.0)
P.set(jobclass[1], jobclass[2], queue, queue, 1.0)
P.set(jobclass[2], jobclass[0], queue, queue, 1.0)
model.link(P)
sn = model.struct()
# %%
ncAvgTable = NC(model).avg_table()
print(ncAvgTable)
# %%
ncAvgSysTable = NC(model).avg_sys_table()
print(ncAvgSysTable)
# %%
jobclass[0].completes = False
jobclass[1].completes = False
ncAvgSysTable2 = NC(model).avg_sys_table()
print(ncAvgSysTable2)