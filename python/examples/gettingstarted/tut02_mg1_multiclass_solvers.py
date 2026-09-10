# %%
from line_solver import *
import os
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
model = Network('M/G/1')
source = Source(model, 'Source')
queue = Queue(model, 'Queue', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')
jobclass1 = OpenClass(model, 'Class1')
jobclass2 = OpenClass(model, 'Class2')
source.set_arrival(jobclass1, Exp(0.5))
source.set_arrival(jobclass2, Exp(0.5))
queue.set_service(jobclass1, Erlang.fit_mean_and_scv(1.0, 1 / 3))
# First use raw Replayer for JMT (matches MATLAB)
tracePath = lineRootFolder() + "/examples/gettingstarted/example_trace.txt"
queue.set_service(jobclass2, Replayer(tracePath))
# %%
P = model.init_routing_matrix()
P.set(jobclass1, Network.serial_routing(source,queue,sink))
P.set(jobclass2, Network.serial_routing(source,queue,sink))
model.link(P)
# %%
jmtAvgTable = JMT(model, seed=23000, samples=10000, verbose=True).avg_table()
print(jmtAvgTable)
# %%
# Now switch to fitted APH for CTMC and MAM (matches MATLAB)
queue.set_service(jobclass2, Replayer(tracePath).fit_aph())
# %%
# Set config.nonmkv = 'none' to disable automatic non-Markovian conversion (matches MATLAB)
ctmcAvgTable2 = CTMC(model, cutoff=2, verbose=True, config={'nonmkv': 'none'}).avg_table()
print(ctmcAvgTable2)
# %%
ctmcAvgTable4 = CTMC(model, cutoff=4, verbose=True, config={'nonmkv': 'none'}).avg_table()
print(ctmcAvgTable4)
# %%
mamAvgTable = MAM(model, verbose=True, config={'nonmkv': 'none'}).avg_table()
print(mamAvgTable)
