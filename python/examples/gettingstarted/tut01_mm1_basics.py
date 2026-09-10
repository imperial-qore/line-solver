# %%
from line_solver import *
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
model = Network('M/M/1')
# Block 1: nodes
source = Source(model, 'Source')
queue = Queue(model, 'Queue', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')
# Block 2: classes
jobclass = OpenClass(model, 'Class1')
source.set_arrival(jobclass, Exp(1))
queue.set_service(jobclass, Exp(2))
# Block 3: topology
model.link(Network.serial_routing(source, queue, sink))
# %%
# Block 4: solution
AvgTable = JMT(model, seed=23000, samples=10000).avg_table()
# %%
# select a particular table row
print(tget(AvgTable, queue, jobclass))
# %%
# select a particular table row by node and class label
print(tget(AvgTable, 'Queue', 'Class1'))
# %%
print(AvgTable['RespT'].tolist())
# %%
print(tget(AvgTable,'Queue','Class1')['RespT'].tolist())