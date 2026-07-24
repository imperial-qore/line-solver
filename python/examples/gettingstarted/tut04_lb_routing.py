# %%
from line_solver import *
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
model = Network('RRLB')
source = Source(model, 'Source')
lb = Router(model, 'LB')
queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
sink = Sink(model, 'Sink')
oclass = OpenClass(model, 'Class1')
source.set_arrival(oclass, Exp(1))
queue1.set_service(oclass, Exp(2))
queue2.set_service(oclass, Exp(2))
# %%
model.add_link(source, lb)
model.add_link(lb, queue1)
model.add_link(lb, queue2)
model.add_link(queue1, sink)
model.add_link(queue2, sink)
lb.set_routing(oclass, RoutingStrategy.RAND)
jmtAvgTable = JMT(model, seed=23000, samples=10000).avg_table()
print(jmtAvgTable)
# %%
model.reset()
lb.set_routing(oclass, RoutingStrategy.RROBIN)
jmtAvgTableRR = JMT(model, seed=23000, samples=10000).avg_table()
print(jmtAvgTableRR)