model = Network('M[2]/M[2]/1-Exhaustive');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.POLLING);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass1 = OpenClass(model, 'myClass1');
source.setArrival(oclass1, Exp(1.0));
queue.setService(oclass1, Exp(4.0));

oclass2 = OpenClass(model, 'myClass2');
source.setArrival(oclass2, Exp(0.8));
queue.setService(oclass2, Exp(1.5));

queue.setPollingType(PollingType.GATED)
queue.setSwitchover(oclass1, Exp(1.0))
queue.setSwitchover(oclass2, Exp(0.5))
%% Block 3: topology
P = model.initRoutingMatrix;
P{1} = Network.serialRouting(source,queue,sink);
P{2} = Network.serialRouting(source,queue,sink);
model.link(P);

MVA(model).getAvgTable() % solution is exact
JMT(model,'seed',23000,'samples',1e6).getAvgTable()