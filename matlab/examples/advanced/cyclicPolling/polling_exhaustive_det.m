model = Network('M[2]/M[2]/1-Gated');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.POLLING);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass1 = OpenClass(model, 'myClass1');
source.setArrival(oclass1, Det.fitMean(1.0));
queue.setService(oclass1, Det.fitMean(0.001));

oclass2 = OpenClass(model, 'myClass2');
source.setArrival(oclass2, Det.fitMean(1.0));
queue.setService(oclass2, Det.fitMean(0.001));

queue.setPollingType(PollingType.EXHAUSTIVE)
queue.setSwitchover(oclass2, Immediate())
queue.setSwitchover(oclass1, Immediate())
%% Block 3: topology
P = model.initRoutingMatrix;
P{1} = Network.serialRouting(source,queue,sink);
P{2} = Network.serialRouting(source,queue,sink);
model.link(P);

MVA(model).getAvgTable() % solution is approximate in general
JMT(model,'seed',23000,'samples',1e5).getAvgTable()