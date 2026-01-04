model = Network('M[2]/M[2]/1-Gated');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass1 = OpenClass(model, 'myClass1');
source.setArrival(oclass1, Exp(0.2));
queue.setService(oclass1, Exp(0.1));

oclass2 = OpenClass(model, 'myClass2');
source.setArrival(oclass2, Exp(0.8));
queue.setService(oclass2, Exp(1.5));

queue.setSwitchover(oclass1, oclass2, Exp(1))
queue.setSwitchover(oclass2, oclass1, Erlang(1,2))
%% Block 3: topology
P = model.initRoutingMatrix;
P{1} = Network.serialRouting(source,queue,sink);
P{2} = Network.serialRouting(source,queue,sink);
model.link(P);

%MVA(model).getAvgTable()
JMT(model,'seed',23000,'keep',true).getAvgTable()