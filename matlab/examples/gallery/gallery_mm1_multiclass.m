function model = gallery_mm1_multiclass
model = Network('M[2]/M[2]/1');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass1 = OpenClass(model, 'myClass1');
source.setArrival(oclass1, Exp(1));
queue.setService(oclass1, Exp(4));

oclass2 = OpenClass(model, 'myClass2');
source.setArrival(oclass2, Exp(0.5));
queue.setService(oclass2, Exp(4));
%% Block 3: topology
P = model.initRoutingMatrix;
P{1} = Network.serialRouting(source,queue,sink);
P{2} = Network.serialRouting(source,queue,sink);
model.link(P);
end