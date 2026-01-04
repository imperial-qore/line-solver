function model = gallery_coxm1
model = Network('Cox/M/1');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, Coxian.fitCentral(1,0.99,1.999));
queue.setService(oclass, Exp(2));
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
end