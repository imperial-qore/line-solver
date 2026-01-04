function model = gallery_mhyp1()
model = Network('M/H/1');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, Exp(1));
queue.setService(oclass, HyperExp.fitMeanAndSCV(0.5,4));
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
end
