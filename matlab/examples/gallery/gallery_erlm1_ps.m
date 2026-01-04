function model = gallery_erlm1_ps
model = Network('Er/M/1-PS');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.PS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, Erlang.fitMeanAndOrder(1,5));
queue.setService(oclass, Exp(2));
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
end
