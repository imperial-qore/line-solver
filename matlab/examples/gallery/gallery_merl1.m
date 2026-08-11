function [model,source,queue,sink,oclass]=gallery_me1
model = Network('M/E/1');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, Exp(1));
queue.setService(oclass, Erlang.fitMeanAndOrder(0.5,2));
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
end