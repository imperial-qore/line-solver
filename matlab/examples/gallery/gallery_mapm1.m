function model = gallery_mapm1(map)
if nargin < 1
    map = MAP.rand;
end
model = Network('MAP/M/1');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, map);
queue.setService(oclass, Exp(2));
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
end