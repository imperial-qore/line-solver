function model = gallery_mapmk(map,k)
if nargin < 1
    map = MAP.rand;
end
if nargin<2
    k = 2;
end
model = Network('MAP/M/k');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, map);
queue.setService(oclass, Exp(2));
queue.setNumberOfServers(k);
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
end