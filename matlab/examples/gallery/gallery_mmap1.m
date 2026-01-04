function model = gallery_mmap1(map)
if nargin < 1
    map = MAP.rand.setMean(0.5);
end
model = Network('M/MAP/1');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, Exp(1));
queue.setService(oclass, map);
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
end