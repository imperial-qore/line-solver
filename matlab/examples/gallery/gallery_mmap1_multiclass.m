function model = gallery_mmap1_multiclass(map1,map2)
if nargin < 1
    map1 = MAP.rand(2).setMean(0.5);
    map2 = MAP.rand(3).setMean(0.5);
end
model = Network('M/MAP/1');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass1 = OpenClass(model, 'myClass1');
source.setArrival(oclass1, Exp(0.35/map1.getMean));
queue.setService(oclass1, map1);
oclass2 = OpenClass(model, 'myClass2');
source.setArrival(oclass2, Exp(0.15/map2.getMean));
queue.setService(oclass2, map2);
%% Block 3: topology
P = model.initRoutingMatrix; 
P{1} = Network.serialRouting(source,queue,sink); 
P{2} = Network.serialRouting(source,queue,sink);
model.link(P);
end