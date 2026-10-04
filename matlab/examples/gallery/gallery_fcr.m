function model = gallery_fcr(K)
% GALLERY_FCR Finite capacity region with dropping around a single queue
% (equivalent to an M/M/1/K queue where K is the region capacity)
if nargin < 1
    K = 3;
end
model = Network('FCR-Dropping');
%% Block 1: nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');
%% Block 2: classes
oclass = OpenClass(model, 'Class1', 0);
source.setArrival(oclass, Exp(0.8));
queue.setService(oclass, Exp(1.0));
%% Block 3: topology
P = model.initRoutingMatrix();
P.set(oclass, oclass, source, queue, 1.0);
P.set(oclass, oclass, queue, sink, 1.0);
model.link(P);
%% Block 4: finite capacity region with dropping
fcr = model.addRegion({queue});
fcr.setGlobalMaxJobs(K);
fcr.setDropRule(oclass, true);
end
