function model = gallery_mm1k(K)
% GALLERY_MM1K M/M/1/K queue with finite capacity K (blocking / loss)
if nargin < 1
    K = 3;
end
model = Network('M/M/1/K');
%% Block 1: nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
queue.setNumberOfServers(1);
queue.setCapacity(K);
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
end
