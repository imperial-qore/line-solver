function model = gallery_fj_open()
% GALLERY_FJ_OPEN Open fork-join network (single class, two parallel tasks)
model = Network('Fork-Join-Open');
%% Block 1: nodes
source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
fork = Fork(model, 'Fork');
join = Join(model, 'Join', fork);
sink = Sink(model, 'Sink');
%% Block 2: classes
oclass = OpenClass(model, 'class1');
source.setArrival(oclass, Exp(0.05));
queue1.setService(oclass, Exp(1.0));
queue2.setService(oclass, Exp(2.0));
%% Block 3: topology
P = model.initRoutingMatrix();
P{oclass,oclass}(source,fork) = 1.0;
P{oclass,oclass}(fork,queue1) = 1.0;
P{oclass,oclass}(fork,queue2) = 1.0;
P{oclass,oclass}(queue1,join) = 1.0;
P{oclass,oclass}(queue2,join) = 1.0;
P{oclass,oclass}(join,sink) = 1.0;
model.link(P);
end
