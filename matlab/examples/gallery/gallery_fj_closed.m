function model = gallery_fj_closed()
% GALLERY_FJ_CLOSED Closed fork-join network (single class, two parallel tasks)
model = Network('Fork-Join-Closed');
%% Block 1: nodes
delay = Delay(model, 'Delay');
queue1 = Queue(model, 'Queue1', SchedStrategy.PS);
queue2 = Queue(model, 'Queue2', SchedStrategy.PS);
fork = Fork(model, 'Fork');
join = Join(model, 'Join', fork);
%% Block 2: classes
oclass = ClosedClass(model, 'class1', 5, delay);
delay.setService(oclass, Exp(1.0));
queue1.setService(oclass, Exp(1.0));
queue2.setService(oclass, Exp(1.0));
%% Block 3: topology
P = model.initRoutingMatrix();
P{oclass,oclass}(delay,fork) = 1.0;
P{oclass,oclass}(fork,queue1) = 1.0;
P{oclass,oclass}(fork,queue2) = 1.0;
P{oclass,oclass}(queue1,join) = 1.0;
P{oclass,oclass}(queue2,join) = 1.0;
P{oclass,oclass}(join,delay) = 1.0;
model.link(P);
end
