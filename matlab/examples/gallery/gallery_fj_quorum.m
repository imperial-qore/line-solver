function model = gallery_fj_quorum()
% GALLERY_FJ_QUORUM Closed fork-join network with a 2-of-3 quorum join
%
% The join fires on the SECOND of the three sibling tasks; the third is
% discarded when it arrives. Solved exactly by SolverLDES and SolverJMT;
% SolverMVA and SolverNC charge the second order statistic of the branch
% completion times at the join (fj_ordstat_exp).
model = Network('Fork-Join-Quorum');
%% Block 1: nodes
delay = Delay(model, 'Delay');
queue1 = Queue(model, 'Queue1', SchedStrategy.PS);
queue2 = Queue(model, 'Queue2', SchedStrategy.PS);
queue3 = Queue(model, 'Queue3', SchedStrategy.PS);
fork = Fork(model, 'Fork');
join = Join(model, 'Join', fork);
%% Block 2: classes
oclass = ClosedClass(model, 'class1', 5, delay);
delay.setService(oclass, Exp(1.0));
queue1.setService(oclass, Exp(2.0));
queue2.setService(oclass, Exp(2.0));
queue3.setService(oclass, Exp(2.0));
join.setStrategy(oclass, JoinStrategy.PARTIAL);
join.setRequired(oclass, 2);
%% Block 3: topology
P = model.initRoutingMatrix();
P{oclass,oclass}(delay,fork) = 1.0;
P{oclass,oclass}(fork,queue1) = 1.0;
P{oclass,oclass}(fork,queue2) = 1.0;
P{oclass,oclass}(fork,queue3) = 1.0;
P{oclass,oclass}(queue1,join) = 1.0;
P{oclass,oclass}(queue2,join) = 1.0;
P{oclass,oclass}(queue3,join) = 1.0;
P{oclass,oclass}(join,delay) = 1.0;
model.link(P);
end
