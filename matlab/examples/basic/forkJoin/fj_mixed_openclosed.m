clear solver AvgTable;

% Flat counterpart of lqn_fork_open_arrival.m: one fork-join traversed by an
% open class (Poisson rate 0.1) and by a closed class (1 job, think time 1), so
% the same fork carries exogenous and closed traffic at once. The MVA fork-join
% transform solves this; the layered builder refuses the analogous LQN because
% the transform's auxiliary Source collides with the layer's open-stream Source.
%
% Reference (SolverJMT, seed 23000, 2e5 samples): Branch1 open RespT 0.362,
% Branch1 closed RespT 0.310, Join open RespT 0.215, closed tput 0.540.

model = Network('ForkJoinOpenClosed');

source = Source(model,'Source');
client = Delay(model,'Client');
prefork = Queue(model,'PreFork',SchedStrategy.FCFS);
branch1 = Queue(model,'Branch1',SchedStrategy.FCFS);
branch2 = Queue(model,'Branch2',SchedStrategy.FCFS);
postjoin = Queue(model,'PostJoin',SchedStrategy.FCFS);
fork = Fork(model,'Fork');
join = Join(model,'Join',fork);
sink = Sink(model,'Sink');

oclass = OpenClass(model, 'Open');
cclass = ClosedClass(model, 'Closed', 1, client);

source.setArrival(oclass, Exp(0.1));
client.setService(cclass, Exp.fitMean(1.0));
client.setService(oclass, Disabled());
prefork.setService(oclass, Exp.fitMean(0.2));
prefork.setService(cclass, Exp.fitMean(0.2));
branch1.setService(oclass, Exp.fitMean(0.3));
branch1.setService(cclass, Exp.fitMean(0.3));
branch2.setService(oclass, Exp.fitMean(0.4));
branch2.setService(cclass, Exp.fitMean(0.4));
postjoin.setService(oclass, Exp.fitMean(0.1));
postjoin.setService(cclass, Exp.fitMean(0.1));

P = model.initRoutingMatrix();
P{oclass,oclass}(source,prefork) = 1.0;
P{oclass,oclass}(prefork,fork) = 1.0;
P{oclass,oclass}(fork,branch1) = 1.0;
P{oclass,oclass}(fork,branch2) = 1.0;
P{oclass,oclass}(branch1,join) = 1.0;
P{oclass,oclass}(branch2,join) = 1.0;
P{oclass,oclass}(join,postjoin) = 1.0;
P{oclass,oclass}(postjoin,sink) = 1.0;

P{cclass,cclass}(client,prefork) = 1.0;
P{cclass,cclass}(prefork,fork) = 1.0;
P{cclass,cclass}(fork,branch1) = 1.0;
P{cclass,cclass}(fork,branch2) = 1.0;
P{cclass,cclass}(branch1,join) = 1.0;
P{cclass,cclass}(branch2,join) = 1.0;
P{cclass,cclass}(join,postjoin) = 1.0;
P{cclass,cclass}(postjoin,client) = 1.0;

model.link(P);

solver = {};
solver{end+1} = JMT(model,'seed',23000,'samples',2e5);
solver{end+1} = MVA(model);

AvgTable = {};
for s=1:length(solver)
    AvgTable{end+1} = solver{s}.getAvgTable;
    AvgTable{s}
end
