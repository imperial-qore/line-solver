clear node jobclass solver AvgTable;
% comparison of product-form scheduling policies
c=1;
psmodel = Network('PS scheduling model');

node{1} = Delay(psmodel, 'Delay');
node{2} = Queue(psmodel, 'Queue1', SchedStrategy.PS);
node{2}.setNumberOfServers(c);

jobclass{1} = ClosedClass(psmodel, 'Class1', 2, node{1}, 0);
jobclass{2} = ClosedClass(psmodel, 'Class2', 2, node{1}, 0);

node{1}.setService(jobclass{1}, Erlang(3,2));
node{1}.setService(jobclass{2}, HyperExp(0.5,3.0,10.0));

node{2}.setService(jobclass{1}, Exp(1));
node{2}.setService(jobclass{2}, Exp(1));

P = psmodel.initRoutingMatrix;
P{1} = Network.serialRouting(node);
P{2} = Network.serialRouting(node);
psmodel.link(P);
%%
fcfsmodel = Network('FCFS scheduling model');

node{1} = Delay(fcfsmodel, 'Delay');
node{2} = Queue(fcfsmodel, 'Queue1', SchedStrategy.FCFS);
node{2}.setNumberOfServers(c);

jobclass{1} = ClosedClass(fcfsmodel, 'Class1', 2, node{1}, 0);
jobclass{2} = ClosedClass(fcfsmodel, 'Class2', 2, node{1}, 0);

node{1}.setService(jobclass{1}, Erlang(3,2));
node{1}.setService(jobclass{2}, HyperExp(0.5,3.0,10.0));

node{2}.setService(jobclass{1}, Exp(1));
node{2}.setService(jobclass{2}, Exp(1));

P = fcfsmodel.initRoutingMatrix;
P{1} = Network.serialRouting(node);
P{2} = Network.serialRouting(node);
fcfsmodel.link(P);
%%
lcfsprmodel = Network('LCFS-PR scheduling model');

node{1} = Delay(lcfsprmodel, 'Delay');
node{2} = Queue(lcfsprmodel, 'Queue1', SchedStrategy.LCFSPR);
node{2}.setNumberOfServers(c);

jobclass{1} = ClosedClass(lcfsprmodel, 'Class1', 2, node{1}, 0);
jobclass{2} = ClosedClass(lcfsprmodel, 'Class2', 2, node{1}, 0);

node{1}.setService(jobclass{1}, Erlang(3,2));
node{1}.setService(jobclass{2}, HyperExp(0.5,3.0,10.0));

node{2}.setService(jobclass{1}, Exp(1));
node{2}.setService(jobclass{2}, Exp(1));

P = lcfsprmodel.initRoutingMatrix;
P{1} = Network.serialRouting(node);
P{2} = Network.serialRouting(node);
lcfsprmodel.link(P);
%%
% This part illustrates the execution of different solvers
solver = {};
solver{end+1} = CTMC(psmodel);
solver{end+1} = CTMC(fcfsmodel);
solver{end+1} = CTMC(lcfsprmodel);
for s=1:length(solver)
    fprintf(1,'MODEL: %s\n',solver{s}.model.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end

