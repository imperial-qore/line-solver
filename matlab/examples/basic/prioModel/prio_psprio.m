clear node jobclass solver AvgTable;

model = Network('MyNetwork');

node{1} = Delay(model, 'SlowDelay');
node{2} = Queue(model, 'PSPRIOQueue', SchedStrategy.PSPRIO);

jobclass{1} = ClosedClass(model, 'Class1', 2, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 2, node{1}, 1);

node{1}.setService(jobclass{1}, Erlang(3,2));
node{1}.setService(jobclass{2}, HyperExp(0.5,3.0,10.0));

node{2}.setService(jobclass{1}, HyperExp(0.1,1.0,10.0));
node{2}.setService(jobclass{2}, Exp(1));

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = model.initRoutingMatrix;
P{1} = Network.serialRouting(node);
P{2} = Network.serialRouting(node);
model.link(P);
%%
% This part illustrates the execution of different solvers
solver = {};
solver{end+1} = CTMC(model);
solver{end+1} = JMT(model,'seed',23000,'verbose',true,'samples',5e3);
solver{end+1} = SSA(model,'seed',23000,'verbose',true,'samples',5e3);
% solver{end+1} = FLD(model);
% solver{end+1} = MVA(model,'exact');
% solver{end+1} = NC(model,'exact');
% solver{end+1} = MAM(model);
% solver{end+1} = LINE(model);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',strrep(solver{s}.getName(),'Solver',''));
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
