clear node jobclass solver AvgTable;

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);

jobclass{1} = ClosedClass(model, 'Class1', 2, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 2, node{1}, 0);

node{1}.setService(jobclass{1}, Erlang(3,2));
node{1}.setService(jobclass{2}, HyperExp(0.5,3.0,10.0));

node{2}.setService(jobclass{1}, HyperExp(0.1,1.0,10.0));
node{2}.setService(jobclass{2}, Exp(1));

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = model.initRoutingMatrix;
P{1,1} = [0.3,0.1; 0.2,0];
P{1,2} = [0.6,0.0; 0.8,0];
P{2,2} = [0,1; 0,0];
P{2,1} = [0,0; 1,0];

model.link(P);
%%
% This part illustrates the execution of different solvers
solver = {};
solver{end+1} = SolverCTMC(model);
solver{end+1} = SolverJMT(model,'seed',23000,'verbose',true,'samples',5e3);
solver{end+1} = SolverSSA(model,'seed',23000,'verbose',true,'samples',5e3);
solver{end+1} = SolverFluid(model);
solver{end+1} = SolverMVA(model,'exact');
solver{end+1} = SolverNC(model,'exact');
solver{end+1} = SolverMAM(model);
solver{end+1} = LINE(model);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
