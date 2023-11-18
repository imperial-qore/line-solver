clear node jobclass

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
node{3} = Queue(model, 'Queue2', SchedStrategy.DPS);

jobclass{1} = ClosedClass(model, 'Class1', 2, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 1, node{1}, 0);

model.addLink(node{1}, node{2});
model.addLink(node{1}, node{3});
model.addLink(node{2}, node{1});
model.addLink(node{3}, node{1});

node{1}.setService(jobclass{1}, Exp(3));
node{1}.setService(jobclass{2}, Exp(0.5));

% this is ignored since the node is PS
w1 = 5; node{2}.setService(jobclass{1}, Exp(0.1), w1);
w2 = 1; node{2}.setService(jobclass{2}, Exp(1), w2);

% this is not ignored since the node is DPS
w1 = 1; node{3}.setService(jobclass{1}, Exp(0.1), w1);
w2 = 5; node{3}.setService(jobclass{2}, Exp(1), w2);

node{1}.setProbRouting(jobclass{1}, node{2}, 0.3)
node{1}.setProbRouting(jobclass{1}, node{3}, 0.7)
node{1}.setProbRouting(jobclass{2}, node{2}, 0.7)
node{1}.setProbRouting(jobclass{2}, node{3}, 0.3)
node{2}.setProbRouting(jobclass{1}, node{1}, 1.0)
node{2}.setProbRouting(jobclass{2}, node{1}, 1.0)
node{3}.setProbRouting(jobclass{1}, node{1}, 1.0)
node{3}.setProbRouting(jobclass{2}, node{1}, 1.0)

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

% This part illustrates the execution of different solvers

solver={};
options = Solver.defaultOptions;
options.verbose=1;
options.samples=1e4;
options.seed = 23000;
solver{end+1} = SolverCTMC(model,options);
solver{end+1} = SolverJMT(model,options);
%solver{end+1} = SolverSSA(model,options);
solver{end+1} = SolverFluid(model,options);
solver{end+1} = SolverMVA(model,options);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end