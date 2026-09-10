clear node jobclass solver AvgTable
model = Network('model');

node{1} = Delay(model, 'InfiniteServer');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
node{2}.setNumServers(2);

jobclass{1} = ClosedClass(model, 'Class1', 3, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 1, node{1}, 0);

node{1}.setService(jobclass{1}, Exp(3));
node{1}.setService(jobclass{2}, Exp(0.5));

node{2}.setService(jobclass{1}, Exp(0.1));
node{2}.setService(jobclass{2}, Exp(1));

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = cell(K,K);
P{1,1} = [0.3,0.1; 0.2,0];
P{1,2} = [0.6,0; 0.8,0];
P{2,2} = [0,1; 0,0];
P{2,1} = [0,0; 1,0];

model.link(P);
model.initFromMarginalAndStarted([2,1;1,0],[0,0;1,0]);

options = Solver.defaultOptions;
options.verbose = 1;
options.seed = 23000;
options.samples = 1e5;

disp('This example shows the execution of the solver on a 2-class 2-node class-switching model with specified initial state.')
% This part illustrates the execution of different solvers
solver{1} = CTMC(model, 'exact',options);
solver{end+1} = JMT(model,options);
solver{end+1} = SSA(model,options);
solver{end+1} = FLD(model,options);
% Solver.defaultOptions carries the GENERIC iter_tol (1e-4), which overrides
% SolverMVA's own default (1e-6) and stops the AMVA fixed point early: on this
% suite that costs up to 2.3e-4 of relative accuracy and makes the MATLAB row
% disagree with the Python row, which passes no options and keeps 1e-6.
mvaOptions = options;
mvaOptions.iter_tol = SolverMVA.defaultOptions.iter_tol;
solver{end+1} = MVA(model,mvaOptions);
solver{end+1} = NC(model,options);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',strrep(solver{s}.getName(),'Solver',''));
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end