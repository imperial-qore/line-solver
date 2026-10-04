clear node jobclass solver AvgTable

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{3} = Source(model,'Source');
node{4} = Sink(model,'Sink');

jobclass{1} = OpenClass(model, 'Class1', 0);

node{1}.setService(jobclass{1}, HyperExp(0.5,3.0,10.0));
node{2}.setService(jobclass{1}, Exp(1));
node{3}.setArrival(jobclass{1}, Exp(0.1));

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = cell(K,K);
P{1,1} = [0,1,0,0; 0,0,0,1; 1,0,0,0; 0,0,0,0];

model.link(P);
%%
options = Solver.defaultOptions;
options.keep=true;
options.verbose=1;
options.cutoff = 10;
options.seed = 23000;
options.iter_max = 200;
%options.samples=2e4;

disp('This example shows the execution of the solver on a 1-class 2-node open model.')
% This part illustrates the execution of different solvers
solver={};
solver{end+1} = CTMC(model, 'exact',options); % CTMC is infinite on this model
solver{end+1} = FLD(model,options);
% Solver.defaultOptions carries the GENERIC iter_tol (1e-4), which overrides
% SolverMVA's own default (1e-6) and stops the AMVA fixed point early: on this
% suite that costs up to 2.3e-4 of relative accuracy and makes the MATLAB row
% disagree with the Python row, which passes no options and keeps 1e-6.
mvaOptions = options;
mvaOptions.iter_tol = SolverMVA.defaultOptions.iter_tol;
solver{end+1} = MVA(model,mvaOptions);
solver{end+1} = MAM(model,options);
solver{end+1} = NC(model,options);
solver{end+1} = JMT(model,options);
% SSA needs a larger budget than the other solvers here: at 1e4 events the
% standard error of Util/Tput on this rho=0.1 queue is ~4%, which is the same
% size as the cross-codebase parity tolerance, so no two RNG streams can be
% expected to agree. 5e4 brings it to ~1.4%. Only SSA is raised, so the JMT and
% LDES sample budgets (and their baselines) are unchanged.
ssaOptions = options;
ssaOptions.samples = 5e4;
solver{end+1} = SSA(model,ssaOptions);
% Solver.defaultOptions carries the GENERIC samples budget (1e4), which
% overrides SolverLDES's own default (2e5) and leaves the simulation 20x
% shorter than the solver asks for. Restore the solver's own default.
ldesOptions = options;
ldesOptions.samples = SolverLDES.defaultOptions.samples;
solver{end+1} = LDES(model,ldesOptions);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',strrep(solver{s}.getName(),'Solver',''));
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
