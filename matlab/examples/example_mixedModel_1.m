clear node jobclass solver AvgTable

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
node{3} = Source(model,'Source');
node{4} = Sink(model,'Sink');

jobclass{1} = ClosedClass(model, 'ClosedClass', 2, node{1}, 0);
jobclass{2} = OpenClass(model, 'OpenClass', 0);

node{1}.setService(jobclass{1}, Erlang(3,2));
node{1}.setService(jobclass{2}, HyperExp(0.5,3.0,10.0));

node{2}.setService(jobclass{1}, HyperExp(0.1,1.0,10.0));
node{2}.setService(jobclass{2}, Exp(1));

node{3}.setArrival(jobclass{2}, Exp(0.1));

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = model.initRoutingMatrix;
P{1,1} = zeros(4); P{1,1}(1:2,1:2) = circul(2);
P{1,2} = zeros(4);
P{2,2} = [0,1,0,0; 0,0,0,1; 1,0,0,0; 0,0,0,0];
P{2,1} = zeros(4);

model.link(P);
%%
options = Solver.defaultOptions;
options.keep=true;
options.verbose=1;
options.cutoff = 3;
options.seed = 23000;
%options.samples=2e4;

disp('This example shows the execution of the solver on a 2-class 2-node mixed model.')
% This part illustrates the execution of different solvers
solver={};
solver{end+1} = SolverCTMC(model,options); % CTMC is infinite on this model
solver{end+1} = SolverJMT(model,options);
solver{end+1} = SolverSSA(model,options);
%solver{end+1} = SolverFluid(model,options);
solver{end+1} = SolverMVA(model,options);
%solver{end+1} = SolverNC(model,options);
solver{end+1} = SolverMAM(model,options);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
