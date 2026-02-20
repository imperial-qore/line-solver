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
solver{end+1} = CTMC(model,options); % CTMC is infinite on this model
solver{end+1} = FLD(model,options);
solver{end+1} = MVA(model,options);
solver{end+1} = MAM(model,options);
solver{end+1} = NC(model,options);
solver{end+1} = JMT(model,options);
solver{end+1} = SSA(model,options);
solver{end+1} = DES(model,options);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',strrep(solver{s}.getName(),'Solver',''));
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
