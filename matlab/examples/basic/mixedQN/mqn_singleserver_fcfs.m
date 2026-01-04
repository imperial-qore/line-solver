clear node jobclass solver AvgTable

% this is a difficult sparse LDMX model
model = Network('model');

M = 4;
node{1} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{2} = Queue(model, 'Queue2', SchedStrategy.FCFS);
node{3} = Queue(model, 'Queue3', SchedStrategy.FCFS);
node{4} = Queue(model, 'Queue4', SchedStrategy.FCFS); % only closed classes
source = Source(model,'Source');
sink = Sink(model,'Sink');

jobclass{1} = ClosedClass(model, 'ClosedClass', 100, node{1}, 0);
jobclass{2} = OpenClass(model, 'OpenClass', 0);

for i=1:M
    node{i}.setService(jobclass{1}, Exp(i));
    node{i}.setService(jobclass{2}, Exp(sqrt(i)));
end

source.setArrival(jobclass{2}, APH.fitMeanAndSCV(3,64));

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = model.initRoutingMatrix;
P{1,1} = zeros(M+1);
P{1,1} = Network.serialRouting(node{1},node{2},node{3},node{4});
P{1,2} = zeros(M+1);
P{2,1} = zeros(M+1);
P{2,2} = Network.serialRouting(source,node{1},node{2},node{3},sink);

model.link(P);
%%
options = Solver.defaultOptions;
options.keep=false;
options.verbose=1;
options.cutoff = 3;
options.seed = 23000;
options.samples = 2e4;
optionssa = options; optionssa.cutoff = Inf;

disp('This example shows the execution of the solver on a 2-class mixed model with 4 single server nodes.')
% This part illustrates the execution of different solvers
solver={};
%solver{end+1} = CTMC(model,options); % CTMC is infinite on this model
solver{end+1} = JMT(model,options);
solver{end+1} = FLD(model,options);
solver{end+1} = MVA(model,'lin');
solver{end+1} = NC(model,options);
solver{end+1} = DES(model,options);
%solver{end+1} = MAM(model,options);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}    
end
