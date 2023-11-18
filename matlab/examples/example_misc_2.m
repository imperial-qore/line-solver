clear node jobclass

% This example illustrates the automatic checks of the model features
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
node{3} = Queue(model, 'Queue2', SchedStrategy.PS);
node{4} = Queue(model, 'Queue3', SchedStrategy.PS);
node{5} = Queue(model, 'Queue4', SchedStrategy.PS);

jobclass{1} = ClosedClass(model, 'Class1', 1, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 2, node{1}, 0);

M = model.getNumberOfNodes;
K = model.getNumberOfClasses;

for i=1:M
    node{i}.setService(jobclass{1}, Exp(3));
    node{i}.setService(jobclass{2}, Exp(10.0));
end

P = model.initRoutingMatrix;
for r=1:K P{r} = circul(M); end
model.link(P);

fprintf(1,'This example shows how to update a model parameters and its solution.\n');

%pause
options = Solver.defaultOptions;
options.verbose = 0;
options.seed = 23000;
solver = SolverCTMC(model,options);

I=10;

T0=tic;
for it=1:I
    QN0 = solver.getAvg;
end
fprintf(1, 'Time for solution without model updates: %f sec\n',toc(T0));

T1=tic;
for it=1:I
    node{2}.setService(jobclass{1}, Exp(it));
    model.refreshStruct % this forces to regenerate the internal data structures
    QN1 = solver.getAvg;
end
fprintf(1, 'Time for solution with updates and full model refreshes: %f sec\n',toc(T1));

T2=tic;
for it=1:I
    node{2}.setService(jobclass{1}, Exp(it));
    model.refreshService; % this regenerates the internal data structures pertaining to the service processes
    QN2 = solver.getAvg;
end
fprintf(1, 'Time for solution with updates and partial model refreshes: %f sec\n',toc(T2));
