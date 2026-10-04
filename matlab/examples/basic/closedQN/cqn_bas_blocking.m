clear node jobclass;
% Closed network with Blocking-After-Service (BAS) finite-buffer blocking.
% Queue1 uses the BAS drop rule; Queue2 has a finite buffer (capacity 1).
% The 'sqd' (Smith Queue Decomposition) approximation handles this model;
% MVA also auto-routes BAS models to 'sqd' under its default method.
model = Network('cqn_bas_blocking');

node{1} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{2} = Queue(model, 'Queue2', SchedStrategy.FCFS);

jobclass{1} = ClosedClass(model, 'Class1', 2, node{1}, 0);

node{1}.setService(jobclass{1}, Exp(1.0));
node{2}.setService(jobclass{1}, Exp(0.8));

node{2}.setCap(1);
node{1}.setDropRule(jobclass{1}, DropStrategy.BAS);

model.link(Network.serialRouting(node{1}, node{2}));

% Run solver with the SQD method explicitly
solver = MVA(model, 'sqd');
avgTable = solver.getAvgTable();
fprintf('SOLVER: MVA (method=sqd)\n');
disp(avgTable);
