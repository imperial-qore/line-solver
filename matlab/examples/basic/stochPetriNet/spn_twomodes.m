clear P T R solver AvgTable

model = Network('model');

P1 = Place(model, 'P1');
P2 = Place(model, 'P2');
T1 = Transition(model, 'T1');
T2 = Transition(model, 'T2');

jobclass = ClosedClass(model, 'Class1', 10, P1, 0); % automatically added to P1 at initialization

% T1
mode = T1.addMode('Mode1');
%T1.setNumberOfServers(mode,Inf)
T1.setDistribution(mode, Exp(2));
T1.setEnablingConditions(mode, jobclass, P1, 4);
T1.setFiringOutcome(mode, jobclass, P2, 4);

% T2
mode = T2.addMode('Mode2');
%T2.setNumberOfServers(mode,Inf)
T2.setDistribution(mode, Exp(3));
T2.setEnablingConditions(mode, jobclass, P2, 2);
T2.setFiringOutcome(mode, jobclass, P1, 2);

routingMatrix = model.initRoutingMatrix();
routingMatrix.set(jobclass, jobclass, P1, T1, 1.0);
routingMatrix.set(jobclass, jobclass, P2, T2, 1.0);
routingMatrix.set(jobclass, jobclass, T1, P2, 1.0);
routingMatrix.set(jobclass, jobclass, T2, P1, 1.0);

model.link(routingMatrix);

%% Solver
options = Solver.defaultOptions;
options.keep=2;
options.verbose=1;
options.cutoff = 10;
options.seed = 23000;

solver = {};
solver{1} = JMT(model,options);
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}
