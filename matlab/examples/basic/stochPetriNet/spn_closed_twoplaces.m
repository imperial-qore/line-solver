clear P T R solver AvgTable

model = Network("model");

P1 = Place(model, "P1");
P2 = Place(model, "P2");

T1 = Transition(model, "T1");
T2 = Transition(model, "T2");
T3 = Transition(model, "T3");

jobclass1 = ClosedClass(model, "Class1", 10, P1, 0);
jobclass2 = ClosedClass(model, "Class2", 7, P1, 0);

% T1
mode1 = T1.addMode("Mode1");
%T1.setNumberOfServers(mode1,Inf);
T1.setDistribution(mode1, Exp(2));
T1.setEnablingConditions(mode1, jobclass1, P1, 2);
T1.setFiringOutcome(mode1, jobclass1, P2, 2);

% T1
mode2 = T1.addMode("Mode2");
%T1.setNumberOfServers(mode2,Inf);
T1.setDistribution(mode2, Exp(3));
T1.setEnablingConditions(mode2, jobclass2, P1, 1);
T1.setFiringOutcome(mode2, jobclass2, P2, 1);

% T2
mode3 = T2.addMode("Mode3");
%T2.setNumberOfServers(mode3,Inf);
T2.setDistribution(mode3, Erlang(1.5, 2));
T2.setEnablingConditions(mode3, jobclass1, P2, 1);
T2.setFiringOutcome(mode3, jobclass1, P1, 1);

% T3
mode4 = T3.addMode("Mode4");
%T3.setNumberOfServers(mode4,Inf);
T3.setDistribution(mode4, Exp(0.5));
T3.setEnablingConditions(mode4, jobclass2, P2, 4);
T3.setFiringOutcome(mode4, jobclass2, P1, 4);

routingMatrix = model.initRoutingMatrix();
routingMatrix.set(jobclass1, jobclass1, P1, T1, 1.0);
routingMatrix.set(jobclass2, jobclass2, P1, T1, 1.0);
routingMatrix.set(jobclass1, jobclass1, P2, T2, 1.0);
routingMatrix.set(jobclass2, jobclass2, P2, T3, 1.0);
routingMatrix.set(jobclass1, jobclass1, T1, P2, 1.0);
routingMatrix.set(jobclass2, jobclass2, T1, P2, 1.0);
routingMatrix.set(jobclass1, jobclass1, T2, P1, 1.0);
routingMatrix.set(jobclass2, jobclass2, T3, P1, 1.0);

model.link(routingMatrix);
%% Set Initial State    
P1.setState([jobclass1.population,jobclass2.population]);
P2.setState([0,0]);

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
