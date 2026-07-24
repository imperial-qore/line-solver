clear P T R solver AvgTable

model = Network("model");

P1 = Place(model, "P1");
P2 = Place(model, "P2");
P3 = Place(model, "P3");
P4 = Place(model, "P4");

T1 = Transition(model, "T1");
T2 = Transition(model, "T2");
T3 = Transition(model, "T3");
T4 = Transition(model, "T4");

jobclass = ClosedClass(model, "Class1", 2, P1, 0);

mode1 = T1.addMode("Mode1");
T1.setDistribution(mode1, Exp(2));
T1.setEnablingConditions(mode1, jobclass, P1, 2);
T1.setFiringOutcome(mode1, jobclass, P2, 2);

mode2 = T2.addMode("Mode2");
T2.setDistribution(mode2, Erlang(3, 4));
T2.setEnablingConditions(mode2, jobclass, P2, 2);
T2.setFiringOutcome(mode2, jobclass, P3, 2);

mode3 = T3.addMode("Mode3");
T3.setDistribution(mode3, HyperExp(0.7, 3, 1.5));
T3.setEnablingConditions(mode3, jobclass, P3, 2);
T3.setFiringOutcome(mode3, jobclass, P4, 2);

mu0 = [1.0, 2.0];
phi0 = [0.6,1.0];

mode4 = T4.addMode("Mode4");
T4.setDistribution(mode4, Coxian(mu0, phi0));
T4.setEnablingConditions(mode4, jobclass, P4, 2);
T4.setFiringOutcome(mode4, jobclass, P1, 2);

routingMatrix = model.initRoutingMatrix();
routingMatrix.set(jobclass, jobclass, P1, T1, 1.0);
routingMatrix.set(jobclass, jobclass, P2, T2, 1.0);
routingMatrix.set(jobclass, jobclass, P3, T3, 1.0);
routingMatrix.set(jobclass, jobclass, P4, T4, 1.0);

routingMatrix.set(jobclass, jobclass, T1, P2, 1.0);
routingMatrix.set(jobclass, jobclass, T2, P3, 1.0);
routingMatrix.set(jobclass, jobclass, T3, P4, 1.0);
routingMatrix.set(jobclass, jobclass, T4, P1, 1.0);

model.link(routingMatrix);
%% Set Initial State    
P1.setState(jobclass.population);
P2.setState(0);
P3.setState(0);
P4.setState(0);

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
