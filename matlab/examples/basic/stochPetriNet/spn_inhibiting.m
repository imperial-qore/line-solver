clear P T R solver AvgTable

model = Network("model");

P1 = Place(model, "P1");
P2 = Place(model, "P2");
P3 = Place(model, "P3");

T1 = Transition(model, "T1");
T2 = Transition(model, "T2");
T3 = Transition(model, "T3");

jobclass = ClosedClass(model, "Class1", 4, P1, 0);

mode1 = T1.addMode("Mode1");
T1.setDistribution(mode1, Exp(2));
T1.setEnablingConditions(mode1, jobclass, P1, 2);
T1.setFiringOutcome(mode1, jobclass, P2, 2);

mode2 = T1.addMode("Mode2");
T1.setDistribution(mode2, Exp(1));
T1.setEnablingConditions(mode2, jobclass, P1, 1);
T1.setFiringOutcome(mode2, jobclass, P3, 1);

mode3 = T2.addMode("Mode3");
T2.setDistribution(mode3, Exp(4));
T2.setEnablingConditions(mode3, jobclass, P2, 1);
T2.setFiringOutcome(mode3, jobclass, P1, 1);

mode4 = T3.addMode("Mode4");
T3.setDistribution(mode4, Exp(1));
T3.setEnablingConditions(mode4, jobclass, P3, 3);
T3.setInhibitingConditions(mode4, jobclass, P2, 1);
T3.setFiringOutcome(mode4, jobclass, P1, 3);

routingMatrix = model.initRoutingMatrix();
routingMatrix.set(jobclass, jobclass, P1, T1, 1.0);
routingMatrix.set(jobclass, jobclass, P2, T2, 1.0);
routingMatrix.set(jobclass, jobclass, P2, T3, 1.0);
routingMatrix.set(jobclass, jobclass, P3, T3, 1.0);

routingMatrix.set(jobclass, jobclass, T1, P2, 1.0);
routingMatrix.set(jobclass, jobclass, T1, P3, 1.0);
routingMatrix.set(jobclass, jobclass, T2, P1, 1.0);
routingMatrix.set(jobclass, jobclass, T3, P1, 1.0);

model.link(routingMatrix);

%% Set Initial State    
P1.setState(jobclass.population);
P2.setState(0);
P3.setState(0);

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
