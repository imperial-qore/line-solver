clear P T R solver AvgTable jobclass

model = Network('model');

%% Nodes
source = Source(model,'Source');
sink = Sink(model,'Sink');

P{1} = Place(model, 'P1');

T{1} = Transition(model, 'T1');

% Source
jobclass{1} = OpenClass(model, 'Class1', 0);
source.setArrival(jobclass{1}, Exp(1));

%% Parameterisation 

% T1
mode = T{1}.addMode('Mode1');
T{1}.setNumberOfServers(mode,Inf);
T{1}.setDistribution(mode,Exp(4));
T{1}.setEnablingConditions(mode,jobclass{1},P{1},1);
T{1}.setFiringOutcome(mode,jobclass{1},sink,1);

% Routing 
M = model.getNumberOfStations();
K = model.getNumberOfClasses();

R = model.initRoutingMatrix(); % initialize routing matrix 

R{1,1}(source,P{1}) = 1; % (Source,Class1) -> (P1,Class1)

R{1,1}(P{1},T{1}) = 1;
R{1,1}(T{1},sink) = 1;

model.link(R);

%% Solver
options = Solver.defaultOptions;
options.keep=2;
options.verbose=1;
options.cutoff = 10;
options.seed = 23000;
%options.samples = 100;

% options.hide_immediate=1;
% options.is_pn=1;
% options.samples=2e4;

% All stations must be initialised.     
% initial_state = [0;2;0;0;0;1;0;0];

% solver = CTMC(model, options);
% solver.getAvgTable();

solver = {};
solver{1} = JMT(model,options);
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}
