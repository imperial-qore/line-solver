clear P T R solver AvgTable jobclass

model = Network('model');

%% Nodes
source = Source(model,'Source');
sink = Sink(model,'Sink');

P{1} = Place(model, 'P1');

T{1} = Transition(model, 'T1');

% Source
jobclass{1} = OpenClass(model, 'Class1', 0);
source.setArrival(jobclass{1}, Exp(0.5)); % arrival rate 0.5

%% Parameterisation

% T1 with Pareto service time
% Pareto(shape, scale) - shape must be >= 2
% With shape=3 and scale=1, mean service time = 3*1/(3-1) = 1.5
mode = T{1}.addMode('Mode1');
T{1}.setNumberOfServers(mode, 1);
T{1}.setDistribution(mode, Pareto(3, 1)); % Pareto with shape=3, scale=1
T{1}.setEnablingConditions(mode, jobclass{1}, P{1}, 1);
T{1}.setFiringOutcome(mode, jobclass{1}, sink, 1);

% Routing
R = model.initRoutingMatrix(); % initialize routing matrix

R{1,1}(source, P{1}) = 1; % (Source,Class1) -> (P1,Class1)
R{1,1}(P{1}, T{1}) = 1;
R{1,1}(T{1}, sink) = 1;

model.link(R);

%% Solver
options = Solver.defaultOptions;
options.verbose = 1;
options.seed = 23000;
options.samples = 1e4; % limit samples for faster execution

solver = {};
solver{1} = JMT(model, options);
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}
