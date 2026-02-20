clear P T R solver AvgTable jobclass

model = Network('model');

%% Nodes
P{1} = Place(model, 'P1');
T{1} = Transition(model, 'T1');

% Source
jobclass{1} = ClosedClass(model, 'Class1', 1, P{1});

%% Parameterisation

% T1
mode1 = T{1}.addMode('Mode1');
T{1}.setDistribution(mode1,Exp.fitMean(1));
T{1}.setEnablingConditions(mode1,jobclass{1},P{1},1);
T{1}.setFiringOutcome(mode1,jobclass{1},P{1},1);

mode2 = T{1}.addMode('Mode2');
T{1}.setDistribution(mode2,Erlang.fitMeanAndOrder(1,2));
T{1}.setEnablingConditions(mode2,jobclass{1},P{1},1);
T{1}.setFiringOutcome(mode2,jobclass{1},P{1},1);

mode3 = T{1}.addMode('Mode3');
T{1}.setDistribution(mode3,HyperExp.fitMeanAndSCV(1,4));
T{1}.setEnablingConditions(mode3,jobclass{1},P{1},1);
T{1}.setFiringOutcome(mode3,jobclass{1},P{1},1);

% Routing
M = model.getNumberOfStations();
K = model.getNumberOfClasses();

R = model.initRoutingMatrix(); % initialize routing matrix

R{1,1}(P{1},T{1}) = 1;
R{1,1}(T{1},P{1}) = 1;

model.link(R);

%% Set Initial State
P{1}.setState(jobclass{1}.population);

state = model.getState;

%% Solver
options = Solver.defaultOptions;
options.keep=2;
options.verbose=1;
%options.cutoff = 10;
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
if ~isempty(AvgTable{1})
    AvgTable{1}
end
