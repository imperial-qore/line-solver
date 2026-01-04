clear P T R solver AvgTable

model = Network('model');

%% Nodes
source = Source(model,'Source');
sink = Sink(model,'Sink');

P{1} = Place(model, 'P1');
P{2} = Place(model, 'P2');
P{3} = Place(model, 'P3');
P{4} = Place(model, 'P4');
P{5} = Place(model, 'P5');
P{6} = Place(model, 'P6');
P{7} = Place(model, 'P7');

T{1} = Transition(model, 'T1');
T{2} = Transition(model, 'T2');
T{3} = Transition(model, 'T3');
T{4} = Transition(model, 'T4');
T{5} = Transition(model, 'T5');
T{6} = Transition(model, 'T6');
T{7} = Transition(model, 'T7');
T{8} = Transition(model, 'T8');

% Source
jobclass{1} = OpenClass(model, 'Class1', 0);
source.setArrival(jobclass{1}, Exp.fitMean(1));


%% Parameterization 

% T1
mode = T{1}.addMode('Mode1');
T{1}.setNumberOfServers(mode,Inf);
T{1}.setDistribution(1,Exp(4));
T{1}.setEnablingConditions(1,jobclass{1},P{1},1);
T{1}.setFiringOutcome(1,jobclass{1},P{2},1);

% T2
mode = T{2}.addMode('Mode1');
T{2}.setNumberOfServers(mode,Inf);
T{2}.setEnablingConditions(1,jobclass{1},P{2},1);
T{2}.setFiringOutcome(1,jobclass{1},P{3},1);
T{2}.setTimingStrategy(1,TimingStrategy.IMMEDIATE);
T{2}.setFiringPriorities(1,1);
T{2}.setFiringWeights(1,1);

% T3
mode = T{3}.addMode('Mode1');
T{3}.setNumberOfServers(mode,Inf);
T{3}.setEnablingConditions(1,jobclass{1},P{2},1);
T{3}.setFiringOutcome(1,jobclass{1},P{4},1);
T{3}.setTimingStrategy(1,TimingStrategy.IMMEDIATE);
T{3}.setFiringPriorities(1,1);
% T{3}.setFiringWeights(1,0.6);
T{3}.setFiringPriorities(1,1);

% T4
mode = T{4}.addMode('Mode1');
T{4}.setNumberOfServers(mode,Inf);
T{4}.setEnablingConditions(1,jobclass{1},P{3},1);
T{4}.setEnablingConditions(1,jobclass{1},P{5},1);
T{4}.setFiringOutcome(1,jobclass{1},P{5},1);
T{4}.setFiringOutcome(1,jobclass{1},P{6},1);
T{4}.setTimingStrategy(1,TimingStrategy.IMMEDIATE);
T{4}.setFiringPriorities(1,1);

% T5
mode = T{5}.addMode('Mode1');
T{5}.setNumberOfServers(mode,Inf);
T{5}.setEnablingConditions(1,jobclass{1},P{4},1);
T{5}.setEnablingConditions(1,jobclass{1},P{5},1);
T{5}.setFiringOutcome(1,jobclass{1},P{7},1);
T{5}.setInhibitingConditions(1,jobclass{1},P{6},1);
T{5}.setTimingStrategy(1,TimingStrategy.IMMEDIATE);
T{5}.setFiringPriorities(1,1);

% T6
mode = T{6}.addMode('Mode1');
T{6}.setNumberOfServers(mode,Inf);
T{6}.setDistribution(1,Erlang(2,2));
T{6}.setEnablingConditions(1,jobclass{1},P{6},1);
T{6}.setFiringOutcome(1,jobclass{1},P{1},1);

% T7
mode = T{7}.addMode('Mode1');
T{7}.setNumberOfServers(mode,Inf);
T{7}.setDistribution(1,Exp(2));
T{7}.setEnablingConditions(1,jobclass{1},P{7},1);
T{7}.setFiringOutcome(1,jobclass{1},P{1},1);
T{7}.setFiringOutcome(1,jobclass{1},P{5},1);

% T8
mode = T{8}.addMode('Mode1');
T{8}.setNumberOfServers(mode,Inf);
T{8}.setDistribution(1,Exp(2));
T{8}.setEnablingConditions(1,jobclass{1},P{4},1);
T{8}.setFiringOutcome(1,jobclass{1},sink,1);

%% Routing 
M = model.getNumberOfStations();
K = model.getNumberOfClasses();

R = model.initRoutingMatrix(); % initialize routing matrix 

R{1,1}(source,P{1}) = 1; % (Source,Class1) -> (P1,Class1)

R{1,1}(P{1},T{1}) = 1; % (P1,Class1) -> (T1,Class1)
R{1,1}(P{2},T{2}) = 1; % (P2,Class1) -> (T2,Class1)
R{1,1}(P{2},T{3}) = 1; % (P2,Class1) -> (T3,Class1)
R{1,1}(P{3},T{4}) = 1; % (P3,Class1) -> (T4,Class1)
R{1,1}(P{4},T{5}) = 1; % (P4,Class1) -> (T5,Class1)
R{1,1}(P{5},T{4}) = 1; % (P5,Class1) -> (T4,Class1)
R{1,1}(P{5},T{5}) = 1; % (P5,Class1) -> (T5,Class1)
R{1,1}(P{6},T{5}) = 1; % (P6,Class1) -> (T5,Class1)
R{1,1}(P{6},T{6}) = 1; % (P6,Class1) -> (T6,Class1)
R{1,1}(P{7},T{7}) = 1; % (P7,Class1) -> (T7,Class1)

R{1,1}(T{1},P{2}) = 1; % (T1,Class1) -> (P2,Class1)
R{1,1}(T{2},P{3}) = 1; % (T2,Class1) -> (P3,Class1)
R{1,1}(T{3},P{4}) = 1; % (T3,Class1) -> (P4,Class1)
R{1,1}(T{4},P{5}) = 1; % (T4,Class1) -> (P5,Class1)
R{1,1}(T{4},P{6}) = 1; % (T4,Class1) -> (P6,Class1)
R{1,1}(T{5},P{7}) = 1; % (T5,Class1) -> (P7,Class1)
R{1,1}(T{6},P{1}) = 1; % (T6,Class1) -> (P1,Class1)
R{1,1}(T{7},sink) = 1; % (T7,Class1) -> (Sink,Class1)
R{1,1}(T{7},P{1}) = 1; % (T7,Class1) -> (P1,Class1)
R{1,1}(T{7},P{5}) = 1; % (T7,Class1) -> (P5,Class1)

R{1,1}(P{4},T{8}) = 1; % (P4,Class1) -> (T5,Class1)
R{1,1}(T{8},sink) = 1; % (T8,Class1) -> (Sink,Class1)

model.link(R);


%% Set Initial State    

source.setState(0);
P{1}.setState(2);
P{2}.setState(0);
P{3}.setState(0);
P{4}.setState(0);
P{5}.setState(1);
P{6}.setState(0);
P{7}.setState(0);

state = model.getState;

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
