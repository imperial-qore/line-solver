% Example 11: Random environments and SolverENV
% This tutorial illustrates how to model a queueing system operating in
% a random environment, where system parameters (e.g., service rates)
% change according to an underlying environmental process.
%
% Scenario: A server that alternates between "Fast" and "Slow" modes.
% In Fast mode, service rate is 4.0. In Slow mode, service rate is 1.0.
% The environment switches from Fast->Slow at rate 0.5 and Slow->Fast at rate 1.0.

warning off;

%% Block 1: Create base network model
% First, we define a closed queueing network with a delay and a queue
baseModel = Network('BaseModel');
delay = Delay(baseModel, 'ThinkTime');
queue = Queue(baseModel, 'Fast/Slow Server', SchedStrategy.FCFS);

% Closed class with 5 jobs
N = 5;
jobclass = ClosedClass(baseModel, 'Jobs', N, delay);
delay.setService(jobclass, Exp(1.0));  % Think time = 1.0
queue.setService(jobclass, Exp(2.0));  % Placeholder service rate

% Connect nodes in a cycle
baseModel.link(Network.serialRouting(delay, queue));

%% Block 2: Create the random environment
% Define two stages: "Fast" and "Slow" with different service rates
env = Environment('ServerModes');

% Stage 1: Fast mode (service rate = 4.0)
fastModel = baseModel.copy();
fastQueue = fastModel.getNodeByName('Fast/Slow Server');
fastQueue.setService(fastModel.classes{1}, Exp(4.0));
env.addStage('Fast', 'operational', fastModel);

% Stage 2: Slow mode (service rate = 1.0)
slowModel = baseModel.copy();
slowQueue = slowModel.getNodeByName('Fast/Slow Server');
slowQueue.setService(slowModel.classes{1}, Exp(1.0));
env.addStage('Slow', 'degraded', slowModel);

% Define transitions between stages
% Fast -> Slow at rate 0.5 (mean time in Fast mode = 2.0)
env.addTransition('Fast', 'Slow', Exp(0.5));
% Slow -> Fast at rate 1.0 (mean time in Slow mode = 1.0)
env.addTransition('Slow', 'Fast', Exp(1.0));

%% Block 3: Inspect the environment structure
% Display the environment stages and their probabilities
env.init();
stageTable = env.getStageTable()

%% Block 4: Solve using SolverENV
% SolverENV requires a solver factory that creates solvers for each stage
% We use the Fluid solver (FLD) with transient analysis

options = Solver.defaultOptions;
options.iter_max = 50;
options.iter_tol = 0.01;
options.verbose = false;

% Solver factory: creates a Fluid solver for each stage model
fldOptions = FLD.defaultOptions;
fldOptions.timespan = [0, 100];
fldOptions.verbose = false;
solverFactory = @(m) FLD(m, fldOptions);

% Create and run the ENV solver
envSolver = ENV(env, solverFactory, options);
[QN, UN, TN] = envSolver.getAvg();

% Display average results weighted by environment probabilities
fprintf('\n--- Environment-Averaged Results ---\n');
envAvgTable = envSolver.getAvgTable()

%% Block 5: Compare with individual stage analysis
% Analyze each stage network in steady-state using MVA
fprintf('\n--- Individual Stage Analysis (MVA) ---\n');
for e = 1:height(env.envGraph.Nodes)
    stageName = env.envGraph.Nodes.Name{e};
    stageModel = env.envGraph.Nodes.Model{e};
    fprintf('\nStage: %s (prob = %.4f)\n', stageName, env.probEnv(e));
    mvaSolver = MVA(stageModel);
    stageTable = mvaSolver.getAvgTable();
    disp(stageTable);
end
