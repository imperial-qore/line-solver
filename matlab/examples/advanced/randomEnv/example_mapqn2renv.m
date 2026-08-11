% Example: Transform MMPP service queue to random environment model
%
% This example shows how to transform a closed queueing network with MMPP2
% service into a random environment model with exponential services
% modulated by the MMPP phases.
%
% MMPP2 (2-phase Markov Modulated Poisson Process) is a point process
% where service completions occur at rates modulated by a 2-state Markov chain.
% This transformation converts such a network into an equivalent random
% environment model where:
%   - Environment has 2 stages (one per MMPP phase)
%   - Service rates are exponential with rates from MMPP D1 diagonal
%   - Environment transitions follow MMPP D0 matrix structure

clear; close all;

% Create closed queueing network with MMPP2 service
model = Network('MMPP_ClosedQN');

% Create nodes
delay = Delay(model, 'Delay');      % Think time station
queue = Queue(model, 'Queue', SchedStrategy.FCFS);  % Service station with MMPP2

% Create closed job class with N=5 jobs
N = 5;  % Number of jobs in the system
jobclass = ClosedClass(model, 'Class1', N, delay);

% Set think time (exponential with mean 1.0)
delay.setService(jobclass, Exp(1.0));

% Set MMPP2 service distribution
% Parameters: lambda0, lambda1, sigma01, sigma10
% lambda0, lambda1 = service rates in phases 0, 1
% sigma01, sigma10 = phase transition rates
%
% Using very different service rates to show distinct phase behavior:
%   Phase 0: slow service (rate 1.0) -> high queue length
%   Phase 1: fast service (rate 10.0) -> low queue length
%
% This creates an MMPP2 with:
%   D0 = [-0.2-1.0,  0.2 ]   = [-1.2,  0.2]
%        [ 0.3,  -0.3-10.0]    [ 0.3, -10.3]
%
%   D1 = [1.0,   0 ]
%        [ 0,  10.0]
queue.setService(jobclass, MMPP2(1.0, 10.0, 0.2, 0.3));

% Create cyclic routing: Delay -> Queue -> Delay
model.link(Network.serialRouting(delay, queue));

fprintf('=== MMPP2 Closed QN to Random Environment Transformation ===\n\n');

fprintf('Original Network:\n');
fprintf('  Topology: Delay -> MMPP2 Queue -> Delay (cyclic)\n');
fprintf('  Population: N = %d jobs\n', N);
fprintf('  Think time: Exp(1.0)\n');
fprintf('  Service: MMPP2 with phases 0,1\n\n');

% Transform to random environment model
fprintf('Transforming to random environment model...\n');
envModel = MAPQN2RENV(model);

fprintf('Transformation complete!\n\n');

fprintf('Environment Model Structure:\n');
fprintf('  Environment: MAPQN_Env with 2 stages\n');
fprintf('  Stage 0 (Phase 0): Closed QN with Exp service rate 1.0 (SLOW)\n');
fprintf('  Stage 1 (Phase 1): Closed QN with Exp service rate 10.0 (FAST)\n');
fprintf('  Transitions:\n');
fprintf('    Phase 0 -> Phase 1: Rate 0.2\n');
fprintf('    Phase 1 -> Phase 0: Rate 0.3\n\n');

% Display environment model structure
fprintf('Environment stages:\n');
stageNames = envModel.envGraph.Nodes.Name;
stageTypes = envModel.envGraph.Nodes.Type;
for i = 1:length(stageNames)
    fprintf('  %s: %s\n', stageNames{i}, stageTypes{i});
end
fprintf('\n');

% Display the stage network structure
fprintf('Stage network details:\n');
stageModels = envModel.envGraph.Nodes.Model;
for i = 1:length(stageModels)
    stageModel = stageModels{i};
    if ~isempty(stageModel)
        fprintf('  %s: Network with %d nodes\n', stageNames{i}, length(stageModel.nodes));
    end
end
fprintf('\n');

% Solve the original MMPP2 model using JMT
try
    fprintf('Solving original MMPP2 model using JMT...\n');
    solver = JMT(model);
    avgTable = solver.getAvgTable();
    fprintf('Original MMPP2 Queue results:\n');
    disp(avgTable);
catch ME
    fprintf('JMT solver error: %s\n', ME.message);
end

% Solve the environment model using ENV with FLD
fprintf('\nSolving environment model using ENV (with FLD)...\n');
try
    % Setup solver options (similar to renv_twostages_repairmen.m)
    options = Solver.defaultOptions;
    options.timespan = [0, Inf];
    options.iter_max = 100;
    options.iter_tol = 0.01;
    options.method = 'default';
    options.verbose = false;

    sfoptions = FLD.defaultOptions;
    sfoptions.timespan = [0, 100];
    sfoptions.verbose = false;

    envSolver = ENV(envModel, @(m) FLD(m, sfoptions), options);
    [QN, UN, TN] = envSolver.getAvg();
    AvgTable = envSolver.getAvgTable();

    fprintf('\nEnvironment Model AvgTable:\n');
    disp(AvgTable);

    fprintf('Note: The random environment model captures MMPP dynamics through\n');
    fprintf('switching between phases with exponential service rates.\n\n');
catch ME
    fprintf('ENV error: %s\n\n', ME.message);
    fprintf('Note: ENV requires transient analysis support.\n');
    fprintf('The environment model was created successfully.\n\n');
end

% Solve the environment model using ENV with CTMC
fprintf('\nSolving environment model using ENV (with CTMC, cutoff=100)...\n');
try
    % Setup solver options
    options = Solver.defaultOptions;
    options.timespan = [0, Inf];
    options.iter_max = 100;
    options.iter_tol = 0.01;
    options.method = 'default';
    options.verbose = false;

    ctmcoptions = CTMC.defaultOptions;
    ctmcoptions.timespan = [0, 100];
    ctmcoptions.cutoff = 100;
    ctmcoptions.verbose = false;

    envSolverCTMC = ENV(envModel, @(m) CTMC(m, ctmcoptions), options);
    [QN_ctmc, UN_ctmc, TN_ctmc] = envSolverCTMC.getAvg();
    AvgTableCTMC = envSolverCTMC.getAvgTable();

    fprintf('\nEnvironment Model AvgTable (CTMC):\n');
    disp(AvgTableCTMC);
catch ME
    fprintf('ENV/CTMC error: %s\n\n', ME.message);
end

% Solve individual stage networks using MVA (steady-state)
fprintf('Solving stage networks individually (steady-state with MVA)...\n');
try
    stageModels = envModel.envGraph.Nodes.Model;
    for i = 1:length(stageModels)
        stageModel = stageModels{i};
        stageName = envModel.envGraph.Nodes.Name{i};
        if ~isempty(stageModel)
            fprintf('\n  Stage %s:\n', stageName);
            stageSolver = MVA(stageModel);
            stageAvgTable = stageSolver.getAvgTable();
            disp(stageAvgTable);
        end
    end
catch ME2
    fprintf('Stage solver error: %s\n', ME2.message);
end

fprintf('\n=== Transformation Complete ===\n');
