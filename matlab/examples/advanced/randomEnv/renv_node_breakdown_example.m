% Example demonstrating the node breakdown/repair API for random environments
%
% This example shows how to easily model a server that can break down and be
% repaired using the new Env convenience methods.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% Create base queueing network model
model = Network('ServerWithFailures');

% Define nodes
source = Source(model, 'Arrivals');
queue = Queue(model, 'Server', SchedStrategy.FCFS);
sink = Sink(model, 'Departures');

% Define job class
jobclass = OpenClass(model, 'Jobs');

% Set service and arrival rates (UP state)
source.setArrival(jobclass, Exp(0.8));  % Arrival rate
queue.setService(jobclass, Exp(2.0));   % Service rate when UP
queue.setNumberOfServers(1);

% Set routing
P = model.initRoutingMatrix();
P{jobclass,jobclass} = [0, 1, 0;   % Source -> Queue
                        0, 0, 1;   % Queue -> Sink
                        0, 0, 0];  % Sink -> (absorbing)
model.link(P);

%% Create environment with breakdown/repair using new API

% Method 1: Using addNodeFailureRepair (recommended) - passing node object
fprintf('Example 1: Using addNodeFailureRepair with node object\n');
env1 = Environment('ServerEnv1');

% Add failure and repair for the server node in one call
% Parameters: baseModel, nodeOrName, breakdownDist, repairDist, downServiceDist
% Note: Can pass either node object (queue) or node name ('Server')
env1.addNodeFailureRepair(model, queue, Exp(0.1), Exp(1.0), Exp(0.5));

% Initialize and print stage table
env1.init();
fprintf('\nStage table for env1:\n');
env1.getStageTable()

%% Method 2: Using separate breakdown and repair calls
fprintf('\nExample 2: Using separate breakdown and repair calls with node object\n');
env2 = Environment('ServerEnv2');

% Add breakdown (creates UP and DOWN stages) - using node object
env2.addNodeBreakdown(model, queue, Exp(0.1), Exp(0.5));

% Add repair transition - using node object
env2.addNodeRepair(queue, Exp(1.0));

% Initialize and print stage table
env2.init();
fprintf('\nStage table for env2:\n');
env2.getStageTable()

%% Method 3: With custom reset policies
fprintf('\nExample 3: With custom reset policies using node object\n');
env3 = Environment('ServerEnv3');

% Reset policy: clear all queues on breakdown
resetBreakdown = @(q) 0*q;  % Clear queues when server breaks down

% Reset policy: keep all jobs on repair
resetRepair = @(q) q;  % Keep jobs when server is repaired

% Using node object instead of string name
env3.addNodeFailureRepair(model, queue, Exp(0.1), Exp(1.0), Exp(0.5), ...
                         resetBreakdown, resetRepair);

env3.init();
fprintf('\nStage table for env3:\n');
env3.getStageTable()

%% Method 4: Modifying reset policies after creation
fprintf('\nExample 4: Modifying reset policies after creation using node object\n');
env4 = Environment('ServerEnv4');

% Create environment with default reset policies - using node object
env4.addNodeFailureRepair(model, queue, Exp(0.1), Exp(1.0), Exp(0.5));

% Update breakdown reset policy to clear queues - using node object
env4.setBreakdownResetPolicy(queue, @(q) 0*q);

% Update repair reset policy (keep jobs) - using node object
env4.setRepairResetPolicy(queue, @(q) q);

env4.init();
fprintf('\nStage table for env4:\n');
env4.getStageTable()

%% Solve the environment model
fprintf('\nSolving environment model with ENV solver...\n');

% Create solver options
options = Solver.defaultOptions;
options.timespan = [0, Inf];
options.iter_max = 100;
options.iter_tol = 0.01;
options.method = 'default';
options.verbose = true;

% Create fluid solver for each stage
sfoptions = FLD.defaultOptions;
sfoptions.timespan = [0, 1e3];
sfoptions.verbose = false;

% Solve using ENV (Ensemble eNVironment) solver
envSolver = ENV(env1, @(model) FLD(model, sfoptions), options);

% Get and display results
fprintf('\nAverage Performance Metrics:\n');
[QN, UN, TN] = envSolver.getAvg();
AvgTable = envSolver.getAvgTable()

fprintf('\nInterpretation:\n');
fprintf('- The system alternates between UP (operational) and DOWN (failed) states\n');
fprintf('- UP state: Server processes jobs at rate 2.0\n');
fprintf('- DOWN state: Server processes jobs at reduced rate 0.5\n');
fprintf('- Breakdown occurs at rate 0.1 (mean time to failure = 10 time units)\n');
fprintf('- Repair occurs at rate 1.0 (mean time to repair = 1 time unit)\n');
fprintf('- Results show averaged performance across both states\n');

%% Compute Reliability Metrics
fprintf('\nReliability Metrics:\n');
ReliabilityTable = env1.getReliabilityTable()

fprintf('\nReliability Interpretation:\n');
fprintf('- MTTF (Mean Time To Failure): %.2f time units\n', ReliabilityTable.Value(1));
fprintf('  Expected time in UP state before breakdown occurs\n');
fprintf('- MTTR (Mean Time To Repair): %.2f time units\n', ReliabilityTable.Value(2));
fprintf('  Expected time in DOWN state before repair completes\n');
fprintf('- MTBF (Mean Time Between Failures): %.2f time units\n', ReliabilityTable.Value(3));
fprintf('  Average cycle time (MTTF + MTTR)\n');
fprintf('- Availability: %.4f (%.2f%%)\n', ReliabilityTable.Value(4), ReliabilityTable.Value(4)*100);
fprintf('  Fraction of time the system is operational\n');
