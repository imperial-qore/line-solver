%% Example: LQN2QN - Convert Layered Network to Queueing Network
%
% This example demonstrates how to convert a Layered Queueing Network (LQN)
% model to an equivalent Queueing Network (QN) using the LQN2QN converter.
%
% The LQN2QN conversion approximates synchronous call blocking by aggregating
% service times: when a task makes a synchCall, its effective service time
% includes its own processing plus all downstream processing (since it blocks
% until the reply).
%
% This approximation correctly models throughput and provides reasonable
% estimates of response times by aggregating them into entry tasks.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear; clc;

%% Create a simple 2-tier LQN model
% This models a basic client-server architecture where:
% - RefTask: reference (think) task with N jobs
% - ClientTask: processes requests, makes synchCall to server
% - ServerTask: handles requests from client
%
% Architecture:
%   [RefTask] <--> [ClientTask] <--> [ServerTask]

fprintf('Creating simple 2-tier LQN model...\n');

% Create the LQN
lqn = LayeredNetwork('TwoTierLQN');

% Define processors
thinkProc = Processor(lqn, 'ThinkProc', Inf, SchedStrategy.INF);
clientProc = Processor(lqn, 'ClientProc', 1, SchedStrategy.FCFS);
serverProc = Processor(lqn, 'ServerProc', 1, SchedStrategy.FCFS);

% Define tasks
refTask = Task(lqn, 'RefTask', 2, SchedStrategy.REF).on(thinkProc);
refTask.setThinkTime(Exp(0.5));

clientTask = Task(lqn, 'ClientTask', 1, SchedStrategy.FCFS).on(clientProc);
serverTask = Task(lqn, 'ServerTask', 1, SchedStrategy.FCFS).on(serverProc);

% Define entries
refEntry = Entry(lqn, 'RefEntry').on(refTask);
clientEntry = Entry(lqn, 'ClientEntry').on(clientTask);
serverEntry = Entry(lqn, 'ServerEntry').on(serverTask);

% Define activities
refActivity = Activity(lqn, 'RefActivity', Immediate()).on(refTask);
refActivity.boundTo(refEntry).synchCall(clientEntry, 1.0);

clientActivity = Activity(lqn, 'ClientActivity', Exp(0.2)).on(clientTask);
clientActivity.boundTo(clientEntry).synchCall(serverEntry, 1.0).repliesTo(clientEntry);

serverActivity = Activity(lqn, 'ServerActivity', Exp(0.3)).on(serverTask);
serverActivity.boundTo(serverEntry).repliesTo(serverEntry);

fprintf('LQN model created with 3 tasks and 3 entries.\n\n');

%% Display LQN structure
fprintf('=== LQN Structure ===\n');
fprintf('Processors:\n');
fprintf('  - ThinkProc (infinite capacity)\n');
fprintf('  - ClientProc (1 server)\n');
fprintf('  - ServerProc (1 server)\n\n');

fprintf('Tasks:\n');
fprintf('  - RefTask: N=2 jobs on ThinkProc\n');
fprintf('  - ClientTask: 1 task on ClientProc\n');
fprintf('  - ServerTask: 1 task on ServerProc\n\n');

fprintf('Service times:\n');
fprintf('  - RefTask think time: 0.5\n');
fprintf('  - ClientTask processing: 0.2 (plus blocked on server)\n');
fprintf('  - ServerTask processing: 0.3\n\n');

%% Convert LQN to QN
fprintf('=== Converting LQN to QN ===\n');
try
    qn = LQN2QN(lqn);
    fprintf('Conversion successful!\n\n');

    % Display converted model
    fprintf('Converted QN has:\n');
    fprintf('  - %d nodes\n', length(qn.nodes));
    fprintf('  - %d job classes\n', length(qn.classes));

    % List nodes
    fprintf('\nQN Nodes:\n');
    for n = 1:length(qn.nodes)
        nodeObj = qn.nodes{n};
        fprintf('  - %s (%s)\n', nodeObj.name, class(nodeObj));
    end

    fprintf('\nJob Classes:\n');
    for c = 1:length(qn.classes)
        classObj = qn.classes{c};
        if isa(classObj, 'ClosedClass')
            fprintf('  - %s (N=%d)\n', classObj.name, classObj.population);
        else
            fprintf('  - %s\n', classObj.name);
        end
    end
catch ME
    fprintf('Error during conversion: %s\n', ME.message);
    return;
end

%% Solve converted QN with DES
fprintf('\n=== Solving Converted QN with DES ===\n');
try
    options = Solver.defaultOptions;
    options.seed = 23000;
    options.samples = 10000;
    options.verbose = 0;

    solver = SolverDES(qn, options);
    avgTable = solver.getAvgTable();

    fprintf('DES solution obtained.\n');
    fprintf('Average metrics:\n');
    disp(avgTable);

catch ME
    fprintf('Error during DES solution: %s\n', ME.message);
end

%% Solve original LQN with LQNS (if available) for comparison
fprintf('\n=== Solving Original LQN with LQNS ===\n');
try
    solverLQNS = SolverLQNS(lqn, 'method', 'lqns');
    avgTableLQNS = solverLQNS.getAvgTable();

    fprintf('LQNS solution obtained.\n');
    fprintf('Average metrics:\n');
    disp(avgTableLQNS);

    % Compare key metrics
    fprintf('\n=== Comparison: LQN2QN+DES vs LQNS ===\n');
    fprintf('This comparison shows the accuracy of the blocking approximation.\n');
    fprintf('Note: DES results vary due to simulation randomness.\n');

catch ME
    fprintf('LQNS not available: %s\n', ME.message);
end

fprintf('\n=== Example Complete ===\n');
fprintf('The LQN2QN converter successfully converted the 2-tier LQN to a QN.\n');
fprintf('The converted model approximates synchronous call blocking by\n');
fprintf('aggregating service times, allowing analysis with QN solvers.\n');
