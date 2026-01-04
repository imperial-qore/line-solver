% Example: Global Class-Level Patience
%
% This example demonstrates setting patience at the job class level rather
% than per-queue. When patience is set globally on a class, it applies to
% all queues that the class visits (unless overridden by queue-specific settings).
%
% Model:
% - Three queues in series
% - Two job classes:
%   * Class 1: Global exponential patience (applies to all queues)
%   * Class 2: Global deterministic patience + queue-specific override at Queue2
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear node jobclass solver AvgTable

model = Network('Global Patience Example');

% Define nodes
node{1} = Source(model, 'Source');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{3} = Queue(model, 'Queue2', SchedStrategy.FCFS);
node{4} = Queue(model, 'Queue3', SchedStrategy.FCFS);
node{5} = Sink(model, 'Sink');

% Define job classes
jobclass{1} = OpenClass(model, 'GlobalExpClass', 0);
jobclass{2} = OpenClass(model, 'OverrideClass', 0);

% Set GLOBAL patience at the class level
% This applies to ALL queues unless overridden
jobclass{1}.setPatience(Exp(0.4));   % Mean patience = 2.5 everywhere
jobclass{2}.setPatience(Det(4.0));   % Fixed patience = 4.0 everywhere (by default)

% Set service times
for c = 1:2
    node{2}.setService(jobclass{c}, Exp(1.0));
    node{3}.setService(jobclass{c}, Exp(1.2));
    node{4}.setService(jobclass{c}, Exp(0.8));
end

% Set arrival rates
node{1}.setArrival(jobclass{1}, Exp(0.3));
node{1}.setArrival(jobclass{2}, Exp(0.2));

% Override global patience for Class 2 at Queue2 only
% This demonstrates per-queue settings taking precedence over global settings
node{3}.setPatience(jobclass{2}, Exp(0.2));  % More impatient at Queue2

% Define routing
K = model.getNumberOfClasses();
P = cell(K, K);
for c = 1:K
    P{c,c} = [0, 1, 0, 0, 0;   % Source -> Queue1
              0, 0, 1, 0, 0;   % Queue1 -> Queue2
              0, 0, 0, 1, 0;   % Queue2 -> Queue3
              0, 0, 0, 0, 1;   % Queue3 -> Sink
              0, 0, 0, 0, 0];  % Sink -> nowhere
end

model.link(P);

%% Solve with JMT
options = Solver.defaultOptions;
options.keep = true;
options.verbose = 1;
options.samples = 1e5;
options.seed = 23000;

disp('This example demonstrates global class-level patience configuration.');
disp(' ');
disp('Class 1 (GlobalExpClass):');
disp('  - Has exponential patience (mean = 2.5) at ALL queues (global setting)');
disp(' ');
disp('Class 2 (OverrideClass):');
disp('  - Has deterministic patience (timeout = 4.0) at Queue1 and Queue3 (global setting)');
disp('  - Has exponential patience (mean = 5.0) at Queue2 (per-queue override)');
disp(' ');

solver{1} = JMT(model, options);

fprintf(1, 'SOLVER: %s\n', solver{1}.getName());
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}

disp(' ');
disp('Note: Global patience settings simplify configuration when the same');
disp('      patience behavior applies across multiple queues. Per-queue settings');
disp('      can override global settings when needed.');
