% Example: M/M/1 Queue with Exponential Impatience (Reneging)
%
% This example demonstrates customer impatience in a simple M/M/1 queue.
% Jobs arrive according to a Poisson process and have exponentially
% distributed service times. If a job waits too long in the queue, it
% reneges (abandons) according to an exponential patience distribution.
%
% Model:
% - Single queue with FCFS scheduling
% - Exponential service time (rate = 1.0)
% - Exponential arrival process (rate = 0.8)
% - Exponential patience/reneging time (rate = 0.5)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear node jobclass solver AvgTable

model = Network('M/M/1 with Exponential Impatience');

% Define nodes
node{1} = Source(model, 'Source');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{3} = Sink(model, 'Sink');

% Define job class
jobclass{1} = OpenClass(model, 'Class1', 0);

% Set service times
node{2}.setService(jobclass{1}, Exp(1.0));  % Mean service time = 1.0

% Set arrival rate
node{1}.setArrival(jobclass{1}, Exp(0.8));  % Arrival rate = 0.8

% Set patience distribution (queue-specific)
% Jobs will renege after an exponentially distributed time with rate 0.5
node{2}.setPatience(jobclass{1}, Exp(0.5));  % Mean patience = 2.0

% Define routing
K = model.getNumberOfClasses();
P = cell(K, K);
P{1,1} = [0, 1, 0;   % Source -> Queue1
          0, 0, 1;   % Queue1 -> Sink
          0, 0, 0];  % Sink -> nowhere

model.link(P);

%% Solve with JMT
options = Solver.defaultOptions;
options.keep = true;
options.verbose = 1;
options.samples = 1e5;
options.seed = 23000;

disp('This example shows a simple M/M/1 queue with exponential impatience.');
disp('Jobs will abandon (renege) the queue after waiting for an exponential time.');
disp(' ');
disp('Model parameters:');
disp('  - Arrival rate: 0.8');
disp('  - Service rate: 1.0');
disp('  - Patience rate: 0.5 (mean patience = 2.0)');
disp(' ');

solver{1} = JMT(model, options);

fprintf(1, 'SOLVER: %s\n', solver{1}.getName());
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}

disp(' ');
disp('Note: The queue length and response time should be lower than a standard M/M/1');
disp('      due to jobs reneging when their patience expires.');
