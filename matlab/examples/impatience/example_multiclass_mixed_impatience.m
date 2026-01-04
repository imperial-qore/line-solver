% Example: Multi-Class Queue with Mixed Impatience Distributions
%
% This example demonstrates customer impatience with different patience
% distributions for different job classes in a multi-class queueing system.
%
% Model:
% - Two queues in tandem
% - Three job classes with different characteristics:
%   * Class 1: Exponential patience (impatient customers)
%   * Class 2: Deterministic patience (timeout-based abandonment)
%   * Class 3: Erlang-2 patience (more patient, less variable)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear node jobclass solver AvgTable

model = Network('Multi-Class with Mixed Impatience');

% Define nodes
node{1} = Source(model, 'Source');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{3} = Queue(model, 'Queue2', SchedStrategy.FCFS);
node{4} = Sink(model, 'Sink');

% Define job classes
jobclass{1} = OpenClass(model, 'ImpatientClass', 0);
jobclass{2} = OpenClass(model, 'TimeoutClass', 0);
jobclass{3} = OpenClass(model, 'PatientClass', 0);

% Set service times (same for all classes at both queues)
for c = 1:3
    node{2}.setService(jobclass{c}, Exp(2.0));  % Mean service = 0.5
    node{3}.setService(jobclass{c}, Exp(1.5));  % Mean service = 0.67
end

% Set arrival rates
node{1}.setArrival(jobclass{1}, Exp(0.3));  % Impatient arrivals
node{1}.setArrival(jobclass{2}, Exp(0.2));  % Timeout arrivals
node{1}.setArrival(jobclass{3}, Exp(0.1));  % Patient arrivals

% Set patience distributions (different for each class at Queue1)
% Class 1: Exponential patience - highly variable, mean = 2.0
node{2}.setPatience(jobclass{1}, Exp(0.5));

% Class 2: Deterministic patience - fixed timeout of 3.0 time units
node{2}.setPatience(jobclass{2}, Det(3.0));

% Class 3: Erlang-2 patience - less variable than exponential, mean = 4.0
node{2}.setPatience(jobclass{3}, Erlang(2, 0.5));

% Queue2 also has impatience, but with different parameters
node{3}.setPatience(jobclass{1}, Exp(0.3));     % More impatient at second queue
node{3}.setPatience(jobclass{2}, Det(5.0));     % Longer timeout
node{3}.setPatience(jobclass{3}, Erlang(2, 0.3));  % Even more patient

% Define routing
K = model.getNumberOfClasses();
P = cell(K, K);
for c = 1:K
    P{c,c} = [0, 1, 0, 0;   % Source -> Queue1
              0, 0, 1, 0;   % Queue1 -> Queue2
              0, 0, 0, 1;   % Queue2 -> Sink
              0, 0, 0, 0];  % Sink -> nowhere
end

model.link(P);

%% Solve with JMT
options = Solver.defaultOptions;
options.keep = true;
options.verbose = 1;
options.samples = 1e5;
options.seed = 23000;

disp('This example shows a multi-class system with different patience distributions.');
disp(' ');
disp('Class 1 (Impatient): Exponential patience, mean = 2.0 at Queue1, 3.33 at Queue2');
disp('Class 2 (Timeout):   Deterministic patience, timeout = 3.0 at Queue1, 5.0 at Queue2');
disp('Class 3 (Patient):   Erlang-2 patience, mean = 4.0 at Queue1, 6.67 at Queue2');
disp(' ');

solver{1} = JMT(model, options);

fprintf(1, 'SOLVER: %s\n', solver{1}.getName());
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}

disp(' ');
disp('Note: Different classes will have different reneging rates due to their');
disp('      patience distributions. Class 1 should see the most reneging.');
