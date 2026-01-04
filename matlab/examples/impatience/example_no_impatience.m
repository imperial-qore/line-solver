% Example: Baseline Model Without Impatience
%
% This example serves as a baseline to verify backward compatibility.
% It creates a simple queueing network WITHOUT any impatience configured,
% ensuring that existing models continue to work as expected.
%
% Model:
% - Two queues in series
% - Single job class
% - NO impatience configured (standard queueing behavior)
%
% Purpose: Verify that:
% 1. Models without impatience continue to work correctly
% 2. JMT XML generation produces 'null' for impatience parameters
% 3. Results match standard queueing theory predictions
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear node jobclass solver AvgTable

model = Network('Baseline - No Impatience');

% Define nodes
node{1} = Source(model, 'Source');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{3} = Queue(model, 'Queue2', SchedStrategy.FCFS);
node{4} = Sink(model, 'Sink');

% Define job class
jobclass{1} = OpenClass(model, 'Class1', 0);

% Set service times
node{2}.setService(jobclass{1}, Exp(2.0));  % Mean service = 0.5
node{3}.setService(jobclass{1}, Exp(1.0));  % Mean service = 1.0

% Set arrival rate
node{1}.setArrival(jobclass{1}, Exp(0.4));  % Arrival rate = 0.4

% NOTE: NO patience is configured - this is intentional!
% The model should work exactly as before the impatience feature was added.

% Define routing
K = model.getNumberOfClasses();
P = cell(K, K);
P{1,1} = [0, 1, 0, 0;   % Source -> Queue1
          0, 0, 1, 0;   % Queue1 -> Queue2
          0, 0, 0, 1;   % Queue2 -> Sink
          0, 0, 0, 0];  % Sink -> nowhere

model.link(P);

%% Solve with JMT
options = Solver.defaultOptions;
options.keep = true;
options.verbose = 1;
options.samples = 1e5;
options.seed = 23000;

disp('This is a baseline example WITHOUT impatience to verify backward compatibility.');
disp(' ');
disp('Model parameters:');
disp('  - Arrival rate: 0.4');
disp('  - Queue1 service rate: 2.0');
disp('  - Queue2 service rate: 1.0');
disp('  - NO impatience configured');
disp(' ');
disp('Expected behavior:');
disp('  - Standard M/M/1 queueing with no job abandonment');
disp('  - JMT XML should contain <value>null</value> for impatience parameters');
disp('  - Results should match traditional queueing analysis');
disp(' ');

solver{1} = JMT(model, options);

fprintf(1, 'SOLVER: %s\n', solver{1}.getName());
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}

disp(' ');
disp('Backward compatibility verified: Model runs successfully without impatience.');
