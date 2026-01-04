clear node jobclass solver AvgTable;
%% Example of an open queueing network with 3 chains, 5 classes, and class switching
% Demonstrates job class switching within multiple open chains.
% Chain 1: Interactive requests (2 classes that switch)
% Chain 2: Batch jobs (2 classes that switch)
% Chain 3: Real-time tasks (1 class)

model = Network('oqn_multichain');

%% Block 1: nodes
node{1} = Source(model, 'RequestSource');
node{2} = Queue(model, 'Router', SchedStrategy.FCFS);
node{3} = Queue(model, 'WebCache', SchedStrategy.PS);
node{4} = Queue(model, 'AppServer', SchedStrategy.FCFS);
node{5} = Queue(model, 'DataServer', SchedStrategy.PS);
node{6} = Sink(model, 'ResponseSink');

%% Block 2: job classes
% Chain 1: Interactive requests (can switch between InteractiveA and InteractiveB)
jobclass{1} = OpenClass(model, 'InteractiveA', 0);
jobclass{2} = OpenClass(model, 'InteractiveB', 0);
% Chain 2: Batch jobs (can switch between BatchA and BatchB)
jobclass{3} = OpenClass(model, 'BatchA', 0);
jobclass{4} = OpenClass(model, 'BatchB', 0);
% Chain 3: Real-time tasks
jobclass{5} = OpenClass(model, 'RealTime', 0);

%% Block 3: service times
% Router
node{2}.setService(jobclass{1}, Exp.fitMean(0.15));
node{2}.setService(jobclass{2}, Exp.fitMean(0.18));
node{2}.setService(jobclass{3}, Exp.fitMean(0.2));
node{2}.setService(jobclass{4}, Exp.fitMean(0.22));
node{2}.setService(jobclass{5}, Exp.fitMean(0.1));

% WebCache
node{3}.setService(jobclass{1}, Exp.fitMean(0.3));
node{3}.setService(jobclass{2}, Exp.fitMean(0.35));
node{3}.setService(jobclass{3}, Exp.fitMean(1.0));
node{3}.setService(jobclass{4}, Exp.fitMean(1.2));
node{3}.setService(jobclass{5}, Exp.fitMean(0.2));

% AppServer
node{4}.setService(jobclass{1}, Exp.fitMean(0.5));
node{4}.setService(jobclass{2}, Exp.fitMean(0.6));
node{4}.setService(jobclass{3}, Exp.fitMean(1.5));
node{4}.setService(jobclass{4}, Exp.fitMean(1.8));
node{4}.setService(jobclass{5}, Exp.fitMean(0.3));

% DataServer
node{5}.setService(jobclass{1}, Exp.fitMean(0.4));
node{5}.setService(jobclass{2}, Exp.fitMean(0.45));
node{5}.setService(jobclass{3}, Exp.fitMean(0.8));
node{5}.setService(jobclass{4}, Exp.fitMean(1.0));
node{5}.setService(jobclass{5}, Exp.fitMean(0.25));

%% Block 4: arrival rates
% Interactive requests
node{1}.setArrival(jobclass{1}, Exp(0.8));
node{1}.setArrival(jobclass{2}, Exp(0.3));
% Batch jobs
node{1}.setArrival(jobclass{3}, Exp(0.5));
node{1}.setArrival(jobclass{4}, Exp(0.4));
% Real-time tasks
node{1}.setArrival(jobclass{5}, Exp(0.6));

%% Block 5: routing with class switching
M = model.getNumberOfNodes();
K = model.getNumberOfClasses();

P = model.initRoutingMatrix();

% ===== Chain 1: Interactive requests with switching =====
% InteractiveA routing:
P{1,1}(1,2) = 1;      % Source -> Router
P{1,1}(2,3) = 0.6;    % Router -> WebCache (stay in InteractiveA)
P{1,2}(2,3) = 0.4;    % Class switch to InteractiveB at Router
P{1,1}(3,4) = 1;      % WebCache -> AppServer
P{1,1}(4,5) = 1;      % AppServer -> DataServer
P{1,1}(5,6) = 1;      % DataServer -> Sink

% InteractiveB routing:
P{2,2}(1,2) = 1;      % Source -> Router
P{2,2}(2,3) = 0.7;    % Router -> WebCache (stay in InteractiveB)
P{2,1}(2,3) = 0.3;    % Class switch to InteractiveA at Router
P{2,2}(3,4) = 1;      % WebCache -> AppServer
P{2,2}(4,5) = 1;      % AppServer -> DataServer
P{2,2}(5,6) = 1;      % DataServer -> Sink

% ===== Chain 2: Batch jobs with switching =====
% BatchA routing:
P{3,3}(1,2) = 1;      % Source -> Router
P{3,3}(2,3) = 0.5;    % Router -> WebCache (stay in BatchA)
P{3,4}(2,3) = 0.5;    % Class switch to BatchB at Router
P{3,3}(3,4) = 1;      % WebCache -> AppServer
P{3,3}(4,5) = 1;      % AppServer -> DataServer
P{3,3}(5,6) = 1;      % DataServer -> Sink

% BatchB routing:
P{4,4}(1,2) = 1;      % Source -> Router
P{4,4}(2,3) = 0.6;    % Router -> WebCache (stay in BatchB)
P{4,3}(2,3) = 0.4;    % Class switch to BatchA at Router
P{4,4}(3,4) = 1;      % WebCache -> AppServer
P{4,4}(4,5) = 1;      % AppServer -> DataServer
P{4,4}(5,6) = 1;      % DataServer -> Sink

% ===== Chain 3: Real-time tasks =====
% RealTime routing:
P{5,5}(1,2) = 1;      % Source -> Router
P{5,5}(2,3) = 1;      % Router -> WebCache
P{5,5}(3,4) = 1;      % WebCache -> AppServer
P{5,5}(4,5) = 1;      % AppServer -> DataServer
P{5,5}(5,6) = 1;      % DataServer -> Sink

model.link(P);

%% Block 6: solve with multiple solvers
options = Solver.defaultOptions;
options.keep = true;
options.verbose = 1;
options.cutoff = 10;
options.seed = 23000;

disp('Open QN with 3 chains, 5 classes, and class switching:');
disp('- Chain 1: InteractiveA and InteractiveB (switch at Router)');
disp('- Chain 2: BatchA and BatchB (switch at Router)');
disp('- Chain 3: RealTime');
disp(' ');

solver = {};
solver{end+1} = MVA(model, options);
solver{end+1} = FLD(model, options);
solver{end+1} = SSA(model, options);

for s = 1:length(solver)
    fprintf('SOLVER: %s\n', solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
    fprintf('\n');
end
