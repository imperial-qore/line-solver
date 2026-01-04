clear node jobclass solver AvgTable;
%% Example of a mixed queueing network with multiple chains and class switching
% Demonstrates class switching within both open and closed job classes.
% Chain 1: Closed class with 2 classes that can switch (interactive users)
% Chain 2: Open classes with 2 classes that can switch (batch requests)
% Chain 3: Open class (external load)

model = Network('mqn_multichain');

%% Block 1: nodes
node{1} = Delay(model, 'ThinkingTime');
node{2} = Queue(model, 'WebServer', SchedStrategy.FCFS);
node{3} = Queue(model, 'AppServer', SchedStrategy.PS);
node{4} = Queue(model, 'DataServer', SchedStrategy.FCFS);
node{5} = Source(model, 'Source');
node{6} = Sink(model, 'Sink');

%% Block 2: job classes
% Chain 1: Closed classes (interactive users that can switch types)
jobclass{1} = ClosedClass(model, 'InteractiveA', 3, node{1});
jobclass{2} = ClosedClass(model, 'InteractiveB', 2, node{1});
% Chain 2: Open classes (batch jobs that can switch types)
jobclass{3} = OpenClass(model, 'BatchA', 0);
jobclass{4} = OpenClass(model, 'BatchB', 0);
% Chain 3: Open class (external load)
jobclass{5} = OpenClass(model, 'ExternalLoad', 0);

%% Block 3: service times
% ThinkingTime
node{1}.setService(jobclass{1}, Exp.fitMean(1.5));
node{1}.setService(jobclass{2}, Exp.fitMean(2.0));

% WebServer
node{2}.setService(jobclass{1}, Exp.fitMean(0.4));
node{2}.setService(jobclass{2}, Exp.fitMean(0.5));
node{2}.setService(jobclass{3}, Exp.fitMean(0.3));
node{2}.setService(jobclass{4}, Exp.fitMean(0.35));
node{2}.setService(jobclass{5}, Exp.fitMean(0.2));

% AppServer
node{3}.setService(jobclass{1}, Exp.fitMean(0.8));
node{3}.setService(jobclass{2}, Exp.fitMean(1.0));
node{3}.setService(jobclass{3}, Exp.fitMean(0.6));
node{3}.setService(jobclass{4}, Exp.fitMean(0.7));
node{3}.setService(jobclass{5}, Exp.fitMean(0.5));

% DataServer
node{4}.setService(jobclass{1}, Exp.fitMean(0.5));
node{4}.setService(jobclass{2}, Exp.fitMean(0.6));
node{4}.setService(jobclass{3}, Exp.fitMean(1.2));
node{4}.setService(jobclass{4}, Exp.fitMean(1.5));
node{4}.setService(jobclass{5}, Exp.fitMean(0.8));

%% Block 4: arrival rates
node{5}.setArrival(jobclass{3}, Exp(0.4));
node{5}.setArrival(jobclass{4}, Exp(0.2));
node{5}.setArrival(jobclass{5}, Exp(0.3));

%% Block 5: routing with class switching
M = model.getNumberOfNodes();
K = model.getNumberOfClasses();

P = model.initRoutingMatrix();

% ===== Chain 1: Closed classes with switching =====
% InteractiveA routing:
P{1,1}(1,2) = 1;      % ThinkingTime -> WebServer
P{1,1}(2,3) = 0.7;    % WebServer -> AppServer (stay in InteractiveA)
P{1,2}(2,3) = 0.3;    % Class switch to InteractiveB at WebServer
P{1,1}(3,4) = 1;      % AppServer -> DataServer
P{1,1}(4,1) = 1;      % DataServer -> ThinkingTime

% InteractiveB routing:
P{2,2}(1,2) = 1;      % ThinkingTime -> WebServer
P{2,2}(2,3) = 0.8;    % WebServer -> AppServer (stay in InteractiveB)
P{2,1}(2,3) = 0.2;    % Class switch to InteractiveA at WebServer
P{2,2}(3,4) = 1;      % AppServer -> DataServer
P{2,2}(4,1) = 1;      % DataServer -> ThinkingTime

% ===== Chain 2: Open classes with switching =====
% BatchA routing:
P{3,3}(5,2) = 1;      % Source -> WebServer
P{3,3}(2,3) = 0.6;    % WebServer -> AppServer (stay in BatchA)
P{3,4}(2,3) = 0.4;    % Class switch to BatchB at WebServer
P{3,3}(3,4) = 1;      % AppServer -> DataServer
P{3,3}(4,6) = 1;      % DataServer -> Sink

% BatchB routing:
P{4,4}(5,2) = 1;      % Source -> WebServer
P{4,4}(2,3) = 0.7;    % WebServer -> AppServer (stay in BatchB)
P{4,3}(2,3) = 0.3;    % Class switch to BatchA at WebServer
P{4,4}(3,4) = 1;      % AppServer -> DataServer
P{4,4}(4,6) = 1;      % DataServer -> Sink

% ===== Chain 3: External load =====
% ExternalLoad routing:
P{5,5}(5,2) = 1;      % Source -> WebServer
P{5,5}(2,3) = 1;      % WebServer -> AppServer
P{5,5}(3,4) = 1;      % AppServer -> DataServer
P{5,5}(4,6) = 1;      % DataServer -> Sink

model.link(P);

%% Block 6: solve with multiple solvers
options = Solver.defaultOptions;
options.keep = true;
options.verbose = 1;
options.cutoff = 5;
options.seed = 23000;

disp('Mixed QN with 3 chains, 5 classes, and class switching:');
disp('- Chain 1 (Closed): InteractiveA and InteractiveB (switch at WebServer)');
disp('- Chain 2 (Open): BatchA and BatchB (switch at WebServer)');
disp('- Chain 3 (Open): ExternalLoad');
disp(' ');

solver = {};
solver{end+1} = MVA(model, options);
solver{end+1} = SSA(model, options);

for s = 1:length(solver)
    fprintf('SOLVER: %s\n', solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgChainTable();
    AvgTable{s}
    fprintf('\n');
end
