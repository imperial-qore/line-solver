clear node jobclass solver AvgTable;
%% Example of a closed queueing network with 3 chains, 5 classes, and class switching
% Product-form (BCMP) network with class switching.
% Chain 1: High-priority users (2 classes that can switch)
% Chain 2: Regular users (2 classes that can switch)
% Chain 3: Background tasks (1 class)
%
% The network maintains product-form property:
% - All service times are exponential
% - Mix of FCFS and PS scheduling
% - All nodes are conventional queueing stations or infinite servers

model = Network('cqn_multichain');

%% Block 1: nodes
node{1} = Delay(model, 'ThinkingTime');
node{2} = Queue(model, 'WebServer', SchedStrategy.FCFS);
node{3} = Queue(model, 'AppServer', SchedStrategy.PS);
node{4} = Queue(model, 'DataServer', SchedStrategy.FCFS);

%% Block 2: job classes
% Chain 1: High-priority users (can switch between HighPriorityA and HighPriorityB)
jobclass{1} = ClosedClass(model, 'HighPriorityA', 2, node{1});
jobclass{2} = ClosedClass(model, 'HighPriorityB', 1, node{1});
% Chain 2: Regular users (can switch between RegularA and RegularB)
jobclass{3} = ClosedClass(model, 'RegularA', 3, node{1});
jobclass{4} = ClosedClass(model, 'RegularB', 2, node{1});
% Chain 3: Background tasks
jobclass{5} = ClosedClass(model, 'Background', 2, node{1});

%% Block 3: service times (all exponential for product-form)
% ThinkingTime (infinite server node)
node{1}.setService(jobclass{1}, Exp.fitMean(1.5));
node{1}.setService(jobclass{2}, Exp.fitMean(1.6));
node{1}.setService(jobclass{3}, Exp.fitMean(2.0));
node{1}.setService(jobclass{4}, Exp.fitMean(2.1));
node{1}.setService(jobclass{5}, Exp.fitMean(3.0));

% WebServer (FCFS, exponential)
node{2}.setService(jobclass{1}, Exp.fitMean(0.4));
node{2}.setService(jobclass{2}, Exp.fitMean(0.42));
node{2}.setService(jobclass{3}, Exp.fitMean(0.5));
node{2}.setService(jobclass{4}, Exp.fitMean(0.52));
node{2}.setService(jobclass{5}, Exp.fitMean(0.8));

% AppServer (PS, exponential)
node{3}.setService(jobclass{1}, Exp.fitMean(0.7));
node{3}.setService(jobclass{2}, Exp.fitMean(0.75));
node{3}.setService(jobclass{3}, Exp.fitMean(0.9));
node{3}.setService(jobclass{4}, Exp.fitMean(0.95));
node{3}.setService(jobclass{5}, Exp.fitMean(1.2));

% DataServer (FCFS, exponential)
node{4}.setService(jobclass{1}, Exp.fitMean(0.5));
node{4}.setService(jobclass{2}, Exp.fitMean(0.52));
node{4}.setService(jobclass{3}, Exp.fitMean(0.6));
node{4}.setService(jobclass{4}, Exp.fitMean(0.62));
node{4}.setService(jobclass{5}, Exp.fitMean(1.0));

%% Block 4: routing with class switching
M = model.getNumberOfNodes();
K = model.getNumberOfClasses();

P = model.initRoutingMatrix();

% ===== Chain 1: High-priority users with switching =====
% HighPriorityA routing:
P{1,1}(1,2) = 1;      % ThinkingTime -> WebServer
P{1,1}(2,3) = 0.7;    % WebServer -> AppServer (stay in HighPriorityA)
P{1,2}(2,3) = 0.3;    % Class switch to HighPriorityB at WebServer
P{1,1}(3,4) = 1;      % AppServer -> DataServer
P{1,1}(4,1) = 1;      % DataServer -> ThinkingTime

% HighPriorityB routing:
P{2,2}(1,2) = 1;      % ThinkingTime -> WebServer
P{2,2}(2,3) = 0.75;   % WebServer -> AppServer (stay in HighPriorityB)
P{2,1}(2,3) = 0.25;   % Class switch to HighPriorityA at WebServer
P{2,2}(3,4) = 1;      % AppServer -> DataServer
P{2,2}(4,1) = 1;      % DataServer -> ThinkingTime

% ===== Chain 2: Regular users with switching =====
% RegularA routing:
P{3,3}(1,2) = 1;      % ThinkingTime -> WebServer
P{3,3}(2,3) = 0.65;   % WebServer -> AppServer (stay in RegularA)
P{3,4}(2,3) = 0.35;   % Class switch to RegularB at WebServer
P{3,3}(3,4) = 1;      % AppServer -> DataServer
P{3,3}(4,1) = 1;      % DataServer -> ThinkingTime

% RegularB routing:
P{4,4}(1,2) = 1;      % ThinkingTime -> WebServer
P{4,4}(2,3) = 0.7;    % WebServer -> AppServer (stay in RegularB)
P{4,3}(2,3) = 0.3;    % Class switch to RegularA at WebServer
P{4,4}(3,4) = 1;      % AppServer -> DataServer
P{4,4}(4,1) = 1;      % DataServer -> ThinkingTime

% ===== Chain 3: Background tasks =====
% Background routing:
P{5,5}(1,2) = 1;      % ThinkingTime -> WebServer
P{5,5}(2,3) = 1;      % WebServer -> AppServer
P{5,5}(3,4) = 1;      % AppServer -> DataServer
P{5,5}(4,1) = 1;      % DataServer -> ThinkingTime

model.link(P);

%% Block 5: solve with multiple solvers
options = Solver.defaultOptions;
options.keep = true;
options.verbose = 1;

disp('Closed QN with 3 chains, 5 classes, and class switching (product-form):');
disp('- Chain 1 (N=3): HighPriorityA and HighPriorityB (switch at WebServer)');
disp('- Chain 2 (N=5): RegularA and RegularB (switch at WebServer)');
disp('- Chain 3 (N=2): Background');
disp('Network maintains BCMP product-form property (exponential, FCFS/PS)');
disp(' ');

solver = {};
solver{end+1} = MVA(model, options);
solver{end+1} = CTMC(model, options);
solver{end+1} = SSA(model, options);

for s = 1:length(solver)
    fprintf('SOLVER: %s\n', solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgChainTable();
    AvgTable{s}
    fprintf('\n');
end
