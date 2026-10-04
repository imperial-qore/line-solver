%% Finite Capacity Region - Constraint Types
% This example demonstrates different FCR constraint types:
% - Global max jobs: limits total jobs across all classes
% - Per-class max jobs: limits jobs of a specific class
% - Drop rule: determines blocking vs dropping behavior

clear; clc;

%% Create a multiclass open network
model = Network('FCR Constraints Demo');

% Nodes
source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

% Job classes
class1 = OpenClass(model, 'HighPriority', 0);
class2 = OpenClass(model, 'LowPriority', 1);

% Arrival and service rates
source.setArrival(class1, Exp(0.3));  % High priority: lower arrival rate
source.setArrival(class2, Exp(0.5));  % Low priority: higher arrival rate
queue1.setService(class1, Exp(1.0));
queue1.setService(class2, Exp(0.8));
queue2.setService(class1, Exp(1.2));
queue2.setService(class2, Exp(1.0));

% Routing: both classes go through both queues
P = model.initRoutingMatrix();
P.set(class1, class1, source, queue1, 1.0);
P.set(class1, class1, queue1, queue2, 1.0);
P.set(class1, class1, queue2, sink, 1.0);
P.set(class2, class2, source, queue1, 1.0);
P.set(class2, class2, queue1, queue2, 1.0);
P.set(class2, class2, queue2, sink, 1.0);
model.link(P);

%% Add Finite Capacity Region with multiple constraints
fcr = model.addRegion({queue1, queue2});

% Global constraint: max 6 jobs total in the region
fcr.setGlobalMaxJobs(2);

% Per-class constraints: high priority gets more space
fcr.setClassMaxJobs(class1, 2);  % HighPriority: max 1 jobs
fcr.setClassMaxJobs(class2, 2);  % LowPriority: max 2 jobs
% Note: per-class limits must sum to >= global limit for consistent behavior

% Drop rules: all classes use the same rule (required by LINE)
% true = drop jobs when limit reached, false = block (wait)
fcr.setDropRule(class1, true);   % true = drop
fcr.setDropRule(class2, true);   % true = drop

%% Solve with JMT
solver = JMT(model, 'seed', 23000, 'samples', 100000);
avgTable = solver.getAvgTable()
