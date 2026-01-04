%% Example: Heterogeneous Servers in LINE
%
% This example demonstrates the heterogeneous server feature in LINE,
% which allows defining queues with multiple server types having different
% service rates and class compatibilities.
%
% Features demonstrated:
% 1. Creating ServerType objects with different server counts
% 2. Setting class compatibility for server types
% 3. Setting different service rates per (class, server type) pair
% 4. Using heterogeneous scheduling policies (FSF, ALIS, etc.)

clear;

%% Example 1: Basic Heterogeneous Server Queue
%
% A queue with two server types:
% - Fast servers (2 servers, rate 2.0)
% - Slow servers (3 servers, rate 1.0)

fprintf('=== Example 1: Basic Heterogeneous Server Queue ===\n\n');

model = Network('HeteroBasic');

% Create nodes
source = Source(model, 'Source');
queue = Queue(model, 'HeteroQueue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

% Create job class
jobClass = OpenClass(model, 'Jobs');

% Set arrival rate
source.setArrival(jobClass, Exp(3.0));  % 3 jobs/sec

% Create server types
fastServers = ServerType('Fast', 2);  % 2 fast servers
slowServers = ServerType('Slow', 3);  % 3 slow servers

% Add server types to queue
queue.addServerType(fastServers);
queue.addServerType(slowServers);

% Set scheduling policy for heterogeneous servers
queue.setHeteroSchedPolicy(HeteroSchedPolicy.FSF);  % Fastest Server First

% Set service rates per server type
queue.setHeteroService(jobClass, fastServers, Exp(2.0));  % Fast: rate 2.0
queue.setHeteroService(jobClass, slowServers, Exp(1.0));  % Slow: rate 1.0

% Routing
model.link(Network.serialRouting({source, queue, sink}));

fprintf('Model configuration:\n');
fprintf('  - Arrival rate: 3.0 jobs/sec\n');
fprintf('  - Fast servers: 2 x rate 2.0 = capacity 4.0\n');
fprintf('  - Slow servers: 3 x rate 1.0 = capacity 3.0\n');
fprintf('  - Total capacity: 7.0 jobs/sec\n');
fprintf('  - Scheduling policy: FSF (Fastest Server First)\n\n');

% Solve using JMT (heterogeneous servers require simulation)
solver = SolverJMT(model, 'seed', 23000);
AvgTable1 = solver.getAvgTable()

%% Example 2: Heterogeneous Servers with Compatibility Restrictions
%
% Two job classes with different server compatibility:
% - ClassA: Can be served by both Fast and Slow servers
% - ClassB: Can only be served by Slow servers

fprintf('\n=== Example 2: Heterogeneous Servers with Compatibility ===\n\n');

model2 = Network('HeteroCompat');

% Create nodes
source2 = Source(model2, 'Source');
queue2 = Queue(model2, 'HeteroQueue', SchedStrategy.FCFS);
sink2 = Sink(model2, 'Sink');

% Create job classes
classA = OpenClass(model2, 'ClassA');
classB = OpenClass(model2, 'ClassB');

% Set arrival rates
source2.setArrival(classA, Exp(1.5));  % 1.5 jobs/sec
source2.setArrival(classB, Exp(1.0));  % 1.0 jobs/sec

% Create server types
fastServers2 = ServerType('Fast', 2);
slowServers2 = ServerType('Slow', 3);

% Set compatibility - Fast only serves ClassA, Slow serves both
fastServers2.setCompatible(classA);
slowServers2.setCompatible({classA, classB});

% Add server types to queue
queue2.addServerType(fastServers2);
queue2.addServerType(slowServers2);

% Set scheduling policy
queue2.setHeteroSchedPolicy(HeteroSchedPolicy.FSF);

% Set service rates per (class, server type)
queue2.setHeteroService(classA, fastServers2, Exp(2.0));
queue2.setHeteroService(classA, slowServers2, Exp(1.0));
queue2.setHeteroService(classB, slowServers2, Exp(1.5));

% Routing (multiclass requires cell arrays)
P = model2.initRoutingMatrix();
P{classA} = Network.serialRouting({source2, queue2, sink2});
P{classB} = Network.serialRouting({source2, queue2, sink2});
model2.link(P);

fprintf('Model configuration:\n');
fprintf('  - ClassA arrival: 1.5 jobs/sec\n');
fprintf('  - ClassB arrival: 1.0 jobs/sec\n');
fprintf('  - Fast servers: 2 (ClassA only, rate 2.0)\n');
fprintf('  - Slow servers: 3 (ClassA rate 1.0, ClassB rate 1.5)\n\n');

solver2 = SolverJMT(model2, 'seed', 23000);
AvgTable2 = solver2.getAvgTable()

fprintf('\nExpected: ClassB has higher response time due to limited server access\n');

%% Example 3: Compare Different Scheduling Policies

fprintf('\n=== Example 3: Comparing Scheduling Policies ===\n\n');

policies = {HeteroSchedPolicy.FSF, HeteroSchedPolicy.ALIS, HeteroSchedPolicy.FAIRNESS};
policyNames = {'FSF', 'ALIS', 'FAIRNESS'};

for p = 1:length(policies)
    fprintf('Policy: %s\n', policyNames{p});

    model3 = Network(sprintf('HeteroPolicy_%s', policyNames{p}));

    source3 = Source(model3, 'Source');
    queue3 = Queue(model3, 'HeteroQueue', SchedStrategy.FCFS);
    sink3 = Sink(model3, 'Sink');

    jobClass3 = OpenClass(model3, 'Jobs');
    source3.setArrival(jobClass3, Exp(2.5));

    fast3 = ServerType('Fast', 2);
    slow3 = ServerType('Slow', 2);

    queue3.addServerType(fast3);
    queue3.addServerType(slow3);
    queue3.setHeteroSchedPolicy(policies{p});

    queue3.setHeteroService(jobClass3, fast3, Exp(2.0));
    queue3.setHeteroService(jobClass3, slow3, Exp(1.0));

    model3.link(Network.serialRouting({source3, queue3, sink3}));

    solver3 = SolverJMT(model3, 'seed', 23000);
    AvgTable3 = solver3.getAvgTable();

    % Display response time from queue
    queueIdx = strcmp(AvgTable3.Station, 'HeteroQueue');
    fprintf('  Response Time: %.4f\n', AvgTable3.RespT(queueIdx));
    fprintf('  Utilization: %.4f\n\n', AvgTable3.Util(queueIdx));
end

fprintf('Examples completed!\n');
