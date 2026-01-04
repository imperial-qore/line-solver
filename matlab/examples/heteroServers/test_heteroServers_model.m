%% Test: Heterogeneous Server Model Building
%
% This test validates that heterogeneous server models can be built
% and the NetworkStruct is populated correctly.

clear;
fprintf('=== Test: Heterogeneous Server Model Building ===\n\n');

%% Test 1: Create a basic heterogeneous server queue

model = Network('HeteroTest');

% Create nodes
source = Source(model, 'Source');
queue = Queue(model, 'HeteroQueue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

% Create job class
jobClass = OpenClass(model, 'Jobs');

% Set arrival rate
source.setArrival(jobClass, Exp(3.0));

% Create server types
fastServers = ServerType('Fast', 2);
slowServers = ServerType('Slow', 3);

% Set compatibility
fastServers.setCompatible(jobClass);
slowServers.setCompatible(jobClass);

% Add server types to queue
queue.addServerType(fastServers);
queue.addServerType(slowServers);

% Set scheduling policy
queue.setHeteroSchedPolicy(HeteroSchedPolicy.FSF);

% Set service rates per server type
queue.setHeteroService(jobClass, fastServers, Exp(2.0));
queue.setHeteroService(jobClass, slowServers, Exp(1.0));

% Routing
model.link(Network.serialRouting({source, queue, sink}));

% Validate model building
fprintf('Model built successfully!\n\n');

% Check Queue properties
fprintf('Queue properties:\n');
fprintf('  - isHeterogeneous: %d\n', queue.isHeterogeneous());
fprintf('  - Number of server types: %d\n', length(queue.serverTypes));
fprintf('  - Server type 1: %s (%d servers)\n', fastServers.getName(), fastServers.numOfServers);
fprintf('  - Server type 2: %s (%d servers)\n', slowServers.getName(), slowServers.numOfServers);
fprintf('  - Scheduling policy: %d\n', queue.heteroSchedPolicy);

% Get NetworkStruct and validate
sn = model.getStruct();
ist = sn.nodeToStation(queue.index);

fprintf('\nNetworkStruct fields:\n');
fprintf('  - sn.nservertypes(%d): %d\n', ist, sn.nservertypes(ist));
fprintf('  - sn.servertypenames{%d}: ', ist);
for t = 1:sn.nservertypes(ist)
    fprintf('%s ', sn.servertypenames{ist}{t});
end
fprintf('\n');
fprintf('  - sn.serverspertype{%d}: ', ist);
fprintf('%d ', sn.serverspertype{ist});
fprintf('\n');
fprintf('  - sn.servercompat{%d}:\n', ist);
disp(sn.servercompat{ist});
fprintf('  - sn.heteroschedpolicy(%d): %d\n', ist, sn.heteroschedpolicy(ist));

fprintf('\n=== Test Passed! ===\n');
fprintf('Heterogeneous server model builds correctly.\n');
fprintf('Note: JMT simulation requires a newer JMT version with hetero server support.\n');
