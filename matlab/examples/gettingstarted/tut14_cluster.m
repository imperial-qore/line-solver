% Example 13: Open cluster
%
% Source -> Dispatcher (Router) -> Server[1..M] -> Sink, with a single open
% class. Two equivalent ways to build the model are shown:
%   (a) the one-liner Network.clusterPs static factory;
%   (b) the chainable Cluster builder, which also exposes helpers to
%       compare dispatching policies on the same cluster.

%% Block 1: one-liner factory
lambda = 0.4;          % arrival rate of the open class
D = ones(3, 1);        % mean service time = 1 at each of the 3 servers
model = Network.clusterPs(lambda, D, RoutingStrategy.RAND);

avgTable = MVA(model).getAvgTable()

%% Block 2: Cluster builder with non-uniform multi-server queues
cluster = Cluster().setNumStations(3).setArrivalRate(0.4).setServiceRate(1.0);
cluster.setScheduling(SchedStrategy.FCFS);
cluster.setStationServers([2; 1; 1]);   % Server1 is M/M/2, Server2/3 are M/M/1
avgTableFcfs = MVA(cluster.build()).getAvgTable()

%% Block 3: cross-check the same FCFS multi-server model under three simulators
%  JMT is a Java-based discrete-event simulator that uses an XML model file;
%  LDES is the LINE Discrete Event Simulator built on the SSJ library and
%  invoked as a subprocess; SSA is LINE's native stochastic simulator
%  using the next-reaction method. All three produce statistically
%  equivalent results on this open-class cluster.
fprintf('\n=== JMT ===\n');
disp(JMT(cluster.build(), 'seed', 23000, 'samples', 20000).getAvgTable());
fprintf('=== LDES ===\n');
disp(LDES(cluster.build(), 'seed', 23000, 'samples', 20000).getAvgTable());
fprintf('=== SSA ===\n');
disp(SSA(cluster.build(), 'seed', 23000, 'samples', 20000).getAvgTable());

%% Block 4: compare dispatching policies via simulation
%  MVA assumes RAND (product-form). For non-product-form policies such as
%  RROBIN we drop to a simulator with a small sample budget.
cluster2 = Cluster().setNumStations(3).setArrivalRate(0.4).setServiceRate(1.0);
cluster2.setScheduling(SchedStrategy.PS);
solverFcn = @(m) JMT(m, 'seed', 23000, 'samples', 5000).getAvgTable();
results = cluster2.compareDispatching(solverFcn, ...
    [RoutingStrategy.RAND, RoutingStrategy.RROBIN]);

keys_ = results.keys;
for k = 1:numel(keys_)
    fprintf('\n=== Dispatching: %s ===\n', keys_{k});
    disp(results(keys_{k}));
end
