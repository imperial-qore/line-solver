% SF_COMPARE  Compare RAND vs RROBIN dispatching on the same cluster.
%
% MVA only handles RAND routing (product-form). For non-product-form
% policies we use SSA (a discrete-event simulator) with a low sample
% budget for a quick comparison.

clear all

cluster = Cluster().setNumStations(4).setArrivalRate(1.0).setServiceRate(0.4);
cluster.setScheduling(SchedStrategy.PS);

solverFcn = @(m) SSA(m, 'seed', 23000, 'samples', 2000).getAvgTable();

policies = [RoutingStrategy.RAND, RoutingStrategy.RROBIN];
results = cluster.compareDispatching(solverFcn, policies);

keys_ = results.keys;
for k = 1:numel(keys_)
    fprintf('\n=== Dispatching: %s ===\n', keys_{k});
    disp(results(keys_{k}));
end
