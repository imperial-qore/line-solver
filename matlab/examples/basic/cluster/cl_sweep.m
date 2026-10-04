% SF_SWEEP  Arrival-rate sweep on a 2-PS-cluster. Shows response time
% growing as lambda approaches saturation.

clear all

cluster = Cluster().setNumStations(2).setArrivalRate(0.1).setServiceRate(1.0);
cluster.setScheduling(SchedStrategy.PS).setDispatching(RoutingStrategy.RAND);

solverFcn = @(m) MVA(m).getAvgTable();
results = cluster.sweepArrivalRate([0.2, 0.5, 0.9, 1.5], solverFcn);

keys_ = sort(keys(results));
for k = 1:numel(keys_)
    fprintf('\n=== lambda = %g ===\n', keys_(k));
    disp(results{keys_(k)});
end
