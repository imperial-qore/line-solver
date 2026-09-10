% CL_STATIONS  Horizontal-scaling sweep: response time against the number of
% stations at a fixed total arrival rate.
%
% Cluster.sweepNumStations replicates the single-station service rate across M
% stations and returns a dictionary keyed by M. Adding stations splits the same
% arrival stream, so per-station utilization falls and response time drops
% towards the bare service time.

clear all

cluster = Cluster().setNumStations(1).setArrivalRate(1.6).setServiceRate(1.0);
cluster.setScheduling(SchedStrategy.PS).setDispatching(RoutingStrategy.RAND);

solverFcn = @(m) MVA(m).getAvgTable();
results = cluster.sweepNumStations([2, 3, 4, 6], solverFcn);

counts = sort(keys(results));
for k = 1:numel(counts)
    T = results{counts(k)};
    % getAvgTable returns an IndexedTable; getTable unwraps the MATLAB table
    tbl = T.getTable();
    isStation = ~strcmp(string(tbl.Station), "Source");
    fprintf('\n=== M = %d  (mean RespT %.4f, max Util %.4f) ===\n', ...
        counts(k), mean(tbl.RespT(isStation)), max(tbl.Util(isStation)));
    disp(T);
end
