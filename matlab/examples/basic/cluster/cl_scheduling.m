% CL_SCHEDULING  Compare scheduling disciplines on the same cluster.
%
% Cluster.compareScheduling rebuilds the model once per discipline and returns a
% dictionary from the discipline to its AvgTable, restoring the original setting
% on exit. With exponential service FCFS and PS agree in the mean; the contrast
% appears once the service SCV departs from 1.

clear all

cluster = Cluster().setNumStations(3).setArrivalRate(0.9).setServiceRate(0.5);
cluster.setDispatching(RoutingStrategy.RAND);
cluster.setServiceSCV(4.0);

solverFcn = @(m) SSA(m, 'seed', 23000, 'samples', 20000).getAvgTable();

disciplines = [SchedStrategy.FCFS, SchedStrategy.PS];
results = cluster.compareScheduling(solverFcn, disciplines);

% compareScheduling keys by the numeric enum value, so map it back to a name
keys_ = keys(results);
for k = 1:numel(keys_)
    fprintf('\n=== Scheduling: %s ===\n', SchedStrategy.toText(str2double(keys_(k))));
    disp(results{keys_(k)});
end
