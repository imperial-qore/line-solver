% DISPATCH_MIXED  Parity test: mixed cluster built via the Cluster builder.
%
% 2 PS server stations shared by one open class (lambda = 0.5) and one
% closed class (N = 3 jobs, think time Z = 1.0), mu = 2.0 for the open
% class and mu = 1.5 for the closed one, RAND dispatching. Solved by MVA
% and re-solved by simulation via LDES for cross-validation. The same
% script runs unchanged under PYTHON, MATLAB, and MATLAB-WRAPPER.

clear all

cluster = Cluster();
cluster.setNumStations(2);
cluster.setMixed(0.5, 3, 1.0);
cluster.setServiceRates([2.0 1.5; 2.0 1.5]);
cluster.setScheduling(SchedStrategy.PS);

avgTable = MVA(cluster.build()).getAvgTable()
ldesTable = LDES(cluster.build(), 'seed', 23000).getAvgTable()
