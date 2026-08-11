% DISPATCH_CLOSED  Parity test: closed cluster built via the Cluster builder.
%
% 3 PS server stations, single closed class with N = 10 jobs, think time
% Z = 1.0, mu = 1.0, RAND dispatching. Solved exactly via MVA, then
% re-solved by stochastic simulation via SSA for cross-validation. The
% same script runs unchanged under PYTHON, PYTHON-WRAPPER, MATLAB, and
% MATLAB-WRAPPER.

clear all

cluster = Cluster();
cluster.setNumStations(3);
cluster.setServiceRate(1.0);
cluster.setScheduling(SchedStrategy.PS);
cluster.setClosed(10, 1.0);

avgTable = MVA(cluster.build()).getAvgTable()
ssaTable = SSA(cluster.build(), 'seed', 23000, 'samples', 20000).getAvgTable()
