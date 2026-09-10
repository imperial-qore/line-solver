% DISPATCH_OPEN  Parity test: open cluster built via the Cluster builder.
%
% 3 PS server stations, single open class, lambda = 0.4, mu = 1.0, RAND
% dispatching. Solved exactly via MVA, then re-solved by stochastic
% simulation via SSA for cross-validation. The same script runs unchanged
% under PYTHON, PYTHON-WRAPPER, MATLAB, and MATLAB-WRAPPER.

clear all

cluster = Cluster();
cluster.setNumStations(3);
cluster.setArrivalRate(0.4);
cluster.setServiceRate(1.0);
cluster.setScheduling(SchedStrategy.PS);

model = cluster.build();
avgTable = MVA(model).getAvgTable()
ssaTable = SSA(model, 'seed', 23000, 'samples', 20000).getAvgTable()
