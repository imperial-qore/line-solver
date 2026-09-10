% SF_CLOSED  Closed cluster: 8 jobs cycling between Think and 3 PS servers.

clear all

N = 8;
Z = 1.0;
D = ones(3, 1);
strategy = {SchedStrategy.PS, SchedStrategy.PS, SchedStrategy.PS};

model = Network.clusterClosed(N, Z, D, strategy, ones(3, 1), RoutingStrategy.RAND);

solver = MVA(model);
disp(solver.getAvgTable());
