% SF_MULTICLASS  Two-class open cluster (e.g. interactive vs batch).

clear all

lambda = [0.3, 0.2];
D = [1.0, 0.5; 1.0, 0.5];   % rows = servers, cols = classes
model = Network.clusterPs(lambda, D, RoutingStrategy.RAND);

solver = MVA(model);
disp(solver.getAvgTable());
