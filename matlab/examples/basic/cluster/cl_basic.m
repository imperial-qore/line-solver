% SF_BASIC  Open cluster with four PS servers and random dispatching.
%
% Source -> Dispatcher -> Server[1..4] -> Sink. The single call to
% Network.clusterPs collapses what would otherwise be ~20 lines of
% manual node/link wiring.

clear all

lambda = 0.4;
D = ones(4, 1);                  % mean service time = 1 at each server
model = Network.clusterPs(lambda, D, RoutingStrategy.RAND);

solver = MVA(model);
disp(solver.getAvgTable());
