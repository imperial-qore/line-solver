clear solver AvgTable

fprintf(1,'This example shows a compact solution of a tandem open queueing network.\n');
lambda = [1,2]/50; % lambda(r) - arrival rate of class r
D = [10,5;5,9];  % D(i,r) - mean demand of class r at station i
Z = [91,92]; % Z(r)  mean service time of class r at delay station i
model = Network.tandemPsInf(lambda,D,Z);
solver{1} = MVA(model);
AvgTable{1} = solver{1}.getAvgTable;
AvgTable{1}

% 1-line version: AvgTable=MVA(Network.tandemPsInf(A,D,Z)).getAvgTable

