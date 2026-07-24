fprintf(1,'This example shows a 1-line exact MVA solution of a cyclic queueing network.\n');
D = [10,5; 5,9]; % S(i,r) - mean service time of class r at station i
N = [1,2]; % N(r) - number of jobs of class r
Z = [91,92; 93,94]; % Z(r)  mean service time of class r at delay station i

model = Network.cyclicPsInf(N,D,Z);
avgTable = MVA(model, 'exact').getAvgTable
