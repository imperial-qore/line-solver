%% Linear Admission Constraints for Finite Capacity Regions
% Open network with a Delay node inside an FCR with DROP policy and
% general linear admission constraints: A*n <= b.
%
% Constraint: 2*n_1 + 1*n_2 <= 5, 1*n_1 + 3*n_2 <= 7
% (general weighted cross-class coupling)

clear; clc;

model = Network('FCR LinCon');

source = Source(model, 'Source');
delay = Delay(model, 'Delay');
sink = Sink(model, 'Sink');

class1 = OpenClass(model, 'Class1', 0);
class2 = OpenClass(model, 'Class2', 0);

source.setArrival(class1, Exp(0.3));
source.setArrival(class2, Exp(0.2));
delay.setService(class1, Exp(1.0));
delay.setService(class2, Exp(0.8));

P = model.initRoutingMatrix();
P.set(class1, class1, source, delay, 1.0);
P.set(class1, class1, delay, sink, 1.0);
P.set(class2, class2, source, delay, 1.0);
P.set(class2, class2, delay, sink, 1.0);
model.link(P);

% Add FCR with DROP policy and linear constraints
fcr = model.addRegion({delay});
fcr.setClassMaxJobs(class1, 1000);
fcr.setClassMaxJobs(class2, 1000);
fcr.setDropRule(class1, true);
fcr.setDropRule(class2, true);
fcr.setGlobalMaxJobs(1000);
fcr.setLinearConstraints([2 1; 1 3], [5; 7]);

solver = LDES(model, 'seed', 23000, 'samples', 500000);
avgTable = solver.getAvgTable()
