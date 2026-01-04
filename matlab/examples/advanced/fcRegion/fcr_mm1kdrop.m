%% Finite Capacity Region with Dropping vs M/M/1/K
% This example shows that an FCR with dropping around a single queue
% behaves like an M/M/1/K queue, where K is the FCR capacity.
% Jobs arriving when the region is full are dropped (lost).

clear; clc;

%% Parameters
arrivalRate = 0.8;
serviceRate = 1.0;
K = 3;

%% Model 1: Queue with FCR (dropping)
model1 = Network('FCR Dropping');

source1 = Source(model1, 'Source');
queue1 = Queue(model1, 'Queue', SchedStrategy.FCFS);
sink1 = Sink(model1, 'Sink');

jobclass1 = OpenClass(model1, 'Class1', 0);

source1.setArrival(jobclass1, Exp(arrivalRate));
queue1.setService(jobclass1, Exp(serviceRate));

P1 = model1.initRoutingMatrix();
P1.set(jobclass1, jobclass1, source1, queue1, 1.0);
P1.set(jobclass1, jobclass1, queue1, sink1, 1.0);
model1.link(P1);

% Add FCR with dropping - jobs are lost when region is full
fcr = model1.addRegion({queue1});
fcr.setGlobalMaxJobs(K);
fcr.setDropRule(jobclass1, true);  % true = drop jobs

%% Model 2: M/M/1/K using queue capacity
model2 = Network('M/M/1/K');

source2 = Source(model2, 'Source');
queue2 = Queue(model2, 'Queue', SchedStrategy.FCFS);
queue2.setNumberOfServers(1);
queue2.setCapacity(K);  % Set queue capacity to K
sink2 = Sink(model2, 'Sink');

jobclass2 = OpenClass(model2, 'Class1', 0);

source2.setArrival(jobclass2, Exp(arrivalRate));
queue2.setService(jobclass2, Exp(serviceRate));

P2 = model2.initRoutingMatrix();
P2.set(jobclass2, jobclass2, source2, queue2, 1.0);
P2.set(jobclass2, jobclass2, queue2, sink2, 1.0);
model2.link(P2);

%% Solve both models and compare results
solver1 = JMT(model1, 'seed', 23000, 'samples', 100000);
avgTable1 = solver1.getAvgTable()

solver2 = JMT(model2, 'seed', 23000, 'samples', 100000);
avgTable2 = solver2.getAvgTable()
