%% Finite Capacity Region with Blocking vs M/M/1
% This example shows that an FCR with blocking (waitq) around a single
% queue behaves identically to a standard M/M/1 queue, since jobs
% simply wait when the region is "full" (no actual capacity limit effect).

clear; clc;

%% Parameters
arrivalRate = 0.5;
serviceRate = 1.0;

%% Model 1: Queue with FCR (blocking)
model1 = Network('FCR Blocking');

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

% Add FCR with blocking - jobs wait when region is full
fcr = model1.addRegion({queue1});
fcr.setGlobalMaxJobs(10);
fcr.setDropRule(jobclass1, false);  % false = block (wait)

%% Model 2: Standard M/M/1 (no FCR)
model2 = Network('M/M/1');

source2 = Source(model2, 'Source');
queue2 = Queue(model2, 'Queue', SchedStrategy.FCFS);
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
avgTable1 = solver1.getAvgNodeTable()

solver2 = JMT(model2, 'seed', 23000, 'samples', 100000);
avgTable2 = solver2.getAvgNodeTable()
