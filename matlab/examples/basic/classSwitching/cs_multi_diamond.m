clear node jobclass solver AvgTable;
%% Example of class switching controlled by a reducible Markov chain
model = Network('mm1cs');

%% Block 1: nodes
node{1} = Source(model, 'Source 1');
node{2} = Queue(model, 'Queue 0', SchedStrategy.FCFS);
node{3} = Queue(model, 'Queue 1', SchedStrategy.FCFS);
node{4} = Queue(model, 'Queue 2', SchedStrategy.FCFS);
node{5} = Sink(model, 'Sink 1');

%% Block 2: classes
jobclass{1} = OpenClass(model, 'Class1', 0);
jobclass{2} = OpenClass(model, 'Class2', 0);
jobclass{3} = OpenClass(model, 'Class3', 0);

node{1}.setArrival(jobclass{1}, Exp(1.000000)); % (Source 1,Class1)
node{2}.setService(jobclass{1}, Exp(10.000000)); % (Queue 1,Class2)
node{3}.setService(jobclass{2}, Exp(20.000000)); % (Queue 1,Class2)
node{4}.setService(jobclass{3}, Exp(30.000000)); % (Queue 2,Class3)

P = model.initRoutingMatrix(); % initialize routing matrix 
P{1,1}(1,2) = 1;
P{1,1}(2,2) = 0.2;
P{1,2}(2,3) = 0.3;
P{1,3}(2,4) = 0.5;
P{2,2}(3,5) = 1;
P{3,3}(4,5) = 1;
model.link(P);

model.printRoutingMatrix();

solver{1} = MVA(model);
AvgTable{1} = solver{1}.getAvgChainTable;
AvgTable{1}
