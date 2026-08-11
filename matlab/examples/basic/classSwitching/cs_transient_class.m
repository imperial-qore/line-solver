clear node jobclass solver AvgTable;
%% Example of class switching controlled by a reducible Markov chain
% In this variant the job remains either in class 2 or class 3 forever
%
% This example demonstrates a transient class model where Class1 can
% switch to Class2 or Class3, but those classes never return to Class1.
% This creates a reducible routing chain with transient states (Class1)
% and recurrent states (Class2/Class3 cycles).

model = Network('mm1cs');

%% Block 1: nodes
node{1} = Delay(model, 'Queue 0');
node{2} = Delay(model, 'Queue 1');
node{3} = Delay(model, 'Queue 2');

%% Block 2: classes
jobclass{1} = ClosedClass(model, 'Class1', 1, node{1});
jobclass{2} = ClosedClass(model, 'Class2', 0, node{1});
jobclass{3} = ClosedClass(model, 'Class3', 0, node{1});

node{1}.setService(jobclass{1}, Exp.fitMean(1.000000)); % (Queue 1,Class1)
node{1}.setService(jobclass{2}, Exp.fitMean(1.000000)); % (Queue 1,Class2)
node{1}.setService(jobclass{3}, Exp.fitMean(1.000000)); % (Queue 1,Class3)
node{2}.setService(jobclass{2}, Exp.fitMean(1.000000)); % (Queue 1,Class2)
node{3}.setService(jobclass{3}, Exp.fitMean(1.000000)); % (Queue 2,Class3)

P = model.initRoutingMatrix(); % initialize routing matrix 
P{1,1}(1,1) = 0.2;
P{2,2}(1,2) = 1.0;
P{3,3}(1,3) = 1.0;
P{1,2}(1,2) = 0.3;
P{1,3}(1,3) = 0.5;
P{2,2}(2,1) = 1;
P{3,3}(3,1) = 1;
model.link(P);

model.printRoutingMatrix();

solver{1} = MVA(model);
AvgTable{1} = solver{1}.getAvgChainTable;
AvgTable{1}
