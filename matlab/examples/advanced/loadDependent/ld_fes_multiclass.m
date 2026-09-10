clear node jobclass

% Multiclass FES Aggregation with Norton's Theorem Verification
%
% This example demonstrates Norton's theorem for closed multiclass
% queueing networks: a subset of stations is replaced by a single
% Flow-Equivalent Server (FES) whose class-dependent service rates equal the
% throughputs of the isolated subnetwork.

N1 = 3; % class-1 jobs
N2 = 2; % class-2 jobs

%% Original 4-station tandem network with 2 classes
model = Network('OriginalModel');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Q1', SchedStrategy.PS);
node{3} = Queue(model, 'Q2', SchedStrategy.PS);
node{4} = Queue(model, 'Q3', SchedStrategy.PS);

jobclass{1} = ClosedClass(model, 'Class1', N1, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', N2, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{1}.setService(jobclass{2}, Exp.fitMean(1.5));
node{2}.setService(jobclass{1}, Exp.fitMean(0.5));
node{2}.setService(jobclass{2}, Exp.fitMean(0.8));
node{3}.setService(jobclass{1}, Exp.fitMean(0.3));
node{3}.setService(jobclass{2}, Exp.fitMean(0.6));
node{4}.setService(jobclass{1}, Exp.fitMean(0.4));
node{4}.setService(jobclass{2}, Exp.fitMean(0.7));

P = model.initRoutingMatrix();
P{1,1} = [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0];
P{2,2} = P{1,1};
model.link(P);

%% Solve original model
fprintf('MVA (original):\n');
AvgOrig = MVA(model, 'method', 'exact').getAvgTable;
disp(AvgOrig);

%% Aggregate Q1, Q2, Q3 into FES and solve with NC (convolution)
[fesModel, ~, ~] = ModelAdapter.aggregateFES(model, ...
    {node{2}, node{3}, node{4}}, struct('solver','mva','verbose',false));

fprintf('NC (FES model):\n');
AvgNC = NC(fesModel, 'method', 'exact').getAvgTable;
disp(AvgNC);
