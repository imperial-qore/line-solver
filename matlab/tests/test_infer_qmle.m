% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
rng(1);

%% define model with true demands
trueDemands = [0.2, 0.4];
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(model, 'Class1', 2, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 3, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{1}.setService(jobclass{2}, Exp.fitMean(1.0));

node{2}.setService(jobclass{1}, Exp.fitMean(trueDemands(1)));
node{2}.setService(jobclass{2}, Exp.fitMean(trueDemands(2)));

P = model.initRoutingMatrix;
P{1} = [0,1; 1,0];
P{2} = [0,1; 1,0];
model.link(P);

%% get true steady-state queue lengths from MVA
solver_mva = SolverMVA(model);
trueQLen = solver_mva.getAvgQLen();
trueQLen1 = trueQLen(2,1); % Queue station, Class 1
trueQLen2 = trueQLen(2,2); % Queue station, Class 2

%% generate model-consistent queue length samples
n = 5000;
ts = 1:n;
qlen1_samples = trueQLen1*ones(n,1) + rand(n,1)*0.02 - 0.01;
qlen2_samples = trueQLen2*ones(n,1) + rand(n,1)*0.02 - 0.01;

%% reset service for estimation
node{2}.setService(jobclass{1}, Exp(NaN));
node{2}.setService(jobclass{2}, Exp(NaN));

%% Estimate demands using QMLE
options = ParamEstimator.defaultOptions;
options.method = 'qmle';
se = ParamEstimator(model, options);

ql1 = SampledMetric(MetricType.QLen, ts, qlen1_samples, node{2}, jobclass{1});
ql2 = SampledMetric(MetricType.QLen, ts, qlen2_samples, node{2}, jobclass{2});

se.addSamples(ql1);
se.addSamples(ql2);
se.interpolate();
estVal = se.estimateAt(node{2})

assert(all(abs(estVal - trueDemands)./trueDemands < 0.10), ...
    sprintf('QMLE: estimated [%.4f, %.4f] too far from true [%.4f, %.4f]', ...
    estVal(1), estVal(2), trueDemands(1), trueDemands(2)));
fprintf(1, 'Estimated demands: Class1=%.4f, Class2=%.4f\n', estVal(1), estVal(2));

%% Solve model
solver = {};
solver{end+1} = SolverMVA(model);

AvgTable = cell(1,length(solver));
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
