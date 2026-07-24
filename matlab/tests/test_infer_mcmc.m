% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
rng(1);

%% define model with true demand within MCMC integral range [0, 0.2]
trueDemand = 0.1;
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(model, 'Class1', 1, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{2}.setService(jobclass{1}, Exp.fitMean(trueDemand));

P = model.initRoutingMatrix;
P{1} = [0,1; 1,0];
model.link(P);

%% get true steady-state queue length from MVA
solver_mva = SolverMVA(model);
trueQLen = solver_mva.getAvgQLen();
trueQLen_queue = trueQLen(2); % Queue station

%% generate model-consistent queue length samples
n = 5000;
ts = 1:n;
qlen_samples = trueQLen_queue * ones(n,1) + rand(n,1)*0.005 - 0.0025;

%% reset service for estimation
node{2}.setService(jobclass{1}, Exp(NaN));

%% Estimate demands using MCMC
options = ParamEstimator.defaultOptions;
options.method = 'mcmc';
se = ParamEstimator(model, options);

ql = SampledMetric(MetricType.QLen, ts, qlen_samples, node{2}); % aggregate queue-length

se.addSamples(ql);
se.interpolate();
estVal = se.estimateAt(node{2})

assert(abs(estVal - trueDemand)/trueDemand < 0.10, ...
    sprintf('MCMC: estimated %.4f, relative error %.2f%% exceeds 10%%', estVal, 100*abs(estVal - trueDemand)/trueDemand));
fprintf(1, 'Estimated demand: Class1=%.4f\n', estVal(1));

%% Solve model
solver = {};
solver{end+1} = SolverMVA(model);

AvgTable = cell(1,length(solver));
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
