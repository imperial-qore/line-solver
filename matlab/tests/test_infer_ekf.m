% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
rng(1);

%% define model
trueDemand = 0.3;
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(model, 'Class1', 5, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{2}.setService(jobclass{1}, Exp.fitMean(trueDemand));

P = model.initRoutingMatrix;
P{1} = [0,1; 1,0];
model.link(P);

%% Get true steady-state metrics from MVA
solver_mva = SolverMVA(model);
trueRespT = solver_mva.getAvgRespT();
trueUtil = solver_mva.getAvgUtil();
trueTput = solver_mva.getAvgTput();

stIdx = node{2}.stationIndex;
trueR = trueRespT(stIdx, 1);
trueU = trueUtil(stIdx, 1);
trueX = trueTput(stIdx, 1);

%% Generate noisy dataset consistent with model steady state
n = 1000;
ts = 1:n;
noise_scale = 0.05;
arvr_samples = trueX * ones(n,1) + (rand(n,1)-0.5) * noise_scale * trueX;
respt_samples = trueR * ones(n,1) + (rand(n,1)-0.5) * noise_scale * trueR;
util_samples = trueU * ones(n,1) + (rand(n,1)-0.5) * noise_scale * trueU;

%% Reset service for estimation
node{2}.setService(jobclass{1}, Exp(NaN));

%% Estimate demands
options = ParamEstimator.defaultOptions;
options.method = 'ekf';
se = ParamEstimator(model, options);

lambda1 = SampledMetric(MetricType.ArvR, ts, arvr_samples, node{2}, jobclass{1});
respT1 = SampledMetric(MetricType.RespT, ts, respt_samples, node{2}, jobclass{1});
util = SampledMetric(MetricType.Util, ts, util_samples, node{2});
se.addSamples(lambda1);
se.addSamples(respT1);
se.addSamples(util);

estVal = se.estimateAt(node{2})

assert(abs(estVal - trueDemand)/trueDemand < 0.10, ...
    sprintf('EKF: estimated %.4f, relative error %.2f%% exceeds 10%%', estVal, 100*abs(estVal - trueDemand)/trueDemand));

%% Solve model
solver = {};
solver{end+1} = SolverMVA(model);

AvgTable = cell(1,length(solver));
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
