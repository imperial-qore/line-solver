% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
rng(1);

%% define model
trueDemands = [0.1, 0.3];
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(model, 'Class1', 1, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 3, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{1}.setService(jobclass{2}, Exp.fitMean(1.0));
node{2}.setService(jobclass{1}, Exp.fitMean(trueDemands(1)));
node{2}.setService(jobclass{2}, Exp.fitMean(trueDemands(2)));

P = model.initRoutingMatrix;
P{1} = [0,1; 1,0];
P{2} = [0,1; 1,0];
model.link(P);

%% Get true steady-state metrics from MVA
solver_mva = SolverMVA(model);
trueRespT = solver_mva.getAvgRespT();
trueUtil = solver_mva.getAvgUtil();
trueTput = solver_mva.getAvgTput();

stIdx = node{2}.stationIndex;
trueR1 = trueRespT(stIdx, 1);
trueR2 = trueRespT(stIdx, 2);
trueU = sum(trueUtil(stIdx, :));
trueX1 = trueTput(stIdx, 1);
trueX2 = trueTput(stIdx, 2);

%% Generate noisy dataset consistent with model steady state
n = 1000;
ts = 1:n;
noise_scale = 0.05;
arvr1_samples = trueX1 * ones(n,1) + (rand(n,1)-0.5) * noise_scale * trueX1;
arvr2_samples = trueX2 * ones(n,1) + (rand(n,1)-0.5) * noise_scale * trueX2;
respt1_samples = trueR1 * ones(n,1) + (rand(n,1)-0.5) * noise_scale * trueR1;
respt2_samples = trueR2 * ones(n,1) + (rand(n,1)-0.5) * noise_scale * trueR2;
util_samples = trueU * ones(n,1) + (rand(n,1)-0.5) * noise_scale * trueU;

%% Reset service for estimation
node{2}.setService(jobclass{1}, Exp(NaN));
node{2}.setService(jobclass{2}, Exp(NaN));

%% Estimate demands
options = ParamEstimator.defaultOptions;
options.method = 'mle';
se = ParamEstimator(model, options);

lambda1 = SampledMetric(MetricType.ArvR, ts, arvr1_samples, node{2}, jobclass{1});
lambda2 = SampledMetric(MetricType.ArvR, ts, arvr2_samples, node{2}, jobclass{2});
respT1 = SampledMetric(MetricType.RespT, ts, respt1_samples, node{2}, jobclass{1});
respT2 = SampledMetric(MetricType.RespT, ts, respt2_samples, node{2}, jobclass{2});
util = SampledMetric(MetricType.Util, ts, util_samples, node{2});

se.addSamples(lambda1);
se.addSamples(lambda2);
se.addSamples(respT1);
se.addSamples(respT2);
se.addSamples(util);
se.interpolate();
estVal = se.estimateAt(node{2})

assert(all(abs(estVal - trueDemands)./trueDemands < 0.10), ...
    sprintf('MLE: estimated [%.4f, %.4f] too far from true [%.4f, %.4f]', ...
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
