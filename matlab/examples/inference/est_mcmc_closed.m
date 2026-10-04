clear node jobclass solver AvgTable
%% Example: MCMC estimation on a closed network
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% define model
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(model, 'Class1', 2, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 3, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{1}.setService(jobclass{2}, Exp.fitMean(1.0));

node{2}.setService(jobclass{1}, Exp(NaN));  % NaN = to be estimated
node{2}.setService(jobclass{2}, Exp(NaN));  % NaN = to be estimated

P = model.initRoutingMatrix;
P{1} = [0,1; 1,0];
P{2} = [0,1; 1,0];
model.link(P);

%% Generate synthetic dataset
n = 500;
ts = 1:n;
qlen_samples = 1.5 + rand(n,1)*0.5;  % aggregate queue-length observations

%% Create SampledMetric objects (aggregate, no class specified)
ql = SampledMetric(MetricType.QLen, ts, qlen_samples, node{2});

%% Estimate with MCMC
fprintf(1, '\n=== MCMC Estimator ===\n');
options = ParamEstimator.defaultOptions;
options.method = 'mcmc';
se = ParamEstimator(model, options);
se.addSamples(ql);
se.interpolate();
estVal = se.estimateAt(node{2});
fprintf(1, 'MCMC demands: Class1=%.4f, Class2=%.4f\n', estVal(1), estVal(2));

%% Solve model
solver{1} = MVA(model);
fprintf(1, '\nSOLVER: %s\n', solver{1}.getName());
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}
