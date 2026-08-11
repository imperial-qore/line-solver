clear node jobclass solver AvgTable
%% Example: EKF estimation with autoMethod selection
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% define model
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(model, 'Class1', 2, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{2}.setService(jobclass{1}, Exp(NaN));  % NaN = to be estimated

P = model.initRoutingMatrix;
P{1} = [0,1; 1,0];
model.link(P);

%% Generate synthetic dataset
n = 100;
ts = 1:n;
arvr_samples = 1.5*ones(n,1) - rand(n,1)*0.1;
util_samples = 0.4*arvr_samples;
respt_samples = 0.4./(1 - util_samples);

%% Create SampledMetric objects
lambda1 = SampledMetric(MetricType.ArvR, ts, arvr_samples, node{2}, jobclass{1});
respT1 = SampledMetric(MetricType.RespT, ts, respt_samples, node{2}, jobclass{1});
util = SampledMetric(MetricType.Util, ts, util_samples, node{2});

%% Estimate with EKF
fprintf(1, '\n=== EKF Estimator ===\n');
options = ParamEstimator.defaultOptions;
options.method = 'ekf';
se = ParamEstimator(model, options);
se.addSamples(lambda1);
se.addSamples(respT1);
se.addSamples(util);
se.interpolate();
estVal_ekf = se.estimateAt(node{2});
fprintf(1, 'EKF demand: Class1=%.4f (true=0.4000)\n', estVal_ekf(1));

%% Now demonstrate autoMethod
fprintf(1, '\n=== autoMethod selection ===\n');
node{2}.setService(jobclass{1}, Exp(NaN));
model.reset;

options2 = ParamEstimator.defaultOptions;
options2.method = 'auto';
se2 = ParamEstimator(model, options2);
se2.addSamples(lambda1);
se2.addSamples(respT1);
se2.addSamples(util);
se2.interpolate();
method = se2.autoMethod();
fprintf(1, 'autoMethod selected: %s\n', method);
fprintf(1, 'Required metrics for %s: %s\n', method, ParamEstimator.getRequiredMetrics(method));
estVal_auto = se2.estimateAt(node{2});
fprintf(1, 'Auto demand: Class1=%.4f (true=0.4000)\n', estVal_auto(1));

%% Solve model
solver{1} = MVA(model);
fprintf(1, '\nSOLVER: %s\n', solver{1}.getName());
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}
