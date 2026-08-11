clear node jobclass solver AvgTable
%% Example: Compare UBO and MLE estimators on a closed network
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% define model
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(model, 'Class1', 1, node{1}, 0);
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
n = 1000;
ts = 1:n;
arvr1_samples = 2*ones(n,1) - rand(n,1)*0.15;
arvr2_samples = ones(n,1) - rand(n,1)*0.15;
util_samples = 0.1*arvr1_samples + 0.3*arvr2_samples;
respt1_samples = 0.1./(1 - util_samples);
respt2_samples = 0.3./(1 - util_samples);

%% Create SampledMetric objects
lambda1 = SampledMetric(MetricType.ArvR, ts, arvr1_samples, node{2}, jobclass{1});
lambda2 = SampledMetric(MetricType.ArvR, ts, arvr2_samples, node{2}, jobclass{2});
respT1 = SampledMetric(MetricType.RespT, ts, respt1_samples, node{2}, jobclass{1});
respT2 = SampledMetric(MetricType.RespT, ts, respt2_samples, node{2}, jobclass{2});
util = SampledMetric(MetricType.Util, ts, util_samples, node{2});

%% Estimate with UBO
fprintf(1, '\n=== UBO Estimator ===\n');
options = ParamEstimator.defaultOptions;
options.method = 'ubo';
se = ParamEstimator(model, options);
se.addSamples(lambda1); se.addSamples(lambda2);
se.addSamples(respT1); se.addSamples(respT2);
se.addSamples(util);
se.interpolate();
estVal_ubo = se.estimateAt(node{2});
fprintf(1, 'UBO demands: Class1=%.4f, Class2=%.4f\n', estVal_ubo(1), estVal_ubo(2));

%% Estimate with MLE
fprintf(1, '\n=== MLE Estimator ===\n');
node{2}.setService(jobclass{1}, Exp(NaN));
node{2}.setService(jobclass{2}, Exp(NaN));
model.reset;

options.method = 'mle';
se = ParamEstimator(model, options);
se.addSamples(lambda1); se.addSamples(lambda2);
se.addSamples(respT1); se.addSamples(respT2);
se.addSamples(util);
se.interpolate();
estVal_mle = se.estimateAt(node{2});
fprintf(1, 'MLE demands: Class1=%.4f, Class2=%.4f\n', estVal_mle(1), estVal_mle(2));

%% Compare results
fprintf(1, '\n=== Comparison ===\n');
fprintf(1, 'True demands: Class1=0.1000, Class2=0.3000\n');
fprintf(1, 'UBO:  Class1=%.4f, Class2=%.4f\n', estVal_ubo(1), estVal_ubo(2));
fprintf(1, 'MLE:  Class1=%.4f, Class2=%.4f\n', estVal_mle(1), estVal_mle(2));

%% Solve model with final estimates
solver{1} = MVA(model);
fprintf(1, '\nSOLVER: %s\n', solver{1}.getName());
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}
