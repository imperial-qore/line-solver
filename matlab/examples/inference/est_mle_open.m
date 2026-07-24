clear node jobclass solver AvgTable
%% Example: MLE estimation on an open network
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% define model
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{3} = Source(model, 'Source');
node{4} = Sink(model, 'Sink');

jobclass{1} = OpenClass(model, 'Class1', 0);
jobclass{2} = OpenClass(model, 'Class2', 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{1}.setService(jobclass{2}, Exp.fitMean(1.0));
node{2}.setService(jobclass{1}, Exp(NaN)); % NaN = to be estimated (true demand = 0.1)
node{2}.setService(jobclass{2}, Exp(NaN)); % NaN = to be estimated (true demand = 0.3)
node{3}.setArrival(jobclass{1}, Exp(1.0)); % arrival rate = 1.0
node{3}.setArrival(jobclass{2}, Exp(0.5)); % arrival rate = 0.5

P = model.initRoutingMatrix;
P{1,1} = [0,1,0,0; 0,0,0,1; 1,0,0,0; 0,0,0,0];
P{2,2} = [0,1,0,0; 0,0,0,1; 1,0,0,0; 0,0,0,0];
model.link(P);

%% Generate synthetic dataset consistent with model
% True demands D1=0.1, D2=0.3, arrival rates lambda1=1.0, lambda2=0.5
% True utilization U = lambda1*D1 + lambda2*D2 = 0.1 + 0.15 = 0.25
% True response times R1 = D1/(1-U) = 0.133, R2 = D2/(1-U) = 0.4
n = 100;
ts = 1:n;
arvr1_samples = 1.0*ones(n,1) + randn(n,1)*0.05;
arvr2_samples = 0.5*ones(n,1) + randn(n,1)*0.03;
util_samples = 0.25*ones(n,1) + randn(n,1)*0.02;
util_samples = max(0.01, min(0.95, util_samples));
respt1_samples = 0.1./(1 - util_samples);
respt2_samples = 0.3./(1 - util_samples);

%% Create SampledMetric objects
lambda1 = SampledMetric(MetricType.ArvR, ts, arvr1_samples, node{2}, jobclass{1});
lambda2 = SampledMetric(MetricType.ArvR, ts, arvr2_samples, node{2}, jobclass{2});
respT1 = SampledMetric(MetricType.RespT, ts, respt1_samples, node{2}, jobclass{1});
respT2 = SampledMetric(MetricType.RespT, ts, respt2_samples, node{2}, jobclass{2});
util = SampledMetric(MetricType.Util, ts, util_samples, node{2});

%% Estimate with MLE
fprintf(1, '\n=== MLE Estimator (Open Network) ===\n');
options = ParamEstimator.defaultOptions;
options.method = 'mle';
se = ParamEstimator(model, options);
se.addSamples(lambda1);
se.addSamples(lambda2);
se.addSamples(respT1);
se.addSamples(respT2);
se.addSamples(util);
se.interpolate();
estVal = se.estimateAt(node{2});
fprintf(1, 'MLE demands: Class1=%.4f (true=0.1000), Class2=%.4f (true=0.3000)\n', estVal(1), estVal(2));

%% Solve model
solver{1} = MVA(model);
fprintf(1, '\nSOLVER: %s\n', solver{1}.getName());
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}
