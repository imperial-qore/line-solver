clear node jobclass solver AvgTable
%% Example: EKF estimation on an open network
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% define model
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{3} = Source(model, 'Source');
node{4} = Sink(model, 'Sink');

jobclass{1} = OpenClass(model, 'Class1', 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{2}.setService(jobclass{1}, Exp(NaN)); % NaN = to be estimated (true demand = 0.4)
node{3}.setArrival(jobclass{1}, Exp(1.0)); % arrival rate = 1.0

P = model.initRoutingMatrix;
P{1,1} = [0,1,0,0; 0,0,0,1; 1,0,0,0; 0,0,0,0];
model.link(P);

%% Generate synthetic dataset consistent with model
% True demand D=0.4, arrival rate lambda=1.0
% True utilization U = lambda*D = 0.4
% True response time R = D/(1-U) = 0.4/0.6 = 0.667
n = 100;
ts = 1:n;
arvr_samples = 1.0*ones(n,1) + randn(n,1)*0.05;
util_samples = 0.4*ones(n,1) + randn(n,1)*0.02;
util_samples = max(0.01, min(0.95, util_samples));
respt_samples = 0.4./(1 - util_samples);

%% Create SampledMetric objects
lambda1 = SampledMetric(MetricType.ArvR, ts, arvr_samples, node{2}, jobclass{1});
respT1 = SampledMetric(MetricType.RespT, ts, respt_samples, node{2}, jobclass{1});
util = SampledMetric(MetricType.Util, ts, util_samples, node{2});

%% Estimate with EKF
fprintf(1, '\n=== EKF Estimator (Open Network) ===\n');
options = ParamEstimator.defaultOptions;
options.method = 'ekf';
se = ParamEstimator(model, options);
se.addSamples(lambda1);
se.addSamples(respT1);
se.addSamples(util);
se.interpolate();
estVal = se.estimateAt(node{2});
fprintf(1, 'EKF demand: Class1=%.4f (true=0.4000)\n', estVal(1));

%% Solve model
solver{1} = MVA(model);
fprintf(1, '\nSOLVER: %s\n', solver{1}.getName());
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}
