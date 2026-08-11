clear node jobclass solver AvgTable
%% Example: ERPS estimation on an open network
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% define model
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
node{3} = Source(model, 'Source');
node{4} = Sink(model, 'Sink');

jobclass{1} = OpenClass(model, 'Class1', 0);
jobclass{2} = OpenClass(model, 'Class2', 0);

node{1}.setService(jobclass{1}, HyperExp(0.5, 3.0, 10.0));
node{1}.setService(jobclass{2}, HyperExp(0.5, 2.0, 8.0));
node{2}.setService(jobclass{1}, Exp(NaN)); % NaN = to be estimated
node{2}.setService(jobclass{2}, Exp(NaN)); % NaN = to be estimated
node{3}.setArrival(jobclass{1}, Exp(0.1));
node{3}.setArrival(jobclass{2}, Exp(0.05));

P = model.initRoutingMatrix;
P{1,1} = [0,1,0,0; 0,0,0,1; 1,0,0,0; 0,0,0,0];
P{2,2} = [0,1,0,0; 0,0,0,1; 1,0,0,0; 0,0,0,0];
model.link(P);

%% Generate synthetic dataset
n = 1000;
ts = 1:n;
arvr1_samples = ones(n,1) - rand(n,1)*0.15;
arvr2_samples = 2*ones(n,1) - rand(n,1)*0.15;
util_samples = 0.1*arvr1_samples + 0.3*arvr2_samples;
respt1_samples = 0.1./(1 - util_samples);
respt2_samples = 0.3./(1 - util_samples);
aqlen1_samples = 1 + util_samples./(1 - util_samples);
aqlen2_samples = 1 + util_samples./(1 - util_samples);

%% Create SampledMetric objects
options = ParamEstimator.defaultOptions;
options.method = 'erps';
se = ParamEstimator(model, options);

aql1 = SampledMetric(MetricType.QLen, ts, aqlen1_samples, node{2});
aql1.setConditional(Event(EventType.ARV, node{2}, jobclass{1}));

aql2 = SampledMetric(MetricType.QLen, ts, aqlen2_samples, node{2});
aql2.setConditional(Event(EventType.ARV, node{2}, jobclass{2}));

respT1 = SampledMetric(MetricType.RespT, ts, respt1_samples, node{2}, jobclass{1});
respT2 = SampledMetric(MetricType.RespT, ts, respt2_samples, node{2}, jobclass{2});

se.addSamples(aql1);
se.addSamples(aql2);
se.addSamples(respT1);
se.addSamples(respT2);
se.interpolate();
estVal = se.estimateAt(node{2})

%% Solve model
solver{1} = MVA(model);
fprintf(1, '\nSOLVER: %s\n', solver{1}.getName());
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}
