clear node jobclass solver AvgTable
%% Example: Compare trace-based estimators (MLPS, FMLPS, Gibbs) on a PS station
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% define model
N = 5; % population
model = Network('model');
node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(model, 'Class1', N, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{2}.setService(jobclass{1}, Exp(NaN));  % NaN = to be estimated

P = model.initRoutingMatrix;
P{1} = [0,1; 1,0];
model.link(P);

%% Generate synthetic trace data
n = 200;
arrival_times = sort(rand(n,1) * 100);
response_times = 0.5 + rand(n,1) * 0.3;
tput_samples = ones(n,1) * (N / (1.0 + 0.5)); % approximate throughput

%% Create trace-format SampledMetric objects
arvData = SampledMetric(MetricType.ArvR, arrival_times, arrival_times, node{2}, jobclass{1});
arvData.setTrace();
rtData = SampledMetric(MetricType.RespT, arrival_times, response_times, node{2}, jobclass{1});
rtData.setTrace();
tputData = SampledMetric(MetricType.Tput, arrival_times, tput_samples, node{2}, jobclass{1});

%% Estimate with MLPS
fprintf(1, '\n=== MLPS Estimator ===\n');
options = ParamEstimator.defaultOptions;
options.method = 'mlps';
se = ParamEstimator(model, options);
se.addSamples(arvData);
se.addSamples(rtData);
se.interpolate();
try
    estVal_mlps = se.estimateAt(node{2});
    fprintf(1, 'MLPS demand: Class1=%.4f\n', estVal_mlps(1));
catch e
    fprintf(1, 'MLPS skipped: %s\n', e.message);
    estVal_mlps = NaN;
end

%% Estimate with FMLPS
fprintf(1, '\n=== FMLPS Estimator ===\n');
node{2}.setService(jobclass{1}, Exp(NaN));
model.reset;

options.method = 'fmlps';
se = ParamEstimator(model, options);
se.addSamples(arvData);
se.addSamples(rtData);
se.interpolate();
try
    estVal_fmlps = se.estimateAt(node{2});
    fprintf(1, 'FMLPS demand: Class1=%.4f\n', estVal_fmlps(1));
catch e
    fprintf(1, 'FMLPS skipped: %s\n', e.message);
    estVal_fmlps = NaN;
end

%% Estimate with Gibbs
fprintf(1, '\n=== Gibbs Estimator ===\n');
node{2}.setService(jobclass{1}, Exp(NaN));
model.reset;

options.method = 'gibbs';
se = ParamEstimator(model, options);
se.addSamples(arvData);
se.addSamples(rtData);
se.addSamples(tputData);
se.interpolate();
try
    estVal_gibbs = se.estimateAt(node{2});
    fprintf(1, 'Gibbs demand: Class1=%.4f\n', estVal_gibbs(1));
catch e
    fprintf(1, 'Gibbs skipped: %s\n', e.message);
    estVal_gibbs = NaN;
end

%% Compare results
fprintf(1, '\n=== Comparison (true demand ~ 0.5) ===\n');
fprintf(1, 'MLPS:  %.4f\n', estVal_mlps(1));
fprintf(1, 'FMLPS: %.4f\n', estVal_fmlps(1));
fprintf(1, 'Gibbs: %.4f\n', estVal_gibbs(1));

%% Solve model with final estimates
solver{1} = MVA(model);
fprintf(1, '\nSOLVER: %s\n', solver{1}.getName());
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}
