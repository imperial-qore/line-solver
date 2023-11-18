%% define model

model = quickModel(false, {SchedStrategy.INF, SchedStrategy.PS}, [[1,0.1];[1,0.3]]);
node = model.getNodes();
jobclass = model.classes;

%% Generate test dataset

N = 5;
[avgResponseTime, avgArrivalRates, utilisation] = generateAvgSamples(model, N, size(jobclass,1));

ts = 1:N;

%% Estimate demands
options = ServiceEstimator.defaultOptions;
options.method = 'ubo';
%options.variant = 'BUNDLE';
options.variant = 'NESTED';
se = ServiceEstimator(model, options);

lambda1 = SampledMetric(MetricType.ArvR, ts, avgArrivalRates(:, 2, 1), node{2}, jobclass{1});
lambda2 = SampledMetric(MetricType.ArvR, ts, avgArrivalRates(:, 2, 2), node{2}, jobclass{2});
respT1 = SampledMetric(MetricType.RespT, ts, avgResponseTime(:, 2, 1), node{2}, jobclass{1});
respT2 = SampledMetric(MetricType.RespT, ts, avgResponseTime(:, 2, 2), node{2}, jobclass{2});
util = SampledMetric(MetricType.Util, ts, utilisation(:, 2), node{2});

se.addSamples(lambda1);
se.addSamples(lambda2);
se.addSamples(respT1);
se.addSamples(respT2);
se.addSamples(util);
se.interpolate();
estVal = se.estimateAt({node{2}})

%% Solve model
solver = {};
solver{end+1} = SolverMVA(model);

AvgTable = cell(1,length(solver));
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());    
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
