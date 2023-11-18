%% define model

model = quickModel(false, {SchedStrategy.FCFS, SchedStrategy.FCFS}, [[1, 0.3]]);
node = model.getNodes();
jobclass = model.classes;

%% Generate test dataset

[avgResponseTime, avgArrivalRates, utilisation] = generateAvgSamples(model, 10, size(jobclass,1));

% Reset node service rates to be estimated
node{2}.setService(jobclass{1}, Exp(NaN));

%% Estimate demands

options = ServiceEstimator.defaultOptions;
options.method = 'ekf';
se = ServiceEstimator(model, options);


%% Solve for test dataset of queue length traces

ts = 1:size(utilisation,1);
lambda1 = SampledMetric(MetricType.ArvR, ts, avgArrivalRates(:, 2), node{2}, jobclass{1});
respT1 = SampledMetric(MetricType.RespT, ts, avgResponseTime(:, 2), node{2}, jobclass{1});
util = SampledMetric(MetricType.Util, ts, utilisation(:, 2), node{2});
se.addSamples(lambda1);
se.addSamples(respT1);
se.addSamples(util);

estVal = se.estimateAt(node{2})

%% Solve model
solver = {};
solver{end+1} = SolverMVA(model);

AvgTable = cell(1,length(solver));
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());    
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
