% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
rng(1);

%% define model
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(model, 'Class1', 1, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 2, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{1}.setService(jobclass{2}, Exp.fitMean(1.0));

node{2}.setService(jobclass{1}, Exp(NaN));  % NaN = to be estimated
node{2}.setService(jobclass{2}, Exp(NaN));  % NaN = to be estimated

P = model.initRoutingMatrix;
P{1} = [0,1; 1,0];
P{2} = [0,1; 1,0];
model.link(P);

%% Generate random dataset for arrival rates and aggregate utilization
n = 1000;
ts = 1:n;
arvr1_samples = 1.5 + rand(n,1)*0.1;
arvr2_samples = 2.0 + rand(n,1)*0.1;
util_samples = 0.2*arvr1_samples + 0.4*arvr2_samples;

%% Estimate demands
options = ParamEstimator.defaultOptions;
options.method = 'ubr';
se = ParamEstimator(model, options);

lambda1 = SampledMetric(MetricType.ArvR, ts, arvr1_samples, node{2}, jobclass{1});
lambda2 = SampledMetric(MetricType.ArvR, ts, arvr2_samples, node{2}, jobclass{2});
util = SampledMetric(MetricType.Util, ts, util_samples, node{2});

se.addSamples(lambda1);
se.addSamples(lambda2);
se.addSamples(util);
se.interpolate();
estVal = se.estimateAt(node{2})

trueDemands = [0.2, 0.4];
assert(all(abs(estVal - trueDemands)./trueDemands < 0.10), ...
    sprintf('UBR: estimated [%.4f, %.4f] too far from true [%.4f, %.4f]', ...
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
