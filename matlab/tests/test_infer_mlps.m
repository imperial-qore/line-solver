clear node jobclass solver AvgTable
%% Test for MLPS estimator
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
rng(1);

%% define model with true demand (1 job per class, multiclass)
trueDemand = 0.5;
model = Network('model');
node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(model, 'Class1', 1, node{1}, 0);
node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{2}.setService(jobclass{1}, Exp.fitMean(trueDemand));
P = model.initRoutingMatrix;
P{1} = [0,1; 1,0];
model.link(P);

%% generate trace data from SSA simulation
solver_ssa = SolverSSA(model, 'seed', 1, 'samples', 10000);
samplePath = solver_ssa.sample(node{2}, 10000);

% extract per-class arrival and departure times at the queue
% with 1 job per class, each class-r ARV is followed by a class-r DEP of the same job
arv_times = [];
dep_times = [];
for e = 1:length(samplePath.event)
    ev = samplePath.event{e};
    if ev.node == node{2}.index && ev.class == 1
        if ev.event == EventType.ARV
            arv_times(end+1) = ev.t;
        elseif ev.event == EventType.DEP
            dep_times(end+1) = ev.t;
        end
    end
end
n = min(length(arv_times), length(dep_times));
arv_times = arv_times(1:n);
dep_times = dep_times(1:n);
response_times = (dep_times - arv_times)';
arrival_times = arv_times';

%% reset service for estimation
node{2}.setService(jobclass{1}, Exp(NaN));

%% create trace-format SampledMetric objects
arvData = SampledMetric(MetricType.ArvR, arrival_times, arrival_times, node{2}, jobclass{1});
arvData.setTrace();
rtData = SampledMetric(MetricType.RespT, arrival_times, response_times, node{2}, jobclass{1});
rtData.setTrace();

%% estimate demands
options = ParamEstimator.defaultOptions;
options.method = 'mlps';
se = ParamEstimator(model, options);
se.addSamples(arvData);
se.addSamples(rtData);
se.interpolate();
estVal = se.estimateAt(node{2})

assert(abs(estVal - trueDemand)/trueDemand < 0.10, ...
    sprintf('MLPS: estimated %.4f, relative error %.2f%% exceeds 10%%', estVal, 100*abs(estVal - trueDemand)/trueDemand));

%% solve model
solver{1} = SolverMVA(model);
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}
