clear node jobclass solver AvgTable

%% define model
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(model, 'Class1', 2, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 3, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{1}.setService(jobclass{2}, Exp.fitMean(1.0));

node{2}.setService(jobclass{1}, Exp(NaN));  % NaN = to be estimated
node{2}.setService(jobclass{2}, Exp(NaN));  % NaN = to be estimated

P = model.initRoutingMatrix;
P{1} = [0,1; 1,0];
P{2} = [0,1; 1,0];
model.link(P);

%% Generate random dataset for queue-lengths
n = 1000;
ts = 1:n;
% synthetic queue-length data consistent with PS demands of 0.2 and 0.4
qlen1_samples = 0.3 * ones(n,1) + rand(n,1) * 0.1;
qlen2_samples = 0.8 * ones(n,1) + rand(n,1) * 0.1;

%% Estimate demands using QMLE
options = ParamEstimator.defaultOptions;
options.method = 'qmle';
se = ParamEstimator(model, options);

ql1 = SampledMetric(MetricType.QLen, ts, qlen1_samples, node{2}, jobclass{1});
ql2 = SampledMetric(MetricType.QLen, ts, qlen2_samples, node{2}, jobclass{2});

se.addSamples(ql1);
se.addSamples(ql2);
se.interpolate();
estVal = se.estimateAt(node{2})

%% Solve model
solver = {};
solver{end+1} = MVA(model);

AvgTable = cell(1,length(solver));
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
