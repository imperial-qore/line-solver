clear node jobclass solver AvgTable

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

%% Generate random dataset for utilization and average arrival rate
n = 1000;
ts = 1:n;
ts2 = 1:1:n;
arvr1_samples = 2*ones(n,1)-rand(n,1)*0.15; 
arvr2_samples = 3*ones(length(ts2),1)-rand(length(ts2),1)*0.25; 
util_samples = ones(n,1)-rand(n,1)*0.05; 
util1_samples = 0.4*(2*ones(n,1)-rand(n,1)*0.15); 

%% Estimate demands
options = ServiceEstimator.defaultOptions;
options.method = 'ubr';
se = ServiceEstimator(model, options);

lambda1 = SampledMetric(MetricType.ArvR, ts, arvr1_samples, node{2}, jobclass{1});
lambda2 = SampledMetric(MetricType.ArvR, ts2, arvr2_samples, node{2}, jobclass{2});
util1 = SampledMetric(MetricType.Util, ts, util1_samples, node{2}, jobclass{1});
util = SampledMetric(MetricType.Util, ts, util_samples, node{2});

se.addSamples(lambda1);
se.addSamples(lambda2);
se.addSamples(util);
se.addSamples(util1);
se.interpolate();
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
