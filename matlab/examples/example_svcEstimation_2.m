clear node jobclass solver AvgTable

%% define model
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(model, 'Class1', 1, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 3, node{1}, 0);

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
arvr1_samples = ones(n,1)-rand(n,1)*0.15; 
arvr2_samples = 2*ones(n,1)-rand(n,1)*0.15; 
util_samples = 0.1 * arvr1_samples + 0.3 * arvr2_samples; 
respt1_samples = 0.1./(1-util_samples); 
respt2_samples = 0.3./(1-util_samples); 
aqlen1_samples = 1 + util_samples ./ (1-util_samples); 
aqlen2_samples = 1 + util_samples ./ (1-util_samples); 

%% Estimate demands
options = ServiceEstimator.defaultOptions;
options.method = 'erps';
se = ServiceEstimator(model, options);

aql1 = SampledMetric(MetricType.QLen, ts, aqlen1_samples, node{2}); % aggregate queue-length
aql1.setConditional(Event(EventType.ID_ARV, node{2}, jobclass{1})); % set that this metric is conditional on class-1 arrivals at node 2

aql2 = SampledMetric(MetricType.QLen, ts, aqlen2_samples, node{2});  % aggregate queue-length
aql2.setConditional(Event(EventType.ID_ARV, node{2}, jobclass{2})); % set that this metric is conditional on class-2 arrivals at node 2

respT1 = SampledMetric(MetricType.RespT, ts, respt1_samples, node{2}, jobclass{1});
respT2 = SampledMetric(MetricType.RespT, ts, respt2_samples, node{2}, jobclass{2});

se.addSamples(aql1);
se.addSamples(aql2);
se.addSamples(respT1);
se.addSamples(respT2);
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
