clear node jobclass solver AvgTable

clearvars -except exampleName; 
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{3} = Source(model,'Source');
node{4} = Sink(model,'Sink');

jobclass{1} = OpenClass(model, 'Class1', 0);

node{1}.setService(jobclass{1}, HyperExp(0.5,3.0,10.0));
node{2}.setService(jobclass{1}, Exp(NaN)); % NaN = to be estimated
node{3}.setArrival(jobclass{1}, Exp(0.1));

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = cell(K,K);
P{1,1} = [0,1,0,0; 0,0,0,1; 1,0,0,0; 0,0,0,0];

model.link(P);

%% Generate random dataset for utilization and average arrival rate
n = 1000;
ts = 1:n; % timestamps at which the system is sampled
arvrate_samples = 2*ones(n,1)-rand(n,1)*0.15; % sampled arrival rate values
util_samples = ones(n,1)-rand(n,1)*0.05;  % sampled utilization values
respt_samples = 0.7./(1-util_samples);  % sampled response time values

%% Estimate demands
estoptions = ServiceEstimator.defaultOptions;
estoptions.method = 'ubr';
se = ServiceEstimator(model, estoptions);

lambda1 = SampledMetric(MetricType.ArvR, ts, arvrate_samples, node{2}, jobclass{1});
respT1 = SampledMetric(MetricType.RespT, ts, respt_samples, node{2}, jobclass{1});
util = SampledMetric(MetricType.Util, ts, util_samples, node{2});

se.addSamples(lambda1);
se.addSamples(respT1);
se.addSamples(util);
se.interpolate(); % ensures that the metrics are interpolated at the same time instants all
estVal = se.estimateAt(node{2}) % the change to the NaN value in Exp is applied automatically

%% Solve model
options = Solver.defaultOptions;
options.keep=true;
options.verbose=1;
options.cutoff = 10;
options.seed = 23000;
%options.samples=2e4;

solver = {};
solver{end+1} = SolverMVA(model);

AvgTable = cell(1,length(solver));
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());    
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
