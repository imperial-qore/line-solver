
    
%% Generate test dataset of queue length traces
for n=1:2
 %% define model
 model = quickModelRNN(false, {SchedStrategy.FCFS, SchedStrategy.FCFS, SchedStrategy.FCFS, SchedStrategy.FCFS, SchedStrategy.FCFS}, [[150, 30, 90, 45, 60]], [20, 40, 60, 20, 35], [100]);
 [timeI, QN] = generateQueueLengthTraces(model, 20);
 times{n} = timeI;
 traces{n} = QN;
end

node = model.getNodes();
jobclass = model.classes;
 
%% Estimate demands

options = ServiceEstimator.defaultOptions;
options.method = 'rnn';
se = ServiceEstimator(model, options);

%% Solve for test dataset of queue length traces
for n=1:length(times)
    QN = traces{n};
    timeI = times{n};
    aql1 = SampledMetric(MetricType.QLen, timeI, QN{1, 1}, node{1}, jobclass{1}); % transient queue-length
    aql2 = SampledMetric(MetricType.QLen, timeI, QN{2, 1}, node{2}, jobclass{1}); 
    aql3 = SampledMetric(MetricType.QLen, timeI, QN{3, 1}, node{3}, jobclass{1});
    aql4 = SampledMetric(MetricType.QLen, timeI, QN{4, 1}, node{4}, jobclass{1});
    aql5 = SampledMetric(MetricType.QLen, timeI, QN{5, 1}, node{5}, jobclass{1});

%    aql6 = SampledMetric(MetricType.QLen, timeI, QN{6, 1}, node{6}, jobclass{1});
%    aql7 = SampledMetric(MetricType.QLen, timeI, QN{7, 1}, node{7}, jobclass{1});
%    aql8 = SampledMetric(MetricType.QLen, timeI, QN{8, 1}, node{8}, jobclass{1});
%    aql9 = SampledMetric(MetricType.QLen, timeI, QN{9, 1}, node{9}, jobclass{1});
%    aql10 = SampledMetric(MetricType.QLen, timeI, QN{10, 1}, node{10}, jobclass{1});

    se.addSamples(aql1);
    se.addSamples(aql2);
    se.addSamples(aql3);
    se.addSamples(aql4);
    se.addSamples(aql5);

%    se.addSamples(aql6);
%    se.addSamples(aql7);
%    se.addSamples(aql8);
%    se.addSamples(aql9);
%    se.addSamples(aql10);
end


estVal = se.estimateAt(node)

%% Solve model
solver = {};
solver{end+1} = SolverMVA(model);

AvgTable = cell(1,length(solver));
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());    
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
