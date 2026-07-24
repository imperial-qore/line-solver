clear node jobclass solver AvgTable
%% Example: RNN-based estimation from queue-length traces
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% define model with multiple queues
model = infer_quick_model_rnn(false, ...
    {SchedStrategy.FCFS, SchedStrategy.FCFS, SchedStrategy.FCFS}, ...
    [[150, 30, 90]], [20, 40, 60], [100]);

node = model.getNodes();
jobclass = model.classes;

%% Generate queue-length traces
fprintf(1, '\n=== Generating queue-length traces ===\n');
for n = 1:2
    model_tmp = infer_quick_model_rnn(false, ...
        {SchedStrategy.FCFS, SchedStrategy.FCFS, SchedStrategy.FCFS}, ...
        [[150, 30, 90]], [20, 40, 60], [100]);
    [timeI, QN] = infer_generate_qlen_traces(model_tmp, 20);
    times{n} = timeI;
    traces{n} = QN;
end

% Reset service rates
for i = 1:length(node)
    if ~isa(node{i}, 'Source') && ~isa(node{i}, 'Sink')
        node{i}.setService(jobclass{1}, Exp(NaN));
    end
end

%% Estimate with RNN
fprintf(1, '\n=== RNN Estimator ===\n');
options = ParamEstimator.defaultOptions;
options.method = 'rnn';
se = ParamEstimator(model, options);

for n = 1:length(times)
    QN = traces{n};
    timeI = times{n};
    for i = 1:min(3, size(QN,1))
        ql = SampledMetric(MetricType.QLen, timeI, QN{i, 1}, node{i}, jobclass{1});
        se.addSamples(ql);
    end
end

estVal = se.estimateAt(node);
fprintf(1, 'RNN estimated demands:\n');
for i = 1:size(estVal, 1)
    fprintf(1, '  Node %d: %.4f\n', i, estVal(i, 1));
end

%% Solve model
solver{1} = MVA(model);
fprintf(1, '\nSOLVER: %s\n', solver{1}.getName());
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}
