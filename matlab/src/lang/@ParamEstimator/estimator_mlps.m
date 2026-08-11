function estVal = estimator_mlps(self, nodes)
node = nodes{1};
% ESTIMATORMLPS Maximum Likelihood for Processor Sharing
% Estimates service demands using the MLPS likelihood-based method.
% Supports closed, open, and mixed queueing networks.
% For open/mixed models, builds a closed equivalent using Z_r = N_r / lambda_r.
%
% Requires trace-format SampledMetric objects with per-request arrival
% timestamps and response times. PS stations only.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

sn = self.model.getStruct;
R = sn.nclasses;

if node.schedStrategy ~= SchedStrategy.PS
    error('The MLPS method is available only for processor sharing stations.');
end

% Build closed equivalent model if any class is open
hasOpen = any(sn.njobs == Inf);
if hasOpen
    [eqModel, eqNode] = self.buildClosedEquivalentForPS(node);
else
    eqModel = self.model;
    eqNode = node;
end

% Extract trace data from SampledMetrics
rt_all = [];
class_all = [];
at_all = [];
for r = 1:R
    jc = self.model.classes{r};

    % Arrival timestamps (trace format)
    arvData = self.getArvR(node, jc);
    if isempty(arvData)
        error('Arrival timestamp data for node %s in class %d is missing.', node.name, r);
    end
    if ~arvData.isTrace()
        error('MLPS estimator requires trace-format arrival data. Use setTrace() on the SampledMetric.');
    end

    % Response times (trace format)
    rtData = self.getRespT(node, jc);
    if isempty(rtData)
        error('Response time data for node %s in class %d is missing.', node.name, r);
    end
    if ~rtData.isTrace()
        error('MLPS estimator requires trace-format response time data. Use setTrace() on the SampledMetric.');
    end

    nSamples = length(rtData.data);
    rt_all = [rt_all; rtData.data];
    at_all = [at_all; arvData.data];
    class_all = [class_all; r * ones(nSamples, 1)];
end

% Compute queue lengths at arrival from trace data
n = length(at_all);
jobid = (1:n)';
ql = infer_compute_ql_at_arrival(at_all, jobid, rt_all, jobid, class_all, R);

% Sort by arrival time
[~, sortIdx] = sort(at_all);
rt_sorted = rt_all(sortIdx);
class_sorted = class_all(sortIdx);
ql = ql(sortIdx, :);

% Remove zero response times
valid = rt_sorted > 0;
rt_sorted = rt_sorted(valid);
class_sorted = class_sorted(valid);
ql = ql(valid, :);

estVal = infer_mlps(eqModel, eqNode, rt_sorted, class_sorted, ql);
estVal = estVal(:)';

end
