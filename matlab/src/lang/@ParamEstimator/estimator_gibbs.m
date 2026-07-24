function estVal = estimator_gibbs(self, nodes)
node = nodes{1};
% ESTIMATORGIBBS Gibbs Sampling estimation from trace-level data
% Estimates service demands using Gibbs sampling on arrival timestamps,
% response times, and throughput data. Requires trace-format
% SampledMetric objects.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

sn = self.model.getStruct;
R = sn.nclasses;
nbCores = node.getNumberOfServers();

% build legacy data cell array for infer_gibbs()
% data format: data{3,k} = arrival timestamps (ms), data{4,k} = response
% times (s), data{6,k} = throughput
data = cell(6, R + 1);

for r = 1:R
    jc = self.model.classes{r};

    % arrival timestamps (trace format)
    arvData = self.getArvR(node, jc);
    if isempty(arvData)
        error('Arrival rate/timestamp data for node %s in class %d is missing.', node.name, r);
    end
    if arvData.isTrace()
        data{3, r} = arvData.data * 1000; % convert to ms
    else
        error('Gibbs estimator requires trace-format arrival data. Use setTrace() on the SampledMetric.');
    end

    % response times (trace format)
    rtData = self.getRespT(node, jc);
    if isempty(rtData)
        error('Response time data for node %s in class %d is missing.', node.name, r);
    end
    if rtData.isTrace()
        data{4, r} = rtData.data;
    else
        error('Gibbs estimator requires trace-format response time data. Use setTrace() on the SampledMetric.');
    end

    % throughput
    tputData = self.getTput(node, jc);
    if isempty(tputData)
        error('Throughput data for node %s in class %d is missing.', node.name, r);
    end
    data{6, r} = tputData.data;
end

tol = self.options.tol;
estVal = infer_gibbs(data, nbCores, tol);
estVal = estVal(:)';

end
