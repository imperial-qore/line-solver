function estVal = estimator_qmle(self, nodes)
% ESTIMATORQMLE Quick Maximum Likelihood Estimation
% Estimates service demands from mean queue-lengths using the QMLE
% closed-form formula. Supports closed, open, and mixed queueing networks.
% For open classes, uses the open-to-closed equivalence Z_r = N_r / lambda_r.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

sn = self.model.getStruct;

if ~iscell(nodes)
    nodes = {nodes};
end

R = sn.nclasses;
M = length(nodes);

% Get effective population for open classes
if isfield(self.options, 'openPopulation') && ~isempty(self.options.openPopulation)
    Nopen = self.options.openPopulation;
else
    Nopen = 100;
end

% Determine population per class
N = zeros(1, R);
for r = 1:R
    if sn.njobs(r) < Inf
        N(r) = sn.njobs(r);
    else
        N(r) = Nopen;
    end
end

% Extract think times: from Delay nodes for closed classes,
% from Source arrival rates for open classes (Z_r = N_r / lambda_r)
Z = zeros(1, R);
allNodes = self.model.getNodes;
for n = 1:length(allNodes)
    if isa(allNodes{n}, 'Delay')
        svcProc = allNodes{n}.getService;
        for r = 1:R
            if sn.njobs(r) < Inf
                Z(r) = Z(r) + svcProc{r}.getMean();
            end
        end
    elseif isa(allNodes{n}, 'Source')
        svcProc = allNodes{n}.getService;
        for r = 1:R
            if sn.njobs(r) == Inf
                lambda_r = 1 / svcProc{r}.getMean();
                Z(r) = N(r) / lambda_r;
            end
        end
    end
end

% extract mean queue-lengths Q(M,R)
Q = zeros(M, R);
for n = 1:M
    for r = 1:R
        qlData = self.getQLen(nodes{n}, self.model.classes{r});
        if isempty(qlData)
            error('Queue-length data for node %s in class %d is missing.', nodes{n}.name, r);
        end
        if iscell(qlData)
            qlData = qlData{1};
        end
        Q(n, r) = mean(qlData.data);
    end
end

estVal = infer_qmle(Q, N, Z);

end
