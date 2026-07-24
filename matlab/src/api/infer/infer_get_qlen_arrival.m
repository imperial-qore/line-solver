function ql = infer_get_qlen_arrival(data)
% INFER_GET_QLEN_ARRIVAL Compute queue lengths at arrival from cell data.
%
% Wrapper around infer_compute_ql_at_arrival for the legacy cell-based
% data format. Assumes data is available in standard format where
% data{3,k} contains arrival times (in ms) and data{4,k} contains
% response times for class k.
%
% Inputs:
%   data - cell array in standard format (6 x K+1)
%
% Returns:
%   ql   - 1 x K cell array, each cell is numSamples(k) x K matrix
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = size(data,2) - 1;

% Collect all samples across classes
at = [];
rt = [];
class = [];
numObs = zeros(1, K);
for k = 1:K
    numObs(k) = size(data{3,k}, 1);
    at = [at; data{3,k}/1000];  % convert ms to secs
    rt = [rt; data{4,k}];
    class = [class; k * ones(numObs(k), 1)];
end

% Compute queue lengths using shared utility
n = length(at);
jobid = (1:n)';
ql_unsorted = infer_compute_ql_at_arrival(at, jobid, rt, jobid, class, K);

% Split into per-class cell arrays in original order
ql = cell(1, K);
counter = 0;
for k = 1:K
    ql{k} = ql_unsorted(counter+1:counter+numObs(k), :);
    counter = counter + numObs(k);
end

end
