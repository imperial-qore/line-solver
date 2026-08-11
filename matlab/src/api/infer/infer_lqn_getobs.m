function z = infer_lqn_getobs(names, metrics, obsSpec)
% INFER_LQN_GETOBS Extract an observation vector from solved LQN metrics.
%
%   Z = INFER_LQN_GETOBS(NAMES, METRICS, OBSSPEC) returns the column vector Z
%   of performance measures selected by OBSSPEC from the per-element numeric
%   metric vectors in METRICS. This is the measurement-model half of the LQN
%   parameter identification method (see INFER_LQN): it maps a solved model to
%   the observation vector z = h(a).
%
%   NAMES   : cell array of element names, as in LayeredNetworkStruct.names.
%   METRICS : struct with fields QLen, Util, RespT, Tput (each a vector indexed
%             by element, aligned with NAMES, as returned by getEnsembleAvg).
%   OBSSPEC : struct array; OBSSPEC(i) has fields:
%       .metric : 'RespT' | 'Util' | 'Tput' | 'QLen'
%       .name   : element name (processor for utilization; entry or reference
%                 task for response time; entry/task for throughput)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

no = numel(obsSpec);
z = zeros(no, 1);
for i = 1:no
    idx = find(strcmp(names, obsSpec(i).name), 1);
    if isempty(idx)
        line_error(mfilename, sprintf('Element ''%s'' not found in the LQN.', obsSpec(i).name));
    end
    switch lower(obsSpec(i).metric)
        case 'respt'
            z(i) = metrics.RespT(idx);
        case 'util'
            z(i) = metrics.Util(idx);
        case 'tput'
            z(i) = metrics.Tput(idx);
        case 'qlen'
            z(i) = metrics.QLen(idx);
        otherwise
            line_error(mfilename, sprintf('Unknown metric ''%s''.', obsSpec(i).metric));
    end
end
end
