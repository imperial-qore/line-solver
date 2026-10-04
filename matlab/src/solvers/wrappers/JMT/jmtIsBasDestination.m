function tf = jmtIsBasDestination(sn, ist, r)
% TF = JMTISBASDESTINATION(SN, IST, R)
%
% True when station IST is the RECEIVING side of a true-BAS relation for class
% R, i.e. an arrival of R that finds IST full must block an upstream station
% rather than be lost.
%
% LINE accepts the BAS declaration in two places -- on the blocking (upstream)
% station, as cqn_bas_blocking.m does, or on the full destination, as a model
% read back from JMT does -- and MNetwork/refreshLocalVars resolves both into
% sn.isbasdestination (BUG-83). Reading sn.droprule at the capped station sees
% only the second form, which is what made SolverJMT refuse the first one.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tf = false;
if isnan(ist) || ist < 1 || ~isfield(sn,'isbasdestination') || isempty(sn.isbasdestination)
    return
end
if size(sn.isbasdestination,1) < ist || size(sn.isbasdestination,2) < r
    return
end
tf = logical(sn.isbasdestination(ist, r));
end
