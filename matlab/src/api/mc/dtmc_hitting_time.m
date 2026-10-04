function h = dtmc_hitting_time(P, targetStates)
% H = DTMC_HITTING_TIME(P, TARGETSTATES)
%
% Mean number of steps to reach any state in TARGETSTATES from each state of a
% DTMC with transition matrix P. Target states have zero hitting time; the
% others solve (I - P_NT) h_NT = 1 over the non-target block. Twin of the
% Python api.mc.dtmc_hitting_time; TARGETSTATES is 1-based here and 0-based
% there, as elsewhere between the two codebases.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = size(P,1);
targetStates = unique(reshape(targetStates, 1, []));
isTarget = false(1,n);
isTarget(targetStates) = true;
nonTarget = find(~isTarget);

h = zeros(n,1);
if isempty(nonTarget)
    return
end

A = eye(length(nonTarget)) - P(nonTarget, nonTarget);
b = ones(length(nonTarget),1);
warnState = warning('off','MATLAB:singularMatrix');
hNT = A \ b;
warning(warnState);
if any(~isfinite(hNT))
    % A state that cannot reach the target set has infinite hitting time; the
    % least-squares solution of the singular system is not that answer.
    hNT = lsqminnorm(A, b);
    hNT(~isfinite(hNT)) = Inf;
end
h(nonTarget) = hNT;
end
