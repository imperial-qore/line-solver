function tf = mam_srcproc_is_renewal(sn, jst)
% MAM_SRCPROC_IS_RENEWAL  True if station JST's class-1 process is renewal.
%
% Reads the (D0,D1) pair out of sn.proc and applies MAM_IS_RENEWAL_MAP. When no
% usable pair is present the answer is FALSE, which is the conservative side:
% the caller uses this to decide whether a marginal-only closed form may be
% applied, and applying one to a process whose correlation structure could not
% be established is exactly the failure this guard exists to prevent.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tf = false;
if jst < 1 || jst > numel(sn.proc)
    return;
end
procj = sn.proc{jst};
if isempty(procj) || numel(procj) < 1
    return;
end
pair = procj{1};
if ~iscell(pair) || numel(pair) < 2
    return;
end
D0 = pair{1};
D1 = pair{2};
if isempty(D0) || isempty(D1) || size(D0, 1) ~= size(D1, 1)
    return;
end
tf = mam_is_renewal_map(D0, D1);
end
