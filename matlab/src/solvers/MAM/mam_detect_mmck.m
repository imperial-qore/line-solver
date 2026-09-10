function [isMmck, muRate] = mam_detect_mmck(sn, ist, K, mmapNode)
% MAM_DETECT_MMCK Decide if station ist matches the M/M/c/K assumptions.
%
% Returns isMmck=true (and the shared service rate muRate) only when:
%   - the aggregated arrival MMAP at the node is single-phase (i.e. Poisson
%     superposition of class arrivals)
%   - every active class has Exp service (procid==EXP) at the station
%   - all active classes share the same exponential service rate
%
% INPUT
%   sn       - network struct
%   ist      - station index
%   K        - number of classes
%   mmapNode - cell array {D0, D1, D_c1, D_c2, ...} for the aggregated
%              arrival MMAP at the node
%
% OUTPUT
%   isMmck   - logical, true if the M/M/c/K closed form is exact
%   muRate   - the shared service rate (NaN when isMmck is false)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

isMmck = false;
muRate = NaN;

if isempty(mmapNode) || ~iscell(mmapNode) || size(mmapNode{1}, 1) ~= 1
    return;  % arrivals are not single-phase Poisson superposition
end

muVals = [];
for k=1:K
    if sn.procid(ist, k) ~= ProcessType.EXP
        % Disabled (NaN rate) classes are skipped; everything else must be Exp
        if isnan(sn.rates(ist, k))
            continue;
        end
        return;
    end
    if isnan(sn.rates(ist, k)) || sn.rates(ist, k) <= 0
        continue;  % no inflow for this class
    end
    muVals(end+1) = sn.rates(ist, k); %#ok<AGROW>
end

if isempty(muVals)
    return;
end

if max(muVals) - min(muVals) > 1e-9 * max(1, max(muVals))
    return;  % per-class service rates differ
end

isMmck = true;
muRate = muVals(1);
end
