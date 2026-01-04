function idx = ljd_linearize(nvec, cutoffs)
% IDX = LJD_LINEARIZE(NVEC, CUTOFFS)
%
% Convert per-class population vector to linearized index
%
% nvec: [n1, n2, ..., nK] - per-class populations
% cutoffs: [N1, N2, ..., NK] - per-class cutoffs
%
% Returns: 1-based linearized index
%
% Index formula: idx = 1 + n1 + n2*(N1+1) + n3*(N1+1)*(N2+1) + ...

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = length(nvec);
idx = 1; % 1-indexed for MATLAB
multiplier = 1;

for k = 1:K
    nk = min(nvec(k), cutoffs(k)); % clamp to cutoff
    idx = idx + nk * multiplier;
    multiplier = multiplier * (cutoffs(k) + 1);
end
end
