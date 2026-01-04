function nvec = ljd_delinearize(idx, cutoffs)
% NVEC = LJD_DELINEARIZE(IDX, CUTOFFS)
%
% Convert linearized index back to per-class population vector
%
% idx: 1-based linearized index
% cutoffs: [N1, N2, ..., NK] - per-class cutoffs
%
% Returns: [n1, n2, ..., nK] - per-class populations

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = length(cutoffs);
nvec = zeros(1, K);
idx = idx - 1; % convert to 0-based for computation

for k = 1:K
    base = cutoffs(k) + 1;
    nvec(k) = mod(idx, base);
    idx = floor(idx / base);
end
end
