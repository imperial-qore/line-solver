function nvec = ljd_delinearize(idx, cutoffs)
% NVEC = LJD_DELINEARIZE(IDX, CUTOFFS)
%
% Inverse of LJD_LINEARIZE: recover the per-class population vector from its
% linearized index.
%
% idx: 1-based linearized index
% cutoffs: [N1, N2, ..., NK] - per-class cutoffs
%
% Returns: [n1, n2, ..., nK]
%
% The forward map is idx = 1 + n1 + n2*(N1+1) + n3*(N1+1)*(N2+1) + ..., i.e. a
% mixed-radix numeral with class k in radix (Nk+1), so the inverse is the
% digit-by-digit division that reads that numeral back.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = length(cutoffs);
nvec = zeros(1, K);
rem = idx - 1;
for k = 1:K
    radix = cutoffs(k) + 1;
    nvec(k) = mod(rem, radix);
    rem = floor(rem / radix);
end
end
