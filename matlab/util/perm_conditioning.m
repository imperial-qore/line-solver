function c = perm_conditioning(A)
% C = PERM_CONDITIONING(A)
%
% Log10 magnitude of the largest Ryser term for this orientation of A.
%
% The inclusion-exclusion expansion of the permanent is largest when every
% column is selected, giving prod_i (sum_j a_ij). Since perm(A) = perm(A'),
% the two orientations return the same value but not the same cancellation,
% so this is the quantity to minimize when choosing between them.
%
% Returns Inf when a row sum vanishes, so that such an orientation is never
% preferred.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

lineSums = sum(abs(A), 2);
if any(lineSums <= 0)
    c = Inf;
    return
end
c = sum(log10(lineSums));
end
