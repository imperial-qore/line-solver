function tf = mam_is_markovian_map(D0, D1)
% TF = MAM_IS_MARKOVIAN_MAP(D0, D1)
%
% True when the pair (D0,D1) is a genuine Markovian arrival process: D0 has
% non-negative off-diagonal rates, D1 is non-negative, and (D0+D1) is an
% infinitesimal generator (zero row sums). A RAP or ME violates the sign
% conditions while still defining a valid point process, so a CTMC assembled
% from it is a rational generator whose stationary solution is a signed vector.
% Every consumer that builds a CTMC out of sn.proc must gate on this.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tf = false;
if isempty(D0) || isempty(D1)
    return;
end
ns = size(D0, 1);
if size(D0, 2) ~= ns || size(D1, 1) ~= ns || size(D1, 2) ~= ns
    return;
end
if any(~isfinite(D0(:))) || any(~isfinite(D1(:)))
    return;
end
tol = 1e-9 * max(1, max(abs([D0(:); D1(:)])));
offDiag = D0(~eye(ns) > 0);
tf = all(offDiag >= -tol) && all(D1(:) >= -tol) ...
    && all(abs(sum(D0 + D1, 2)) <= tol);
end
