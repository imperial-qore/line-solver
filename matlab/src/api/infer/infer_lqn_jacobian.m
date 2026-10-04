function [H, h0] = infer_lqn_jacobian(hfun, a, fdStep, fdFloor)
% INFER_LQN_JACOBIAN Forward finite-difference sensitivity matrix of h at a.
%
%   [H, H0] = INFER_LQN_JACOBIAN(HFUN, A, FDSTEP, FDFLOOR) returns the
%   sensitivity matrix H = dh/da evaluated at the parameter vector A by forward
%   finite differences, and H0 = HFUN(A). HFUN maps a parameter vector to a
%   numeric observation vector z = h(a). This is the approximate sensitivity
%   matrix H_k used in the EKF update (Zheng, Yang, Woodside, Litoiu, Iszlai,
%   "Tracking Time-Varying Parameters in Software Systems with Extended Kalman
%   Filters", CASCON 2005).
%
%   Each column is H(:,i) = (HFUN(a + d_i) - H0)/d_i with the per-parameter
%   step d_i = FDSTEP * max(|a_i|, FDFLOOR). HFUN is evaluated numel(A)+1 times.
%
%   FDSTEP  : relative perturbation (default 1e-3)
%   FDFLOOR : minimum absolute perturbation scale (default 1e-6)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(fdStep),  fdStep  = 1e-3; end
if nargin < 4 || isempty(fdFloor), fdFloor = 1e-6; end

a = a(:);
np = numel(a);
h0 = hfun(a);
h0 = h0(:);
no = numel(h0);

H = zeros(no, np);
for i = 1:np
    d = fdStep * max(abs(a(i)), fdFloor);
    ap = a;
    ap(i) = ap(i) + d;
    hi = hfun(ap);
    H(:, i) = (hi(:) - h0) / d;
end
end
