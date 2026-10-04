function [V, info] = fluid_refine_meanfield(x, sigma2, Sigma, terms, covblk, epsrel)
% [V, INFO] = FLUID_REFINE_MEANFIELD(X, SIGMA2, SIGMA, TERMS, COVBLK, EPSREL)
%
% Refined mean field correction of a fluid fixed point (Gast, POMACS 2017).
%
% The mean-field fixed point x* is the leading term of an expansion of the
% true stationary mean in powers of the system size. The next term is
% obtained by carrying the second moment through the drift: writing A for the
% Jacobian at x* and B for its Hessian tensor, the correction V solves the
% linear system
%
%   A*V + (1/2) * sum_{j,k} Sigma_{jk} * d2F/dx_j dx_k = 0
%
% with Sigma the stationary covariance from FLUID_LYAPUNOV. Because Sigma
% scales with the population, V is the O(1/N) term of the expansion written
% directly in job counts, so no explicit density rescaling is needed. The
% Hessian contraction is evaluated without ever forming the tensor: writing
% Sigma = sum_m lam_m*v_m*v_m' by eigendecomposition,
%
%   sum_{jk} Sigma_{jk} d2F/dx_j dx_k = sum_m lam_m * d2F/dv_m^2
%
% and each directional second derivative is one central second difference, so
% the cost is O(rank(Sigma)) drift evaluations rather than O(n^2).
%
% The drift must be twice differentiable for this to mean anything. The
% first-order closure is only piecewise linear -- its second derivative is
% zero away from the kink and a delta at it -- so this function must be
% called on the Gaussian-closed drift, i.e. with the SIGMA2 that
% SOLVER_FLUID_MOMENTS converged to under options.method='refined'. Passing
% SIGMA2 = 0 is rejected rather than silently returning zero.
%
% Parameters:
%   x      - fluid fixed point
%   sigma2 - (M x 1) station population variances defining the smooth drift
%   Sigma  - (n x n) stationary covariance
%   terms  - representation from FLUID_MOMENT_TERMS
%   covblk - per-station covariance blocks closing the DPS share ratio
%   epsrel - relative step of the second difference (default: 1e-4)
%
% Returns:
%   V    - (n x 1) correction to be added to x
%   info - struct with fields rank, stepsize, residual
%
% See also SOLVER_FLUID_MOMENTS, FLUID_LYAPUNOV, FLUID_MIN_CLOSURE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5
    covblk = {};
end
if nargin < 6 || isempty(epsrel)
    epsrel = 1e-4;
end
if ~any(sigma2 > 0)
    line_error(mfilename,['The refined mean field expansion needs a twice-differentiable drift, but the ' ...
        'first-order closure is only piecewise linear. Reach this function through ' ...
        'options.method=''refined'', which converges the Gaussian closure first.']);
end

x = x(:);
n = numel(x);
F = @(z) terms.driftFcn(z, sigma2, covblk);

% eigendecomposition of the covariance, dropping numerically null directions
Sigma = (Sigma + Sigma')/2;
[Vec, Lam] = eig(Sigma);
lam = diag(Lam);
keep = lam > max(lam)*sqrt(eps) & lam > 0;
Vec = Vec(:,keep);
lam = lam(keep);

scale = max(1, norm(x));
step = epsrel*scale;
b = zeros(n,1);
F0 = F(x);
for m = 1:numel(lam)
    d = Vec(:,m);
    Fp = F(x + step*d);
    Fm = F(x - step*d);
    b = b + lam(m) * (Fp - 2*F0 + Fm) / step^2;
end
b = 0.5*b;

% solve A*V = -b on the reachable subspace, where A is invertible
A = terms.jacFcn(x, sigma2, covblk);
Vbasis = orth(terms.D);
Ar = Vbasis'*A*Vbasis;
condAr = cond(Ar);
if ~isfinite(condAr) || condAr > 1/sqrt(eps)
    line_error(mfilename,sprintf(['The fluid Jacobian is numerically singular on the reachable subspace ' ...
        '(condition number %.3g), so the refinement equation A*V = -b has no meaningful solution. The fixed ' ...
        'point sits at a drift kink or the model is marginally stable; use options.method=''minnormal'', which ' ...
        'resums the same correction without inverting A.'], condAr));
end
Vr = -Ar \ (Vbasis'*b);
V = Vbasis*Vr;

% the refinement is the next term of an asymptotic expansion, so it is only
% meaningful while it stays small against the leading term; a correction of
% the same size as the fixed point means the expansion has not kicked in at
% this population, and returning it would be worse than refusing
if norm(V) > 0.5*max(norm(x), sqrt(eps))
    line_error(mfilename,sprintf(['The 1/N refinement (norm %.3g) is not small against the mean-field fixed ' ...
        'point (norm %.3g), so the asymptotic expansion is outside its range of validity at this population. ' ...
        'Use options.method=''minnormal''.'], norm(V), norm(x)));
end

info = struct('rank',numel(lam),'stepsize',step,'residual',norm(A*V + b),'condition',condAr);
end
