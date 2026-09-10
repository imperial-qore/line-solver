function [Sigma, info] = fluid_lyapunov(A, Qdiff, D, tol)
% [SIGMA, INFO] = FLUID_LYAPUNOV(A, QDIFF, D, TOL)
%
% Stationary covariance of the linear noise approximation.
%
% Around a fixed point x* of the fluid drift, the fluctuation process
% Z = (X - x*) obeys the linear stochastic differential equation
% dZ = A*Z*dt + sqrt(Qdiff)*dW, whose stationary covariance solves the
% Lyapunov equation
%
%   A*Sigma + Sigma*A' + Qdiff = 0,   Qdiff = D*diag(r(x*))*D'
%
% A is singular whenever the model conserves population: every closed class
% contributes a left null vector, so the equation has no unique solution on
% the full state space. It does have one on the reachable subspace, which is
% exactly range(D): the state can only move along jump directions, so the
% fluctuation lives there and nowhere else. Both A = D*diag(rateBase)*G and
% Qdiff map into range(D) as well, so restricting to an orthonormal basis V
% of range(D) is an exact reduction, not an approximation, and the reduced
% Lyapunov equation is nonsingular whenever the fixed point is stable.
%
% Parameters:
%   A     - (n x n) drift Jacobian at the fixed point
%   Qdiff - (n x n) diffusion matrix D*diag(r)*D'
%   D     - (n x nevents) jump matrix, spanning the reachable subspace
%   tol   - stability margin; eigenvalues of the reduced A with real part
%           above -tol are reported as non-hyperbolic (default: sqrt(eps))
%
% Returns:
%   Sigma - (n x n) stationary covariance, supported on range(D)
%   info  - struct with fields rank, maxRealEig, stable
%
% See also SOLVER_FLUID_MOMENTS, FLUID_DRIFT_JACOBIAN.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4 || isempty(tol)
    tol = sqrt(eps);
end

n = size(A,1);
V = orth(D); % orthonormal basis of the reachable subspace
if isempty(V)
    Sigma = zeros(n);
    info = struct('rank',0,'maxRealEig',-Inf,'stable',true);
    return
end

Ar = V'*A*V;
Qr = V'*Qdiff*V;
Qr = (Qr + Qr')/2;

ev = eig(Ar);
maxRe = max(real(ev));
info = struct('rank',size(V,2),'maxRealEig',maxRe,'stable',maxRe < -tol);
if ~info.stable
    % Identified so the caller can tell this apart from any other failure:
    % FLUID_RESOLVE_DEFAULT_METHOD cannot see a non-hyperbolic fixed point in
    % advance (it exists only once the mean is solved), so @SolverFLD/runAnalyzer
    % switches a RESOLVED 'minnormal' to the first-order method on this
    % identifier alone. An explicit options.method='minnormal' still fails.
    throw(MException('LINE:FluidNonHyperbolic', ...
        ['[%s.m] The fluid fixed point is not exponentially stable on the reachable ' ...
        'subspace (largest Jacobian eigenvalue has real part %g), so the linear noise approximation has ' ...
        'no stationary covariance. This happens at an unstable model or at a drift kink; use ' ...
        'options.method=''closing'' for the mean only.'], mfilename, maxRe));
end

% Bartels-Stewart via the base-MATLAB Sylvester solver: Ar*W + W*Ar' = -Qr
W = sylvester(Ar, Ar', -Qr);
W = (W + W')/2;

Sigma = V*W*V';
Sigma = (Sigma + Sigma')/2;
end
