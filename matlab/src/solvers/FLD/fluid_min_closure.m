function [h, dh, d2h] = fluid_min_closure(n, c, s2, vc, cov_nc)
% [H, DH, D2H] = FLUID_MIN_CLOSURE(N, C, S2, VC, COV_NC)
%
% Min-normal moment closure of E[min(X,Y)] for jointly normal X, Y.
%
% The mean-field ODEs of SolverFLD close the moment hierarchy at first order,
% replacing E[min(X,Y)] by min(E[X],E[Y]). That closure is exact only where
% min() is locally linear, so its error peaks at the kink, i.e. exactly where
% the two means meet -- rho ~ 1 at a queueing station, the "switch point" of
% the process-algebra literature. This function returns instead the
% expectation under a bivariate normal marginal, together with its derivative
% with respect to the first mean.
%
% This is the min-normal closure of Guenther, Stefanek and Bradley (EPEW/UKPEW
% 2012, LNCS 7587:32-47, eq. 4), implemented in their general two-population
% form:
%
%   E[min(X,Y)] = E[X]*Phi((E[Y]-E[X])/th) + E[Y]*Phi((E[X]-E[Y])/th)
%                 - th*phi((E[Y]-E[X])/th)
%   th = (Var[X] - 2*Cov[X,Y] + Var[Y])^(1/2)
%
% and the derivative with respect to E[X] is Phi((E[Y]-E[X])/th), the
% probability that X is the smaller of the two. SolverFLD only ever needs the
% specialisation Y = c, the deterministic server count: VC = 0 and COV_NC = 0
% give th = sqrt(S2) and the expression collapses to
% N - (N-C)*Phi(d) - sqrt(S2)*phi(d) with d = (N-C)/sqrt(S2), the
% truncated-normal form. The general arguments are kept so the function IS the
% published closure rather than one instance of it, and so a future
% state-dependent capacity (a population rather than a constant) needs no new
% derivation.
%
% With TH = 0 the expressions collapse to min(n,c) and to the indicator
% 1{n < c}, so the first-order closure is recovered exactly and callers share
% a single code path. The derivative at the kink is taken as 0, the right
% derivative of min(), which is the convention already implied by the strict
% inequality test in ODE_RATES_CLOSING_FACTORS.
%
% Parameters:
%   n      - mean of the first argument (scalar or vector)
%   c      - mean of the second argument (server count when deterministic)
%   s2     - variance of the first argument
%   vc     - variance of the second argument (default 0, deterministic)
%   cov_nc - covariance of the two arguments (default 0)
%
% Returns:
%   h   - E[min(X,Y)]
%   dh  - d/dE[X] E[min(X,Y)] = P(X < Y) under the normal marginal
%   d2h - d2/dE[X]^2 E[min(X,Y)] = -phi(d)/th, the density of the kink. min()
%         is piecewise linear, so its second derivative is carried entirely by
%         the atom at X = Y; smoothing over the normal marginal turns that atom
%         into the density. Zero on the degenerate branch, where the closure is
%         the first-order one and the kink is not smoothed at all.
%
% See also FLUID_CAPACITY_CLOSURE, SOLVER_FLUID_MOMENTS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(s2)
    s2 = 0;
end
if nargin < 4 || isempty(vc)
    vc = 0;
end
if nargin < 5 || isempty(cov_nc)
    cov_nc = 0;
end

n = n(:); c = c(:); s2 = s2(:); vc = vc(:); cov_nc = cov_nc(:);
sz = max([numel(n), numel(c), numel(s2), numel(vc), numel(cov_nc)]);
if isscalar(n), n = n*ones(sz,1); end
if isscalar(c), c = c*ones(sz,1); end
if isscalar(s2), s2 = s2*ones(sz,1); end
if isscalar(vc), vc = vc*ones(sz,1); end
if isscalar(cov_nc), cov_nc = cov_nc*ones(sz,1); end

th2 = s2 - 2*cov_nc + vc;
th2(th2 < 0) = 0; % a covariance beyond the Cauchy-Schwarz bound is not admissible

h = zeros(sz,1);
dh = zeros(sz,1);
d2h = zeros(sz,1);

% degenerate marginal: first-order closure
%
% THE INDICATOR CARRIES A BAND, and it is a cross-codebase requirement rather
% than a modelling choice. A saturated fluid fixed point sits exactly AT n = c,
% and each engine's ODE stops on its own residual: MATLAB lands at 1.0004 and
% the C++ port at 1 - 1.8e-13 on the same model, so a strict `n < c` reads
% saturated in one and unsaturated in the other. That flips a whole Jacobian row
% between zero and unit, and with it the hyperbolicity verdict that decides
% whether SolverFLD answers with 'minnormal' or falls back to the first-order
% method -- two visibly different answers from the same model. A population
% within FineTol of the server count IS at the kink, in every codebase.
deg = ~(th2 > 0) | isinf(c);
h(deg) = min(n(deg), c(deg));
dh(deg) = double((c(deg) - n(deg)) > GlobalConstants.FineTol*max(1, n(deg)));

sm = ~deg;
if any(sm)
    nsm = n(sm); csm = c(sm);
    th = sqrt(th2(sm));
    d = (nsm - csm) ./ th;
    % normcdf/normpdf without the Statistics Toolbox
    Phid = 0.5*erfc(-d/sqrt(2));
    phid = exp(-0.5*d.^2)/sqrt(2*pi);
    hsm = nsm.*(1-Phid) + csm.*Phid - th.*phid;
    dhsm = 1 - Phid;
    d2hsm = -phid ./ th;
    % A POPULATION IS NONNEGATIVE AND THE NORMAL MARGINAL IS NOT. For X >= 0
    % pathwise, min(X,c) >= 0, and min() being concave, Jensen puts
    % E[min(X,c)] <= min(E[X],c): the closure's value belongs to [0, min(n,c)].
    % The normal marginal has no such support, and the mass it places on X < 0
    % drags the expectation OUT of that range once the mean falls to about one
    % standard deviation -- n = 0, c = 1, th = 0.664 returns -0.019. A negative
    % expected number in service is a station that CREATES work, and every
    % trajectory starts with the queues empty, so every outer iterate past
    % sigma2 = 0 enters there: on CQN_Cox_CS_9 (Delay + PS + PS(c=5), N=6) the
    % first such window moved 8.7e3 of mass and the drift norm reached 5.4e9,
    % which with NonNegative on the integrator hid by clamping and never
    % returned. Project onto the admissible range, and take the derivatives of
    % the bound that binds so the Jacobian still matches h: at the ceiling
    % min(n,c) that is 1 below the server count and 0 above, at the floor 0.
    hi = min(nsm, csm);
    atLo = hsm < 0;
    atHi = hsm > hi;
    hsm = min(max(hsm, 0), hi);
    dhsm(atLo) = 0;
    dhsm(atHi) = double(nsm(atHi) < csm(atHi));
    d2hsm(atLo | atHi) = 0;
    h(sm) = hsm;
    dh(sm) = dhsm;
    d2h(sm) = d2hsm;
end
end
