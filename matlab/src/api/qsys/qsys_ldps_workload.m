function [F, t, p] = qsys_ldps_workload(lambda, B, alpha, N, t, ngrid)
% [F,t,p] = QSYS_LDPS_WORKLOAD(lambda,B,alpha,N,t) - Stationary distribution of
% the quantity of work in a single-stage load-dependent processor sharing
% station with Poisson arrivals and blocking.
%
% Assumed model (Cohen, 1979, Sect. 9; the model of Sect. 7 with one stage):
%   - a single service stage fed by a Poisson arrival stream of rate lambda;
%   - blocking capacity N: a request arriving when N requests are already
%     present is lost and leaves no trace on the state;
%   - generalized processor sharing: when x requests are present each of them
%     accrues service at rate f(x), so the stage completes work at total rate
%     x*f(x). This function is parametrized by the LINE load-dependent total
%     rate scaling alpha(x)=x*f(x), i.e. the argument of setLoadDependence at
%     a PS station;
%   - required service times are i.i.d. with absolutely continuous
%     distribution B and finite mean beta.
% The result is insensitive to B beyond its shape only through the equilibrium
% residual distribution Psi below; the state probabilities p depend on B only
% through beta.
%
% Let psi_t denote the total amount of service still to be given to the
% requests present at time t. Cohen (1979) eqs. (9.1)-(9.3) give
%
%   Pr{psi_t < psi} = sum_{h=0}^{N} p_h Psi^{h*}(psi),
%   p_h = (rho^h/h!) phi(h) / sum_{k=0}^{N} (rho^k/k!) phi(k),   rho=lambda*beta
%   phi(h) = 1/prod_{k=1}^{h} f(k),   phi(0)=1,
%   Psi(psi) = int_0^psi (1-B(v))/beta dv,
%
% with Psi^{h*} the h-fold convolution of Psi and Psi^{0*} degenerate at zero.
% Substituting f(k)=alpha(k)/k the factorial cancels, leaving
%
%   p_h propto rho^h / prod_{k=1}^{h} alpha(k),
%
% which is the familiar load-dependent birth-death form. Psi is the
% equilibrium (residual life) distribution of B, so psi_t is a mixture of
% h-fold convolutions of residual service times, with an atom p_0 at zero.
%
% Inputs:
%   lambda : Poisson arrival rate (finite, positive)
%   B      : Distribution object for the required service time (continuous,
%            finite positive mean)
%   alpha  : rate scaling alpha(n)=n*f(n) for n=1..N (finite, positive)
%   N      : blocking capacity (finite positive integer)
%   t      : optional grid at which the CDF is returned. Default: an
%            automatically sized grid covering the bulk of the distribution.
%   ngrid  : optional number of points of the internal uniform quadrature
%            grid on which the convolutions are formed. Default 2001.
% Outputs:
%   F : F(j) = Pr{psi_t <= t(j)}. Note F(1)=p(1)=Pr{psi_t=0} when t(1)=0,
%       since the workload has an atom at zero of size p_0.
%   t : grid at which F is reported
%   p : (1,N+1) vector, p(h+1)=Pr{x_t=h}, the stationary number in system
%
% This is the model of Cohen (1979) Sect. 9 only. It is not the weighted
% GPS/DPS discipline of SchedStrategy.GPS, whose per-class weights this
% formula does not represent.
%
% Reference: J.W. Cohen, "The multiple phase service network with generalized
% processor sharing", Acta Informatica 12, 245-284 (1979), Sect. 9.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4
    line_error(mfilename,'Requires at least four inputs: lambda, B, alpha, N.');
end

%% gating: reject anything outside the assumed model
if ~isnumeric(lambda) || ~isscalar(lambda) || ~isfinite(lambda) || lambda <= 0
    line_error(mfilename,'lambda must be a finite positive scalar, the rate of the Poisson arrival stream.');
end
if ~isa(B,'Distribution')
    line_error(mfilename,'B must be a Distribution object giving the required service time.');
end
if B.isDisabled() || B.isImmediate()
    line_error(mfilename,'B must be an active service time distribution, not Disabled or Immediate.');
end
if ~B.isContinuous()
    line_error(mfilename,'B must be a continuous distribution: Cohen (1979) Sect. 9 assumes an absolutely continuous required service time.');
end
beta = B.getMean();
if ~isfinite(beta) || beta <= 0
    line_error(mfilename,'B must have a finite positive mean.');
end
if ~isnumeric(N) || ~isscalar(N) || ~isfinite(N) || N < 1 || N ~= round(N)
    line_error(mfilename,'N must be a finite positive integer, the blocking capacity of the service stage.');
end
if ~isnumeric(alpha) || ~isvector(alpha) || isempty(alpha)
    line_error(mfilename,'alpha must be a numeric vector holding the rate scaling alpha(n)=n*f(n).');
end
alpha = alpha(:).';
if numel(alpha) < N
    line_error(mfilename,sprintf('alpha must supply the rate scaling for n=1..%d, but only %d entries were given.',N,numel(alpha)));
end
alpha = alpha(1:N);
if any(~isfinite(alpha)) || any(alpha <= 0)
    line_error(mfilename,'alpha(n) must be finite and strictly positive for n=1..N, since every request in a busy stage is served at a positive rate.');
end

%% stationary number in system, eq. (9.1)
% p_h propto rho^h/prod_{k<=h} alpha(k), accumulated in logs so that large
% rho or large N do not overflow before normalization.
rho = lambda * beta;
logw = [0, cumsum(log(rho) - log(alpha))];
logw = logw - max(logw);
w = exp(logw);
p = w / sum(w);

%% grid
% Psi has mean m1e, the mean of the equilibrium residual life of B. The
% workload is a mixture of h-fold convolutions of Psi, so sizing the grid on
% the largest h carrying non-negligible mass covers the bulk of the support.
m1e = beta * (1 + B.getSCV()) / 2;
if ~isfinite(m1e) || m1e <= 0
    line_error(mfilename,'B must have a finite second moment: the equilibrium residual service time is otherwise undefined.');
end
userGrid = nargin >= 5 && ~isempty(t);
if userGrid
    if ~isnumeric(t) || ~isvector(t) || any(~isfinite(t)) || any(t < 0)
        line_error(mfilename,'t must be a vector of finite non-negative times.');
    end
    t = t(:).';
    tmax = max(t);
    if tmax <= 0
        tmax = m1e;
    end
else
    hmax = find(p > 1e-12, 1, 'last') - 1;
    if isempty(hmax) || hmax < 1
        hmax = 1;
    end
    tmax = m1e * (hmax + 8*sqrt(hmax));
    tmax = max(tmax, 8*m1e);
end

% see _kb/03-api-layer.md (qsys/ family) for rationale
if nargin < 6 || isempty(ngrid)
    ngrid = 2001;
end
if ~isnumeric(ngrid) || ~isscalar(ngrid) || ~isfinite(ngrid) || ngrid < 2 || ngrid ~= round(ngrid)
    line_error(mfilename,'ngrid must be a finite integer of at least 2.');
end
tg = linspace(0, tmax, ngrid);
dt = tg(2) - tg(1);

%% equilibrium residual service distribution, eq. (9.3)
% Density of Psi is (1-B(v))/beta. evalCDF is not vectorized by every
% Distribution subclass, so it is evaluated pointwise.
e = arrayfun(@(v) (1 - B.evalCDF(v)) / beta, tg);

%% workload distribution, eq. (9.2)
% F = sum_h p_h Psi^{h*}. The h=0 term is degenerate at zero, contributing the
% atom p_0 over the whole non-negative grid.
Fg = p(1) * ones(1, ngrid);
dens = e;
for h = 1:N
    if h > 1
        dens = local_convtrap(dens, e, dt, ngrid);
    end
    Fg = Fg + p(h+1) * cumtrapz(tg, dens);
end

if userGrid
    F = interp1(tg, Fg, t, 'linear');
else
    F = Fg;
    t = tg;
end
end

function c = local_convtrap(f, g, dt, n)
% Convolution of two densities sampled on a uniform grid, using the
% trapezoidal rule rather than the rectangle rule implied by a bare conv().
% The correction matters here because the equilibrium density does not vanish
% at the origin: e(0)=1/beta. Writing t_i=i*dt,
%
%   (f*g)(t_i) = int_0^{t_i} f(v) g(t_i-v) dv
%              ~ dt*[ sum_{j=0}^{i} f_j g_{i-j} - (f_0 g_i + f_i g_0)/2 ],
%
% i.e. conv() less half of each endpoint. Without the correction each
% convolution over-counts by dt*f_0*g_i, which accumulates over h and drives
% the mixture CDF above one.
c = conv(f, g);
c = c(1:n) * dt;
c = c - dt * (f(1)*g(1:n) + f(1:n)*g(1)) / 2;
end
