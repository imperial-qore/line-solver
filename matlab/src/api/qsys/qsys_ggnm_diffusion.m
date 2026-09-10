function result = qsys_ggnm_diffusion(lambda, mu, n, m, ca, cs, varargin)
% QSYS_GGNM_DIFFUSION Diffusion approximation for the G/GI/n/m queue.
%
% RESULT = QSYS_GGNM_DIFFUSION(LAMBDA, MU, N, M, CA, CS) approximates the
% steady state of a queue with a general arrival process of rate LAMBDA and
% variability CA, iid general service of mean 1/MU and variability CS, N servers
% and M extra waiting spaces (M = Inf for an unbounded queue).
%
% THE APPROXIMATION IS ONE DIFFUSION WITH TWO REGIONS. Below the staffing level
% the queue behaves like an infinite-server system, whose limit is NORMAL with
% variance-to-mean ratio the ASYMPTOTIC PEAKEDNESS
%
%   z = 1 + (ca^2 - 1) omega_G,   omega_G = int G^c(x)^2 dx / int G^c(x) dx
%
% (eqs. 1.6-1.7); above it the queue behaves like a single-server queue, whose
% limit is EXPONENTIAL with variability v = (ca^2 + cs^2)/2 (eq. 3.7). The
% steady-state law is a normal piece spliced to an exponential piece and every
% measure below is an integral of that density (eq. 3.14).
%
% WHAT z SAYS. The service-time distribution enters the delay probability ONLY
% through omega_G, which is 1 for deterministic service, 1/2 for exponential,
% and falls toward 0 as service gets more variable. So at ca^2 = 1 the delay
% probability does not depend on the service law at all (z = 1), which is the
% long-standing M/GI/n-by-M/M/n approximation; away from ca^2 = 1 it does, and
% this is how much.
%
% With M = Inf the delay probability is alpha(beta/sqrt(z)) for the Halfin-Whitt
% function alpha (eq. 3.10), so this generalizes QSYS_MMK_QED.
%
% Options:
%   'serviceCcdf', G - G^c(x) = P(S > x); the exponential of rate MU by default
%   'tol', TOL       - service-tail cut for the peakedness integral, 1e-12
%   'panels', P      - Simpson panels for it, default 4000
%
% Returns a struct with fields beta (the QED server slack), gamma (the scaled
% waiting room), peakedness (z), peakednessWeight (omega_G), variability (v),
% probDelay, probBlock, meanQueueLength, meanNumber, meanWait, utilization,
% throughput and trafficIntensity.
%
% Example:
%   res = qsys_ggnm_diffusion(95, 1, 100, 20, 1.5, 0.8);
%
% Reference: W. Whitt (2004). A diffusion approximation for the G/GI/n/m queue.
% Operations Research 52(6), 922-941.
%
% See also QSYS_MMK_QED, QSYS_GIGK_APPROX_WHITT, QSYS_ERLANGA.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

options = struct('serviceccdf', [], 'tol', 1e-12, 'panels', 4000);
for i = 1:2:numel(varargin)
    if i+1 > numel(varargin)
        line_error(mfilename, sprintf('option %s has no value', char(varargin{i})));
    end
    name = lower(char(varargin{i}));
    if ~isfield(options, name)
        line_error(mfilename, sprintf('unknown option %s', char(varargin{i})));
    end
    options.(name) = varargin{i+1};
end

if lambda <= 0 || mu <= 0
    line_error(mfilename, 'The arrival and service rates must be positive.');
end
n = round(n);
if n < 1
    line_error(mfilename, 'The number of servers n must be at least 1.');
end
if m < 0
    line_error(mfilename, 'The number of extra waiting spaces m must be non-negative.');
end

ca2 = ca^2;
cs2 = cs^2;
ES = 1/mu;
rho = lambda/(n*mu);
beta = sqrt(n)*(1 - rho);                       % eq. (0.1)
if isfinite(m)
    gamma = m/sqrt(n);                          % eq. (0.3)
else
    gamma = Inf;
end

if isempty(options.serviceccdf)
    omega = 0.5;                                % exponential service
else
    omega = qsys_ggnm_diffusion_omega(options.serviceccdf, ES, options.tol, options.panels);
end
z = 1 + (ca2 - 1)*omega;                        % eq. (1.6), asymptotic peakedness
if z <= 0
    line_error(mfilename, ['the asymptotic peakedness came out non-positive; check ca and the ' ...
        'service ccdf']);
end
v = (ca2 + cs2)/2;                              % eq. (3.7) with the weight w = 1
b = beta/sqrt(z);
r = beta/v;                                     % rate of the exponential piece

% Mass on the exponential piece. The tail factor is 1-exp(-r*gamma), negative
% together with r when the queue is overloaded, so the ratio stays positive and
% alpha stays in (0,1) on both sides of beta = 0.
if isfinite(gamma)
    tail = -expm1(-r*gamma);
else
    tail = 1;
end
if abs(r) < 1e-14
    % beta = 0: the exponential piece degenerates to a uniform on [0,gamma].
    if ~isfinite(gamma)
        line_error(mfilename, 'with beta = 0 the queue needs a finite waiting room to be stable');
    end
    alpha = 1/(1 + qsys_ggnm_diffusion_Phi(b)/(qsys_ggnm_diffusion_phi(b)*gamma/sqrt(z)));
    meanAbove = gamma/2;
    densityAtTop = alpha/gamma;
else
    alpha = 1/(1 + b*qsys_ggnm_diffusion_Phi(b)/(qsys_ggnm_diffusion_phi(b)*tail));
    if isfinite(gamma)
        e = exp(-r*gamma);
        meanAbove = (1/r - (gamma + 1/r)*e)/tail;
        densityAtTop = alpha*r*e/tail;
    else
        meanAbove = 1/r;
        densityAtTop = 0;
    end
end

% Mean of the normal piece, N(-beta, z) conditioned below 0.
meanBelow = -beta - sqrt(z)*qsys_ggnm_diffusion_phi(b)/qsys_ggnm_diffusion_Phi(b);
meanScaled = (1-alpha)*meanBelow + alpha*meanAbove;

% Eq. (7.5): the loss rate of the diffusion at the upper boundary, divided by
% the arrival rate, is the density there times v over sqrt(n).
probBlock = 0;
if isfinite(gamma)
    probBlock = min(max(densityAtTop*v/sqrt(n), 0), 1);
end
meanQueue = sqrt(n)*alpha*meanAbove;
throughput = lambda*(1 - probBlock);

result.beta = beta;
result.gamma = gamma;
result.peakedness = z;
result.peakednessWeight = omega;
result.variability = v;
result.probDelay = alpha;
result.probBlock = probBlock;
result.meanQueueLength = meanQueue;
result.meanNumber = n + sqrt(n)*meanScaled;
if throughput > 0
    result.meanWait = meanQueue/throughput;
else
    result.meanWait = 0;
end
result.utilization = min(rho, 1);
result.throughput = throughput;
result.trafficIntensity = rho;
end

function y = qsys_ggnm_diffusion_phi(x)
% Standard normal density.
y = exp(-x.^2/2)/sqrt(2*pi);
end

function y = qsys_ggnm_diffusion_Phi(x)
% Standard normal cdf, through erfc so no statistics toolbox is needed.
y = erfc(-x/sqrt(2))/2;
end

function omega = qsys_ggnm_diffusion_omega(ccdf, ES, tol, panels)
% omega_G = int G^c(x)^2 dx / int G^c(x) dx of eq. (1.7), by Simpson on a grid
% cut where the ccdf is negligible. The denominator is E[S], so only the
% numerator is integrated.
hi = 1;
while ccdf(hi) > tol
    hi = 2*hi;
    if hi > 1e12
        line_error(mfilename, 'the service ccdf does not decay, so its peakedness is undefined');
    end
end
x = linspace(0, hi, panels+1);
y = arrayfun(@(xx) ccdf(xx)^2, x);
w = ones(1, panels+1);
w(2:2:end-1) = 4;
w(3:2:end-2) = 2;
omega = (hi/(3*panels)*sum(w.*y))/ES;
end
