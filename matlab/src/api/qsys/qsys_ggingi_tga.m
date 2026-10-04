function result = qsys_ggingi_tga(lambda, mu, n, ca, cs, patienceCcdf, varargin)
% QSYS_GGINGI_TGA Truncated Gaussian approximation for the G/GI/n+GI queue.
%
% RESULT = QSYS_GGINGI_TGA(LAMBDA, MU, N, CA, CS, PATIENCECCDF) approximates the
% steady state of a heavily-loaded multiserver queue with abandonment: a general
% stationary arrival process of rate LAMBDA and variability CA, iid general
% service of mean 1/MU and variability CS, N servers, unlimited waiting room and
% iid general patience with complementary cdf PATIENCECCDF.
%
% THE APPROXIMATION IS A FLUID CENTRE PLUS A GAUSSIAN FLUCTUATION, TRUNCATED. In
% the efficiency-driven regime (rho > 1 held fixed as N grows) the fluid limit
% gives the centre -- all servers busy, waiting time w = F^-1(1-1/rho), queue
% Q = lambda int_0^w F^c -- and the many-server central limit theorem gives a
% NORMAL fluctuation of order sqrt(N) around it. Adding the two directly can
% produce negative queues and negative waits, so both are TRUNCATED at zero,
% which is what makes the formulas usable down to moderate overload; the paper
% reports good accuracy for rho > 1.02 and abandonment rates below 2.
%
% Three independent sources of variability enter separately, which is what lets
% the exponential-service formula be generalized: the service law appears only
% as the factor (cs+1)rho in sigma_W^2 (eq. 24), which is 2rho at cs = 1.
%
% An UNDERLOADED model (rho <= 1) has no queue in the limit; the number in system
% is then normal with the infinite-server variance, whose variance-to-mean ratio
% is the asymptotic peakedness of QSYS_GGNM_DIFFUSION.
%
% Options:
%   'patiencePdf', F  - the patience density; differenced from the ccdf when absent
%   'serviceCcdf', G  - G^c(x) = P(S > x), used only in the underloaded branch
%
% Returns a struct with fields regime, trafficIntensity, fluidWait,
% fluidQueueLength, meanWait, varWait, meanQueueLength, varQueueLength,
% meanNumberInService, meanNumber, probDelay, probAbandon, sigmaW and sigmaX.
%
% Example:
%   fc = @(x) exp(-0.5*x); fp = @(x) 0.5*exp(-0.5*x);
%   res = qsys_ggingi_tga(120, 1, 100, 1, 1, fc, 'patiencePdf', fp);
%
% Reference: Y. Liu, W. Whitt, Y. Yu (2016). Approximations for heavily-loaded
% G/GI/n+GI queues. Naval Research Logistics 63(3), 187-217.
%
% See also QSYS_GGISGI_FLUID, QSYS_ERLANGA, QSYS_GGNM_DIFFUSION.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

options = struct('patiencepdf', [], 'serviceccdf', []);
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
ca2 = ca^2;
rho = lambda/(n*mu);
lamPn = lambda/n;                       % the per-server rate of the reference's scaling

if isempty(options.patiencepdf)
    h = 1e-6;
    pdfFun = @(x) max(0, (patienceCcdf(max(0,x-h)) - patienceCcdf(x+h))/(2*h));
else
    pdfFun = options.patiencepdf;
end

result.trafficIntensity = rho;
if rho <= 1
    % Underloaded: no queue in the limit, and the number in system is normal
    % with the infinite-server variance (eq. 10).
    if isempty(options.serviceccdf)
        omega = 0.5;
    else
        hi = qsys_ggingi_tga_invccdf(options.serviceccdf, 1e-12);
        omega = qsys_ggingi_tga_simpson(@(x) options.serviceccdf(x)^2, 0, hi)*mu;
    end
    result.regime = 'underloaded';
    result.fluidWait = 0;
    result.fluidQueueLength = 0;
    result.meanWait = 0;
    result.varWait = 0;
    result.meanQueueLength = 0;
    result.varQueueLength = 0;
    result.meanNumberInService = lambda/mu;
    result.meanNumber = lambda/mu;
    result.probDelay = 0;
    result.probAbandon = 0;
    result.sigmaW = 0;
    result.sigmaX = sqrt((lambda/mu)*(1 + (ca2-1)*omega));
    return
end

% Overloaded: the fluid centre of Theorem 2.1(b).
w = qsys_ggingi_tga_invccdf(patienceCcdf, 1/rho);
fw = pdfFun(w);
if fw <= 0
    line_error(mfilename, ['the patience density vanishes at the fluid waiting time, so the ' ...
        'Gaussian correction is undefined there']);
end
qPerServer = lamPn * qsys_ggingi_tga_simpson(patienceCcdf, 0, w);

% Eq. (24): the service law enters only through the (cs+1)rho term, which is
% 2rho for exponential service and recovers eq. (11) there.
sigmaW2 = ((ca2 - 1) + (cs + 1)*rho)/(2*mu*rho^2*fw);
sigmaX2 = mu^2*sigmaW2 + lamPn*qsys_ggingi_tga_simpson( ...
    @(x) patienceCcdf(x)*(1 + (ca2-1)*patienceCcdf(x)), 0, w);
sigmaW = sqrt(max(sigmaW2, 0));
sigmaX = sqrt(max(sigmaX2, 0));

aW = sqrt(n)*w/sigmaW;                  % eq. (21)
aX = sqrt(n)*qPerServer/sigmaX;         % eq. (19)
[~, vW] = qsys_ggingi_tga_truncmoments(aW);
[~, vX] = qsys_ggingi_tga_truncmoments(aX);

result.regime = 'overloaded';
result.fluidWait = w;
result.fluidQueueLength = n*qPerServer;
result.meanWait = w*(qsys_ggingi_tga_Phi(aW) + qsys_ggingi_tga_phi(aW)/aW);
result.varWait = (sigmaW^2/n)*vW;
result.meanQueueLength = n*qPerServer*(qsys_ggingi_tga_Phi(aX) + qsys_ggingi_tga_phi(aX)/aX);
result.varQueueLength = n*sigmaX^2*vX;
% E[B] = E[min(X_n,n)]: every server is busy but for the lower tail of the
% Gaussian fluctuation.
result.meanNumberInService = n - sqrt(n)*sigmaX*(qsys_ggingi_tga_phi(aX) - ...
    aX*(1 - qsys_ggingi_tga_Phi(aX)));
result.meanNumber = result.meanNumberInService + result.meanQueueLength;
result.probDelay = qsys_ggingi_tga_Phi(aW);                       % eq. (22)
% Eq. (23): a customer abandons when its patience falls short of its wait.
pa = qsys_ggingi_tga_simpson(@(x) (1 - qsys_ggingi_tga_Phi(aW*(x/w - 1)))*pdfFun(x), 0, ...
    max(20*w, w + 20));
result.probAbandon = min(max(pa, 0), 1);
result.sigmaW = sigmaW;
result.sigmaX = sigmaX;
end

function y = qsys_ggingi_tga_phi(x)
% Standard normal density.
y = exp(-x^2/2)/sqrt(2*pi);
end

function y = qsys_ggingi_tga_Phi(x)
% Standard normal cdf, through erfc so no statistics toolbox is needed.
y = erfc(-x/sqrt(2))/2;
end

function [m1, v] = qsys_ggingi_tga_truncmoments(a)
% Mean and variance of max(Z,-a) for a standard normal Z, the shape every
% truncated Gaussian measure is built from.
m1 = qsys_ggingi_tga_phi(a) - a*(1 - qsys_ggingi_tga_Phi(a));
m2 = qsys_ggingi_tga_Phi(a) - a*qsys_ggingi_tga_phi(a) + a^2*(1 - qsys_ggingi_tga_Phi(a));
v = max(0, m2 - m1^2);
end

function v = qsys_ggingi_tga_simpson(f, a, b)
% Composite Simpson rule on a fixed even panel count.
if b <= a
    v = 0;
    return
end
m = 2000;
x = linspace(a, b, m+1);
y = arrayfun(f, x);
w = ones(1, m+1);
w(2:2:end-1) = 4;
w(3:2:end-2) = 2;
v = (b-a)/(3*m) * sum(w.*y);
end

function w = qsys_ggingi_tga_invccdf(ccdf, target)
% Smallest w with F^c(w) = target, by doubling then bisection.
lo = 0;
hi = 1;
while ccdf(hi) > target
    hi = 2*hi;
    if hi > 1e12
        line_error(mfilename, ['the patience ccdf never falls to 1/rho, so the overloaded model ' ...
            'has no fluid equilibrium']);
    end
end
while hi - lo > 1e-12*max(1, hi)
    mid = (lo + hi)/2;
    if ccdf(mid) > target
        lo = mid;
    else
        hi = mid;
    end
end
w = (lo + hi)/2;
end
