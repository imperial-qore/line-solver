function result = qsys_ggisgi_fluid(lambda, mu, s, patienceCcdf, varargin)
% QSYS_GGISGI_FLUID Steady state of the G/GI/s+GI fluid model.
%
% RESULT = QSYS_GGISGI_FLUID(LAMBDA, MU, S, PATIENCECCDF) computes the unique
% steady state of the deterministic fluid model of a multiserver queue with
% customer abandonment: arrival rate LAMBDA, S servers each of rate MU, general
% service times and general patience times whose complementary cdf is the
% function handle PATIENCECCDF, F^c(t) = P(patience > t).
%
% THE MODEL. Scale the content by S and let S grow. Customers become quanta of
% fluid but their sojourns do not shrink, so the ages survive the limit: the
% state is the density b(x) of fluid that has been IN SERVICE for time x and the
% density q(x) of fluid that has been WAITING for time x. With rho = lambda/(s*mu),
%
%   rho <= 1  b(x) = rho G^c(x),  q = 0,  no abandonment, no wait
%   rho >  1  b(x) = G^c(x),      q(x) = rho F^c(x) on [0,w] and 0 beyond
%
% where the queue boundary w solves F^c(w) = 1/rho (eq. 3.6). That one equation
% carries the whole overloaded regime: fluid that survives its patience for w
% enters service, so the surviving fraction F^c(w) must equal the fraction
% 1/rho the servers can absorb.
%
% WHAT THE DISTRIBUTIONS CONTRIBUTE (Corollary 3.1). The rates and the number in
% service depend on G and F only through their means. The wait w, the queue
% content and its age profile depend on F BEYOND its mean but on G only through
% its mean. The age profile in service depends on G beyond its mean. Neither the
% number of servers nor anything about the arrival process beyond its rate
% appears at all, which is why this model says nothing about the QED regime and
% everything about the overloaded one.
%
% Options:
%   'servingCcdf', G   - the service-time ccdf, needed only for the in-service
%                        age density; defaults to exponential of rate MU
%   'agePoints', X     - ages at which to return the two densities
%   'tol', T           - bisection tolerance for w, default 1e-12
%   'maxTime', TMAX    - largest age searched for w, default grows automatically
%
% Returns a struct with fields:
%   regime           - 'underloaded', 'balanced' or 'overloaded'
%   trafficIntensity - rho = lambda/(s*mu)
%   offeredWait      - w, the wait of every customer who is served, 0 unless overloaded
%   meanWait         - E[W] over all customers, = int_0^w F^c(t)dt = m_a F_e(w)
%   meanWaitServed   - w again, the fluid wait being deterministic
%   meanWaitAbandon  - E[patience | patience <= w]
%   probAbandon      - 1 - 1/rho when overloaded, 0 otherwise
%   meanQueueLength  - Q = lambda * meanWait, in customers
%   meanNumberInService - B = min(lambda/mu, s), in customers
%   meanNumber       - B + Q
%   utilization      - min(rho,1)
%   throughput       - min(lambda, s*mu)
%   abandonRate      - lambda - throughput
%   agePoints, serviceAgeDensity, queueAgeDensity - the densities per server at
%                      the requested ages, b(x) and q(x)
%
% ACCURACY. This is the s -> Inf limit, so it is an approximation at finite S
% that improves with S and with the overload. Whitt (2004) reports it as crude
% at s = 100, rho = 1.02 and good at s = 100, rho = 1.10.
%
% Example:
%   % 100 agents, 10% overload, exponential patience of mean 5
%   res = qsys_ggisgi_fluid(110, 1, 100, @(t) exp(-t/5));
%
% Reference: W. Whitt (2006). Fluid models for multiserver queues with
% abandonments. Operations Research 54(1), 37-54, Theorem 3.1 and Corollary 3.2.
%
% See also QSYS_ERLANGA, QSYS_MGISRGI_WHITT.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

options = struct('servingccdf', [], 'agepoints', [], 'tol', 1e-12, 'maxtime', []);
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

if lambda <= 0
    line_error(mfilename, 'The arrival rate lambda must be positive.');
end
if mu <= 0
    line_error(mfilename, 'The service rate mu must be positive.');
end
if s < 1
    line_error(mfilename, 'The number of servers s must be at least 1.');
end
if ~isa(patienceCcdf, 'function_handle')
    line_error(mfilename, 'The patience ccdf must be a function handle F^c(t) = P(T > t).');
end
servingCcdf = options.servingccdf;
if isempty(servingCcdf)
    servingCcdf = @(x) exp(-mu*x);
end

rho = lambda / (s*mu);
result.trafficIntensity = rho;

if rho <= 1
    % Underloaded and balanced, eq. (3.2): the queue is empty and the model is
    % the infinite-server fluid model, so nothing waits and nothing abandons.
    if abs(rho - 1) <= eps(1)
        result.regime = 'balanced';
    else
        result.regime = 'underloaded';
    end
    w = 0;
    meanWait = 0;
    probAbandon = 0;
    meanWaitAbandon = 0;
else
    result.regime = 'overloaded';
    % Eq. (3.6): F^c(w) = 1/rho. F^c is non-increasing, so bisect on the first
    % bracket found by doubling.
    target = 1/rho;
    w = qsys_ggisgi_fluid_invccdf(patienceCcdf, target, options.tol, options.maxtime);
    % Eq. (3.14): W = int_0^w F^c(t) dt = m_a F_e(w), the mean over ALL fluid,
    % served and abandoning alike.
    meanWait = qsys_ggisgi_fluid_integral(patienceCcdf, 0, w);
    probAbandon = 1 - 1/rho;
    % E[T | T <= w] = (W - w F^c(w)) / F(w) by parts, F^c(w) = 1/rho.
    meanWaitAbandon = (meanWait - w/rho) / probAbandon;
end

result.offeredWait = w;
result.meanWait = meanWait;
result.meanWaitServed = w;
result.meanWaitAbandon = meanWaitAbandon;
result.probAbandon = probAbandon;
result.meanQueueLength = lambda * meanWait;             % eq. (3.11), Little's law
result.meanNumberInService = min(lambda/mu, s);
result.meanNumber = result.meanNumberInService + result.meanQueueLength;
result.utilization = min(rho, 1);
result.throughput = min(lambda, s*mu);
result.abandonRate = lambda - result.throughput;

if ~isempty(options.agepoints)
    x = options.agepoints(:).';
    sigma = min(rho, 1);                                % fluid rate into service, per server
    result.agePoints = x;
    result.serviceAgeDensity = sigma * arrayfun(servingCcdf, x);
    if rho > 1
        result.queueAgeDensity = rho * arrayfun(patienceCcdf, x) .* (x <= w);
    else
        result.queueAgeDensity = zeros(size(x));
    end
end
end

function w = qsys_ggisgi_fluid_invccdf(ccdf, target, tol, maxTime)
% Smallest w with F^c(w) = target, found by doubling then bisection. F^c is
% non-increasing, so the doubling either brackets the crossing or proves that
% the patience law never decays that far.
lo = 0;
if ccdf(0) < target
    line_error(mfilename, 'the patience ccdf is below 1/rho at t = 0, so it is not a ccdf');
end
if isempty(maxTime)
    hi = 1;
    while ccdf(hi) > target
        hi = 2*hi;
        if hi > 1e12
            line_error(mfilename, ['the patience ccdf never falls to 1/rho, so the overloaded ' ...
                'fluid model has no equilibrium: too little of the fluid is willing to abandon']);
        end
    end
else
    hi = maxTime;
    if ccdf(hi) > target
        line_error(mfilename, 'the patience ccdf is still above 1/rho at maxTime');
    end
end
while hi - lo > tol*max(1, hi)
    mid = (lo + hi)/2;
    if ccdf(mid) > target
        lo = mid;
    else
        hi = mid;
    end
end
w = (lo + hi)/2;
end

function v = qsys_ggisgi_fluid_integral(f, a, b)
% Composite Simpson rule on a fixed fine grid: the integrand is a ccdf, hence
% monotone and bounded, so a fixed grid is enough and is reproducible.
if b <= a
    v = 0;
    return
end
n = 2000;
x = linspace(a, b, n+1);
y = arrayfun(f, x);
wgt = ones(1, n+1);
wgt(2:2:end-1) = 4;
wgt(3:2:end-2) = 2;
v = (b - a)/(3*n) * sum(wgt .* y);
end
