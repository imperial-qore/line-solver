function result = qsys_mtginf(lambdaFun, serviceCcdf, ES, tvals, varargin)
% QSYS_MTGINF Exact time-varying analysis of the Mt/G/Inf queue.
%
% RESULT = QSYS_MTGINF(LAMBDAFUN, SERVICECCDF, ES, TVALS) evaluates the
% infinite-server queue with a non-homogeneous Poisson arrival rate LAMBDAFUN
% and iid service times with complementary cdf SERVICECCDF and mean ES, at the
% times TVALS.
%
% THE RESULT IS EXACT, not an approximation. With infinitely many servers
% customers never interact, so the model is a Poisson random measure and the
% number in system at time t is POISSON with mean
%
%   m(t) = E[ int_{t-S}^{t} lambda(u) du ] = ES * E[lambda(t - Se)]
%        = int_0^Inf lambda(t-x) P(S > x) dx
%
% where Se is the STATIONARY-EXCESS (equilibrium) law of the service time, with
% density P(S>x)/ES. Because the law is Poisson, the variance equals the mean
% and every quantile follows from it.
%
% THE PHYSICS. Reading m(t) as ES*E[lambda(t-Se)] says the time-varying load is
% the stationary load ES*lambda(t) subjected to a TIME LAG and a SPACE SHIFT: to
% first order m(t) ~ ES*lambda(t - E[Se]) with E[Se] = E[S^2]/(2*ES), so peak
% congestion LAGS peak arrival rate, and by more than the mean service time when
% the service law is variable. The pointwise stationary approximation
% ES*lambda(t) is the zeroth-order term of the same expansion, which is exactly
% why it misses the lag.
%
% Options:
%   'startTime', T0  - the system started empty at T0; the default -Inf assumes
%                      the arrival rate has been running forever
%   'ES2', M2        - the second moment of the service time, which adds the lag
%                      E[Se] and the first-order lag approximation to the output
%   'servicePdf', G  - the service density, used for the exact departure rate;
%                      without it the departure rate comes from the flow balance
%                      m'(t) = lambda(t) - delta(t) by a central difference
%   'tol', TOL       - service-tail cut for the age integral, default 1e-12
%   'panels', N      - Simpson panels for that integral, default 4000
%   'maxAge', A      - cap on the age integrated over, default 1e12
%
% Returns a struct with fields:
%   times            - the requested times
%   meanNumber       - m(t), the Poisson mean
%   varNumber        - equal to meanNumber, the law being Poisson
%   arrivalRate      - lambda(t)
%   departureRate    - delta(t) = E[lambda(t-S)]
%   offeredLoadPSA   - ES*lambda(t), the pointwise stationary approximation
%   meanLag          - E[Se], when ES2 is given
%   lagApproximation - ES*lambda(t-E[Se]), when ES2 is given
%
% Example:
%   % sinusoidal arrivals, exponential service of rate 2
%   res = qsys_mtginf(@(t) 10+5*sin(t), @(x) exp(-2*x), 1/2, linspace(0,2*pi,50));
%
% Reference: S. G. Eick, W. A. Massey, W. Whitt (1993). The physics of the
% Mt/G/infinity queue. Operations Research 41(4), 731-742.
%
% See also QSYS_MGINF, QSYS_ERLANGA, QSYS_GGISGI_FLUID.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

options = struct('starttime', -Inf, 'es2', [], 'servicepdf', [], 'tol', 1e-12, ...
    'panels', 4000, 'maxage', 1e12);
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

if ES <= 0
    line_error(mfilename, 'The mean service time ES must be positive.');
end
t = tvals(:).';
cut = qsys_mtginf_tailcut(serviceCcdf, options.tol, options.maxage);
unbounded = ~isfinite(options.starttime);

mean_ = qsys_mtginf_mean(lambdaFun, serviceCcdf, t, options.starttime, cut, options.panels, unbounded);
arrival = qsys_mtginf_eval(lambdaFun, t);

result.times = t;
result.meanNumber = mean_;
result.varNumber = mean_;                  % Poisson: the variance is the mean
result.arrivalRate = arrival;
result.offeredLoadPSA = ES * arrival;

if ~isempty(options.servicepdf)
    dep = zeros(1, numel(t));
    for i = 1:numel(t)
        if unbounded
            hi = cut;
        else
            hi = min(cut, max(0, t(i) - options.starttime));
        end
        [x, w] = qsys_mtginf_simpson(0, hi, options.panels);
        dep(i) = sum(w .* qsys_mtginf_eval(lambdaFun, t(i)-x) .* qsys_mtginf_eval(options.servicepdf, x));
    end
    result.departureRate = dep;
else
    % Flow balance m'(t) = lambda(t) - delta(t), differentiated centrally.
    h = 1e-5 * max(1, max(abs(t)));
    up = qsys_mtginf_mean(lambdaFun, serviceCcdf, t+h, options.starttime, cut, options.panels, unbounded);
    dn = qsys_mtginf_mean(lambdaFun, serviceCcdf, t-h, options.starttime, cut, options.panels, unbounded);
    result.departureRate = arrival - (up - dn)/(2*h);
end

if ~isempty(options.es2)
    lag = options.es2 / (2*ES);            % E[Se], the time lag
    result.meanLag = lag;
    result.lagApproximation = ES * qsys_mtginf_eval(lambdaFun, t - lag);
end
end

function m = qsys_mtginf_mean(lambdaFun, serviceCcdf, t, startTime, cut, panels, unbounded)
% The Poisson mean m(t), shared by the public entry point and by the
% finite-difference departure rate so that neither re-derives the other. With an
% infinite past the age grid does not move with t, so the service ccdf is
% evaluated once rather than once per time point.
if unbounded
    [xs, ws] = qsys_mtginf_simpson(0, cut, panels);
    gcs = qsys_mtginf_eval(serviceCcdf, xs);
end
m = zeros(1, numel(t));
for i = 1:numel(t)
    if unbounded
        x = xs; w = ws; gc = gcs;
    else
        hi = min(cut, max(0, t(i) - startTime));
        [x, w] = qsys_mtginf_simpson(0, hi, panels);
        gc = qsys_mtginf_eval(serviceCcdf, x);
    end
    % m(t) = int lambda(t-x) P(S>x) dx: the arrivals of age x still in service.
    m(i) = sum(w .* qsys_mtginf_eval(lambdaFun, t(i)-x) .* gc);
end
end

function y = qsys_mtginf_eval(f, x)
% Evaluate a user handle on a grid, accepting either a vectorized handle or a
% scalar one. Trying the vectorized call first matters: these grids have
% thousands of points.
y = f(x);
if numel(y) == 1 && numel(x) > 1
    y = y * ones(size(x));
elseif numel(y) ~= numel(x)
    y = arrayfun(f, x);
end
y = reshape(y, size(x));
end

function [x, w] = qsys_mtginf_simpson(a, b, n)
% Nodes and weights of the composite Simpson rule on an even panel count.
if mod(n, 2) == 1
    n = n + 1;
end
if b <= a
    x = a;
    w = 0;
    return
end
x = linspace(a, b, n+1);
w = ones(1, n+1);
w(2:2:end-1) = 4;
w(3:2:end-2) = 2;
w = w * (b - a)/(3*n);
end

function x = qsys_mtginf_tailcut(ccdf, tol, cap)
% Smallest doubling point at which the service ccdf is below tol.
x = 1;
while ccdf(x) > tol
    x = 2*x;
    if x > cap
        x = cap;
        return
    end
end
end
