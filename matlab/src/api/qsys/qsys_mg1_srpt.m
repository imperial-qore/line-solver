function [W, rho_total] = qsys_mg1_srpt(lambda, mu, cs)
% QSYS_MG1_SRPT Compute mean response time for M/G/1/SRPT queue
%
% [W, RHO] = QSYS_MG1_SRPT(LAMBDA, MU, CS) computes the mean response time
% for each job class in an M/G/1 queue with Shortest Remaining Processing
% Time (SRPT) scheduling.
%
% SRPT is a size-based policy: it always serves the job with the smallest
% remaining processing time, preempting whenever a shorter job arrives.
% The class-conditional response time is obtained from the Schrage-Miller
% formula (Eqs (1)-(3) of Bansal-Harchol-Balter, SIGMETRICS 2003, citing
% Schrage-Miller 1966). For a job of size x:
%
%   E[T(x)] = E[W(x)] + E[R(x)]
%   E[W(x)] = lambda*(m2(x) + x^2*(1-F(x))) / (2*(1-rho(x))^2)   (waiting)
%   E[R(x)] = integral_0^x dt/(1-rho(t))                          (residence)
%
% where f(t) is the overall (mixture) job-size density, F(t) its CDF,
%   rho(x) = lambda * integral_0^x t*f(t) dt     (load from jobs of size <= x)
%   m2(x)  = integral_0^x t^2 f(t) dt
%
% The per-class mean is E[T_r] = integral_0^inf E[T(x)] f_r(x) dx, with
% f_r the class-r size density and f = sum_r (lambda_r/lambda) f_r the
% mixture. Because E[T(x)] depends only on the job size (SRPT is size-based,
% not class-based), this integral is exact. The integrals are evaluated by
% cumulative trapezoidal quadrature on a common grid.
%
% Each class is represented by a job-size distribution matched to its
% (mean=1/mu_r, scv=cs_r^2): exponential when cs_r=1, a two-phase balanced
% hyperexponential when cs_r>1, and a Tijms mixture of Erlang-(k-1)/Erlang-k
% when cs_r<1. For the fully exponential case this reproduces the exact
% M/M/1/SRPT hyperexponential-mixture result.
%
% PARAMETERS:
%   lambda : Vector of arrival rates per class
%   mu     : Vector of service rates per class
%   cs     : Vector of coefficients of variation per class (cs=1 for exponential)
%
% RETURNS:
%   W   : Vector of mean response times per class (original class order)
%   rho : System load measure Q/(1+Q) with Q = sum(lambda.*W)
%
% REFERENCES:
%   - L. E. Schrage and L. W. Miller, "The queue M/G/1 with the shortest
%     remaining processing time discipline", Operations Research, 14:670-684, 1966.
%   - N. Bansal and M. Harchol-Balter, "Analysis of SRPT scheduling:
%     investigating unfairness", SIGMETRICS 2001, Sec. 4, Eqs (1)-(3).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Ensure inputs are column vectors
lambda = reshape(lambda, [], 1);
mu = reshape(mu, [], 1);
cs = reshape(cs, [], 1);

% Validate input lengths
if ~isequal(length(lambda), length(mu), length(cs))
    error('qsys_mg1_srpt:InvalidInput', ...
        'lambda, mu, and cs must have the same length');
end

% Validate positive values
if any(lambda <= 0) || any(mu <= 0) || any(cs < 0)
    error('qsys_mg1_srpt:InvalidInput', ...
        'lambda and mu must be positive, cs must be non-negative');
end

K = length(lambda);
lambda_total = sum(lambda);
p = lambda / lambda_total;              % size-mixture probabilities

% Overall utilization and stability check
rho_util = sum(lambda ./ mu);
if rho_util >= 1
    error('qsys_mg1_srpt:UnstableSystem', ...
        sprintf('System is unstable: utilization rho = %g >= 1', rho_util));
end

% Build per-class job-size representations and collect phase-rate bounds
fits = cell(K, 1);
rate_min = Inf;
rate_max = 0;
for r = 1:K
    fits{r} = srpt_fit(mu(r), cs(r));
    rate_min = min(rate_min, fits{r}.rate_min);
    rate_max = max(rate_max, fits{r}.rate_max);
end

% Integration grid: span 40 e-foldings of the slowest phase, resolve the
% fastest phase with at least 200 points per rate ratio.
xmax = 40.0 / rate_min;
N = min(2000000, max(20000, ceil(200.0 * rate_max / rate_min)));
x = linspace(0.0, xmax, N + 1)';
dx = x(2) - x(1);

% Mixture density f(x) and tail Fbar(x) = 1 - F(x)
fmix = zeros(N + 1, 1);
Fbar = zeros(N + 1, 1);
for r = 1:K
    fmix = fmix + p(r) * srpt_pdf(fits{r}, x);
    Fbar = Fbar + p(r) * srpt_tail(fits{r}, x);
end

% Truncated moments rho(x) and m2(x) by cumulative trapezoidal integration
rho_x = lambda_total * cumtrapz(x, x .* fmix);
m2_x = cumtrapz(x, x.^2 .* fmix);

% Guard the (1-rho(x)) factors: rho(x) -> rho_util < 1 as x -> inf
denom = max(1.0 - rho_x, 1e-12);

% Schrage-Miller waiting and residence terms, and E[T(x)]
Wait = lambda_total * (m2_x + x.^2 .* Fbar) ./ (2.0 * denom.^2);
Res = cumtrapz(x, 1.0 ./ denom);
ET = Wait + Res;

% Per-class mean response time: integrate E[T(x)] against class density
W = zeros(K, 1);
for r = 1:K
    W(r) = trapz(x, ET .* srpt_pdf(fits{r}, x));
end

% Load measure returned as rhohat = Q/(1+Q) (qsys convention)
Q = sum(lambda .* W);
rho_total = Q / (1 + Q);

end


function fit = srpt_fit(mu, cs)
% Match a job-size distribution to mean 1/mu and scv cs^2.
% Returns a struct with a type tag, parameters, and phase-rate bounds.
mean_x = 1.0 / mu;
c2 = cs^2;
if abs(c2 - 1.0) < 1e-9
    % Exponential
    fit.type = 'exp';
    fit.rate = mu;
    fit.rate_min = mu;
    fit.rate_max = mu;
elseif c2 > 1.0
    % Two-phase balanced-means hyperexponential (matches mean and scv)
    pr = 0.5 * (1.0 + sqrt((c2 - 1.0) / (c2 + 1.0)));
    r1 = 2.0 * pr * mu;
    r2 = 2.0 * (1.0 - pr) * mu;
    fit.type = 'h2';
    fit.p = pr;
    fit.r1 = r1;
    fit.r2 = r2;
    fit.rate_min = min(r1, r2);
    fit.rate_max = max(r1, r2);
else
    % Tijms mixture of Erlang-(k-1) and Erlang-k with common rate
    k = ceil(1.0 / c2);
    pr = (1.0 / (1.0 + c2)) * (k * c2 - sqrt(k * (1.0 + c2) - k * k * c2));
    rate = (k - pr) / mean_x;
    fit.type = 'erlmix';
    fit.k = k;
    fit.p = pr;
    fit.rate = rate;
    fit.rate_min = rate;
    fit.rate_max = rate;
end
end


function y = srpt_pdf(fit, x)
% Job-size probability density evaluated on the grid x.
switch fit.type
    case 'exp'
        y = fit.rate * exp(-fit.rate * x);
    case 'h2'
        y = fit.p * fit.r1 * exp(-fit.r1 * x) ...
            + (1.0 - fit.p) * fit.r2 * exp(-fit.r2 * x);
    otherwise % 'erlmix'
        y = fit.p * erlang_pdf(fit.k - 1, fit.rate, x) ...
            + (1.0 - fit.p) * erlang_pdf(fit.k, fit.rate, x);
end
end


function y = srpt_tail(fit, x)
% Complementary CDF P(X > x) evaluated on the grid x.
switch fit.type
    case 'exp'
        y = exp(-fit.rate * x);
    case 'h2'
        y = fit.p * exp(-fit.r1 * x) + (1.0 - fit.p) * exp(-fit.r2 * x);
    otherwise % 'erlmix'
        y = fit.p * erlang_tail(fit.k - 1, fit.rate, x) ...
            + (1.0 - fit.p) * erlang_tail(fit.k, fit.rate, x);
end
end


function y = erlang_pdf(n, rate, x)
% Erlang-n (shape n, given rate) density, computed in log space to avoid
% overflow when rate*x is large. n=0 is a point mass at 0 (density 0 for x>0).
% f(x) = rate * pois(n-1; rate*x), pois(m;t) = exp(-t) t^m / m!.
if n <= 0
    y = zeros(size(x));
else
    t = rate * x;
    m = n - 1;
    logp = m * log(t) - t - gammaln(m + 1);   % log Poisson(m; t)
    logp(t <= 0) = -Inf;                        % t=0 => pois(m>0)=0
    if m == 0
        logp(t <= 0) = 0;                       % pois(0;0)=1
    end
    y = rate * exp(logp);
end
end


function y = erlang_tail(n, rate, x)
% Erlang-n complementary CDF P(X>x) = sum_{j=0}^{n-1} exp(-rate x)(rate x)^j/j!
% (upper Poisson tail), computed in log space. n=0 tail is 0 for x>0.
if n <= 0
    y = zeros(size(x));
else
    t = rate * x;
    y = zeros(size(x));
    for j = 0:(n - 1)
        logp = j * log(t) - t - gammaln(j + 1);
        logp(t <= 0) = -Inf;
        if j == 0
            logp(t <= 0) = 0;
        end
        y = y + exp(logp);
    end
end
end
