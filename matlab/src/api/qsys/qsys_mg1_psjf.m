function [W, rho_total] = qsys_mg1_psjf(lambda, mu, cs)
% QSYS_MG1_PSJF Compute mean response time for M/G/1/PSJF queue
%
% [W, RHO] = QSYS_MG1_PSJF(LAMBDA, MU, CS) computes the mean response time
% for each job class in an M/G/1 queue with Preemptive Shortest Job First
% (PSJF) scheduling.
%
% Under PSJF, priority is based on a job's original size (not remaining size).
% Jobs with smaller original sizes always preempt jobs with larger sizes.
%
% For PSJF, the mean response time for a job of size x is given by
% (Section 3.2 of Wierman-Harchol-Balter 2003):
%
%   E[T(x)]^PSJF = (lambda * integral_0^x t^2*f(t)dt) / (2*(1-rho(x))^2)
%                  + x / (1 - rho(x))
%
% where:
%   - rho(x) = lambda * integral_0^x t*f(t)dt (truncated load)
%   - f(t) is the service time density (mixture of exponentials)
%
% For exponential service, the expected response time for class k is computed
% by integrating E[T(x)] over the exponential distribution of class k sizes.
%
% CLASSIFICATION (Wierman-Harchol-Balter 2003):
%   PSJF is "Always Unfair" - some job size is treated unfairly under all
%   loads and all service distributions.
%
% PARAMETERS:
%   lambda : Vector of arrival rates per class
%   mu     : Vector of service rates per class
%   cs     : Vector of coefficients of variation per class (cs=1 for exponential)
%
% RETURNS:
%   W   : Vector of mean response times per class
%   rho : Overall system utilization (modified for Little's law)
%
% REFERENCES:
%   - A. Wierman and M. Harchol-Balter, "Classifying scheduling policies with
%     respect to unfairness in an M/GI/1", SIGMETRICS 2003, Section 3.2.
%   - L. Kleinrock, "Queueing Systems, Volume II: Computer Applications",
%     Wiley, 1976.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Ensure inputs are column vectors
lambda = reshape(lambda, [], 1);
mu = reshape(mu, [], 1);
cs = reshape(cs, [], 1);

% Validate input lengths
if ~isequal(length(lambda), length(mu), length(cs))
    error('qsys_mg1_psjf:InvalidInput', ...
        'lambda, mu, and cs must have the same length');
end

% Validate positive values
if any(lambda <= 0) || any(mu <= 0) || any(cs < 0)
    error('qsys_mg1_psjf:InvalidInput', ...
        'lambda and mu must be positive, cs must be non-negative');
end

K = length(lambda);

% Overall utilization
rho_total = sum(lambda ./ mu);

% Stability check
if rho_total >= 1
    error('qsys_mg1_psjf:UnstableSystem', ...
        sprintf('System is unstable: utilization rho = %g >= 1', rho_total));
end

% Check if all classes are exponential (cs = 1)
if all(abs(cs - 1) < 1e-6)
    % Use specialized formula for exponential service
    W = qsys_mg1_psjf_exp(lambda, mu);
else
    % General case: use numerical integration
    W = qsys_mg1_psjf_general(lambda, mu, cs);
end

% Compute rhohat = Q/(1+Q) to match qsys convention
Q = sum(lambda .* W);
rho_total = Q / (1 + Q);

end

function W = qsys_mg1_psjf_exp(lambda, mu)
% Specialized PSJF formula for exponential service times
%
% For exponential service, we integrate E[T(x)] over the exponential
% distribution of job sizes. The service time mixture PDF is:
%   f(t) = sum_i (lambda_i/lambda_total) * mu_i * exp(-mu_i * t)
%
% The expected response time for class k is:
%   E[T_k] = integral_0^inf E[T(x)] * f_k(x) dx
%
% where f_k(x) = mu_k * exp(-mu_k * x) and
%   E[T(x)] = x/(1-rho(x)) + m2(x)/(2*(1-rho(x))^2)

K = length(lambda);
lambda_total = sum(lambda);
p = lambda / lambda_total;  % mixing probabilities

W = zeros(K, 1);

for k = 1:K
    % Integrate E[T(x)] * f_k(x) from 0 to infinity
    % Use numerical integration with appropriate upper limit
    mu_k = mu(k);

    % Upper limit: 20 mean service times covers 99.9999998% of exponential
    x_max = 20 / mu_k;

    % Response time function for a job of size x
    T_of_x = @(x) compute_psjf_response(x, lambda, mu, p, lambda_total);

    % Exponential density for class k
    f_k = @(x) mu_k * exp(-mu_k * x);

    % Integrand: E[T(x)] * f_k(x)
    integrand = @(x) T_of_x(x) .* f_k(x);

    % Numerical integration
    W(k) = integral(integrand, 0, x_max, 'RelTol', 1e-8, 'AbsTol', 1e-10);
end

end

function T = compute_psjf_response(x, lambda, mu, p, lambda_total)
% Compute PSJF response time E[T(x)] for a job of size x
%
% E[T(x)] = x/(1-rho(x)) + m2(x)/(2*(1-rho(x))^2)
%
% where:
%   rho(x) = lambda * integral_0^x t * f(t) dt
%   m2(x) = lambda * integral_0^x t^2 * f(t) dt
%
% For mixture of exponentials:
%   integral_0^x t * mu_i * exp(-mu_i*t) dt = 1/mu_i - (1/mu_i + x)*exp(-mu_i*x)
%   integral_0^x t^2 * mu_i * exp(-mu_i*t) dt = 2/mu_i^2 - (2/mu_i^2 + 2x/mu_i + x^2)*exp(-mu_i*x)

K = length(lambda);

% Handle vectorized x
T = zeros(size(x));

for j = 1:numel(x)
    xj = x(j);

    % Compute truncated first moment: integral_0^x t * f(t) dt
    m1_x = 0;
    for i = 1:K
        mu_i = mu(i);
        % integral_0^x t * mu_i * exp(-mu_i*t) dt
        int_t = 1/mu_i - (1/mu_i + xj) * exp(-mu_i * xj);
        m1_x = m1_x + p(i) * int_t;
    end

    % Truncated load
    rho_x = lambda_total * m1_x;

    % Compute truncated second moment: integral_0^x t^2 * f(t) dt
    m2_x = 0;
    for i = 1:K
        mu_i = mu(i);
        % integral_0^x t^2 * mu_i * exp(-mu_i*t) dt
        int_t2 = 2/mu_i^2 - (2/mu_i^2 + 2*xj/mu_i + xj^2) * exp(-mu_i * xj);
        m2_x = m2_x + p(i) * int_t2;
    end

    % Scale by lambda_total for m2 term
    m2_x_scaled = lambda_total * m2_x;

    % PSJF response time formula
    if rho_x >= 1
        T(j) = Inf;
    else
        T(j) = xj / (1 - rho_x) + m2_x_scaled / (2 * (1 - rho_x)^2);
    end
end

end

function W = qsys_mg1_psjf_general(lambda, mu, cs)
% General PSJF formula for non-exponential service (numerical integration)
% Uses Gamma approximation for service time distributions

K = length(lambda);
lambda_total = sum(lambda);

% For non-exponential, use class-based approximation (original formula)
% with numerical correction for variance

mean_service = 1 ./ mu;
[~, sort_idx] = sort(mean_service, 'ascend');

lambda_sorted = lambda(sort_idx);
mu_sorted = mu(sort_idx);
cs_sorted = cs(sort_idx);

rho_i = lambda_sorted ./ mu_sorted;

W_sorted = zeros(K, 1);

for k = 1:K
    x = 1 / mu_sorted(k);
    rho_x = sum(rho_i(1:k));

    m2_x = 0;
    for i = 1:k
        E_S2_i = (1 + cs_sorted(i)^2) / mu_sorted(i)^2;
        m2_x = m2_x + lambda_sorted(i) * E_S2_i;
    end

    if rho_x >= 1
        W_sorted(k) = Inf;
    else
        waiting_term = m2_x / (2 * (1 - rho_x)^2);
        service_term = x / (1 - rho_x);
        W_sorted(k) = waiting_term + service_term;
    end
end

[~, unsort_idx] = sort(sort_idx);
W = W_sorted(unsort_idx);

end
