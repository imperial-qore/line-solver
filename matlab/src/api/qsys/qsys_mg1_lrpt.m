function [W, rho_total] = qsys_mg1_lrpt(lambda, mu, cs)
% QSYS_MG1_LRPT Compute mean response time for M/G/1/LRPT queue
%
% [W, RHO] = QSYS_MG1_LRPT(LAMBDA, MU, CS) computes the mean response time
% for each job class in an M/G/1 queue with Longest Remaining Processing
% Time (LRPT) scheduling.
%
% Under LRPT, the job with the longest remaining processing time receives
% exclusive service. This is a remaining-size based policy that favors
% large jobs.
%
% For LRPT, the slowdown for a job of size x is given by
% (Section 3.2 of Wierman-Harchol-Balter 2003):
%
%   E[S(x)]^LRPT = 1/(1-rho) + lambda*E[X^2]/(2*x*(1-rho)^2)
%   E[T(x)]^LRPT = x * E[S(x)]^LRPT
%
% For exponential service, the expected response time for class k is computed
% by integrating E[T(x)] over the exponential distribution of class k sizes.
%
% CLASSIFICATION (Wierman-Harchol-Balter 2003):
%   LRPT is "Always Unfair" - it favors large jobs at the expense of
%   small jobs.
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
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Ensure inputs are column vectors
lambda = reshape(lambda, [], 1);
mu = reshape(mu, [], 1);
cs = reshape(cs, [], 1);

% Validate input lengths
if ~isequal(length(lambda), length(mu), length(cs))
    error('qsys_mg1_lrpt:InvalidInput', ...
        'lambda, mu, and cs must have the same length');
end

% Validate positive values
if any(lambda <= 0) || any(mu <= 0) || any(cs < 0)
    error('qsys_mg1_lrpt:InvalidInput', ...
        'lambda and mu must be positive, cs must be non-negative');
end

K = length(lambda);

% Overall utilization
rho_total = sum(lambda ./ mu);

% Stability check
if rho_total >= 1
    error('qsys_mg1_lrpt:UnstableSystem', ...
        sprintf('System is unstable: utilization rho = %g >= 1', rho_total));
end

% Check if all classes are exponential (cs = 1)
if all(abs(cs - 1) < 1e-6)
    % Use specialized formula for exponential service
    W = qsys_mg1_lrpt_exp(lambda, mu);
else
    % General case: use class-based approximation
    W = qsys_mg1_lrpt_general(lambda, mu, cs);
end

% Compute rhohat = Q/(1+Q) to match qsys convention
Q = sum(lambda .* W);
rho_total = Q / (1 + Q);

end

function W = qsys_mg1_lrpt_exp(lambda, mu)
% Specialized LRPT formula for exponential service times
%
% Uses the Wierman-Harchol-Balter formula:
%   E[S(x)]^LRPT = 1/(1-rho) + lambda*E[X^2]/(2*x*(1-rho)^2)
%   E[T(x)]^LRPT = x/(1-rho) + lambda*E[X^2]/(2*(1-rho)^2)
%
% For exponential service, integrate over the job size distribution:
%   E[T_k] = integral_0^inf E[T(x)] * f_k(x) dx

K = length(lambda);
lambda_total = sum(lambda);
rho_total = sum(lambda ./ mu);
p = lambda / lambda_total;  % mixing probabilities

% Compute E[X^2] for the mixture distribution
% E[X^2] = sum_i p_i * E[S_i^2] = sum_i p_i * 2/mu_i^2
E_X2 = sum(p .* 2 ./ (mu.^2));

W = zeros(K, 1);

for k = 1:K
    mu_k = mu(k);

    % Upper limit for integration: 20 mean service times
    x_max = 20 / mu_k;

    % LRPT response time for a job of size x
    T_of_x = @(x) x / (1 - rho_total) + lambda_total * E_X2 / (2 * (1 - rho_total)^2);

    % Exponential density for class k
    f_k = @(x) mu_k * exp(-mu_k * x);

    % Integrand: E[T(x)] * f_k(x)
    integrand = @(x) T_of_x(x) .* f_k(x);

    % Numerical integration
    W(k) = integral(integrand, 0, x_max, 'RelTol', 1e-8, 'AbsTol', 1e-10);
end

end

function W = qsys_mg1_lrpt_general(lambda, mu, cs)
% General LRPT formula using class-based preemptive priority approximation

K = length(lambda);
mean_service = 1 ./ mu;

% Sort classes by mean service time DESCENDING for LRPT priority
[~, sort_idx] = sort(mean_service, 'descend');

lambda_sorted = lambda(sort_idx);
mu_sorted = mu(sort_idx);

rho_i = lambda_sorted ./ mu_sorted;

W_sorted = zeros(K, 1);

for k = 1:K
    if k == 1
        rho_prev = 0;
    else
        rho_prev = sum(rho_i(1:k-1));
    end

    rho_curr = sum(rho_i(1:k));
    E_R_k = sum(lambda_sorted(1:k) ./ (mu_sorted(1:k).^2));
    W_q = E_R_k / ((1 - rho_prev) * (1 - rho_curr));
    W_sorted(k) = W_q + 1 / mu_sorted(k);
end

[~, unsort_idx] = sort(sort_idx);
W = W_sorted(unsort_idx);

end
