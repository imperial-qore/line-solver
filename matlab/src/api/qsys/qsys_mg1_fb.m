function [W, rho_total] = qsys_mg1_fb(lambda, mu, cs)
% QSYS_MG1_FB Compute mean response time for M/G/1/FB (Feedback/LAS) queue
%
% [W, RHO] = QSYS_MG1_FB(LAMBDA, MU, CS) computes the mean response time
% for each job class in an M/G/1 queue with Feedback (FB) scheduling,
% also known as Least Attained Service (LAS) or Shortest Elapsed Time (SET).
%
% Under FB/LAS, the job with the least attained service (smallest age)
% receives priority. This is an age-based policy where priority depends
% on how much service a job has received, not its original or remaining size.
%
% For FB, the mean response time for a job of size x is given by
% (Section 3.3 of Wierman-Harchol-Balter 2003):
%
%   E[T(x)]^FB = (lambda * integral_0^x t*F_bar(t)dt) / (1-rho_x)^2
%                + x / (1 - rho_x)
%
% where:
%   - rho_x = lambda * integral_0^x F_bar(t)dt (load from jobs completing
%             service before reaching age x)
%   - F_bar(t) = 1 - F(t) (complementary CDF / survival function)
%
% For exponential service, the expected response time for class k is computed
% by integrating E[T(x)] over the exponential distribution of class k sizes.
%
% CLASSIFICATION (Wierman-Harchol-Balter 2003):
%   FB is "Always Unfair" - some job size is treated unfairly under all
%   loads and all service distributions. However, FB approximates SRPT
%   for heavy-tailed distributions and is practical since job sizes
%   need not be known in advance.
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
%     respect to unfairness in an M/GI/1", SIGMETRICS 2003, Section 3.3.
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
    error('qsys_mg1_fb:InvalidInput', ...
        'lambda, mu, and cs must have the same length');
end

% Validate positive values
if any(lambda <= 0) || any(mu <= 0) || any(cs < 0)
    error('qsys_mg1_fb:InvalidInput', ...
        'lambda and mu must be positive, cs must be non-negative');
end

K = length(lambda);

% Overall utilization
rho_total = sum(lambda ./ mu);

% Stability check
if rho_total >= 1
    error('qsys_mg1_fb:UnstableSystem', ...
        sprintf('System is unstable: utilization rho = %g >= 1', rho_total));
end

% Check if all classes are exponential (cs = 1)
if all(abs(cs - 1) < 1e-6)
    % Use specialized formula for exponential service
    W = qsys_mg1_fb_exp(lambda, mu);
else
    % General case: use original approximation
    W = qsys_mg1_fb_general(lambda, mu, cs);
end

% Compute rhohat = Q/(1+Q) to match qsys convention
Q = sum(lambda .* W);
rho_total = Q / (1 + Q);

end

function W = qsys_mg1_fb_exp(lambda, mu)
% Specialized FB formula for exponential service times
%
% For exponential service, we integrate E[T(x)] over the exponential
% distribution of job sizes. The service time mixture PDF is:
%   f(t) = sum_i (lambda_i/lambda_total) * mu_i * exp(-mu_i * t)
%
% The survival function for the mixture is:
%   F_bar(t) = sum_i (lambda_i/lambda_total) * exp(-mu_i * t)
%
% The expected response time for class k is:
%   E[T_k] = integral_0^inf E[T(x)] * f_k(x) dx
%
% where f_k(x) = mu_k * exp(-mu_k * x) and
%   E[T(x)] = numerator(x)/(1-rho_x)^2 + x/(1-rho_x)

K = length(lambda);
lambda_total = sum(lambda);
p = lambda / lambda_total;  % mixing probabilities

W = zeros(K, 1);

for k = 1:K
    % Integrate E[T(x)] * f_k(x) from 0 to infinity
    mu_k = mu(k);

    % Upper limit: 20 mean service times covers 99.9999998% of exponential
    x_max = 20 / mu_k;

    % Response time function for a job of size x
    T_of_x = @(x) compute_fb_response(x, lambda, mu, p, lambda_total);

    % Exponential density for class k
    f_k = @(x) mu_k * exp(-mu_k * x);

    % Integrand: E[T(x)] * f_k(x)
    integrand = @(x) T_of_x(x) .* f_k(x);

    % Numerical integration
    W(k) = integral(integrand, 0, x_max, 'RelTol', 1e-8, 'AbsTol', 1e-10);
end

end

function T = compute_fb_response(x, lambda, mu, p, lambda_total)
% Compute FB/LAS response time E[T(x)] for a job of size x
%
% E[T(x)] = numerator(x) / (1-rho_x)^2 + x / (1-rho_x)
%
% where:
%   rho_x = lambda * integral_0^x F_bar(t) dt
%   numerator(x) = lambda * integral_0^x t * F_bar(t) dt
%
% For mixture of exponentials:
%   F_bar(t) = sum_i p_i * exp(-mu_i * t)
%   integral_0^x exp(-mu_i*t) dt = (1 - exp(-mu_i*x)) / mu_i
%   integral_0^x t*exp(-mu_i*t) dt = (1 - exp(-mu_i*x)*(1 + mu_i*x)) / mu_i^2

K = length(lambda);

% Handle vectorized x
T = zeros(size(x));

for j = 1:numel(x)
    xj = x(j);

    % Compute rho_x = lambda * integral_0^x F_bar(t) dt
    rho_x = 0;
    for i = 1:K
        mu_i = mu(i);
        % integral_0^x exp(-mu_i*t) dt = (1 - exp(-mu_i*x)) / mu_i
        int_Fbar = (1 - exp(-mu_i * xj)) / mu_i;
        rho_x = rho_x + p(i) * lambda_total * int_Fbar;
    end

    % Compute numerator = lambda * integral_0^x t * F_bar(t) dt
    numerator = 0;
    for i = 1:K
        mu_i = mu(i);
        % integral_0^x t*exp(-mu_i*t) dt = (1 - exp(-mu_i*x)*(1 + mu_i*x)) / mu_i^2
        int_tFbar = (1 - exp(-mu_i * xj) * (1 + mu_i * xj)) / mu_i^2;
        numerator = numerator + p(i) * lambda_total * int_tFbar;
    end

    % FB response time formula
    if rho_x >= 1
        T(j) = Inf;
    else
        T(j) = numerator / (1 - rho_x)^2 + xj / (1 - rho_x);
    end
end

end

function W = qsys_mg1_fb_general(lambda, mu, cs)
% General FB formula for non-exponential service (deterministic class sizes)

K = length(lambda);

W = zeros(K, 1);

for k = 1:K
    x = 1 / mu(k);

    rho_x = 0;
    for i = 1:K
        if abs(cs(i) - 1) < 1e-6  % Exponential case
            integral_Fbar = (1 - exp(-mu(i)*x)) / mu(i);
        else
            integral_Fbar = min(x, 1/mu(i));
        end
        rho_x = rho_x + lambda(i) * integral_Fbar;
    end

    numerator = 0;
    for i = 1:K
        if abs(cs(i) - 1) < 1e-6  % Exponential case
            mu_i = mu(i);
            integral_tFbar = (1 - exp(-mu_i*x)*(1 + mu_i*x)) / mu_i^2;
        else
            integral_tFbar = min(x^2/2, 1/mu(i)^2);
        end
        numerator = numerator + lambda(i) * integral_tFbar;
    end

    if rho_x >= 1
        W(k) = Inf;
    else
        waiting_term = numerator / (1 - rho_x)^2;
        service_term = x / (1 - rho_x);
        W(k) = waiting_term + service_term;
    end
end

end
