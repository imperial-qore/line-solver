function [W, rho_total] = qsys_mg1_setf(lambda, mu, cs)
% QSYS_MG1_SETF Compute mean response time for M/G/1/SETF queue
%
% [W, RHO] = QSYS_MG1_SETF(LAMBDA, MU, CS) computes the mean response time
% for each job class in an M/G/1 queue with SETF (Shortest Elapsed Time
% First) scheduling.
%
% SETF is the non-preemptive version of FB/LAS (Feedback/Least Attained
% Service). Under SETF:
% - Jobs are ordered by their attained service (elapsed processing time)
% - The job with the least attained service has highest priority
% - However, once a job begins service, it runs to completion (non-preemptive)
%
% The difference from FB/LAS is that arriving jobs with less attained service
% must wait until the currently serving job completes, rather than preempting it.
%
% For SETF, the mean response time follows a modified FB formula that accounts
% for the non-preemptive nature. The formula includes an additional residual
% service time term due to the non-preemptive blocking:
%
%   E[T(x)]^SETF = E[T(x)]^FB + E[R] / (1 - rho_x)
%
% where E[R] is the mean residual service time and rho_x is the truncated load.
%
% PARAMETERS:
%   lambda : Vector of arrival rates per class
%   mu     : Vector of service rates per class
%   cs     : Vector of coefficients of variation per class
%
% RETURNS:
%   W   : Vector of mean response times per class
%   rho : Overall system utilization (modified for Little's law)
%
% REFERENCES:
%   - M. Nuyens and A. Wierman, "The Foreground-Background queue: A survey",
%     Performance Evaluation, 2008.
%   - A. Wierman and M. Harchol-Balter, "Classifying scheduling policies with
%     respect to unfairness in an M/GI/1", SIGMETRICS 2003.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Ensure inputs are column vectors
lambda = reshape(lambda, [], 1);
mu = reshape(mu, [], 1);
cs = reshape(cs, [], 1);

% Validate input lengths
if ~isequal(length(lambda), length(mu), length(cs))
    error('qsys_mg1_setf:InvalidInput', ...
        'lambda, mu, and cs must have the same length');
end

% Validate positive values
if any(lambda <= 0) || any(mu <= 0) || any(cs < 0)
    error('qsys_mg1_setf:InvalidInput', ...
        'lambda and mu must be positive, cs must be non-negative');
end

K = length(lambda);

% Compute per-class utilizations
rho_i = lambda ./ mu;

% Overall utilization
rho_total = sum(rho_i);

% Stability check
if rho_total >= 1
    error('qsys_mg1_setf:UnstableSystem', ...
        sprintf('System is unstable: utilization rho = %g >= 1', rho_total));
end

% Compute mean residual service time for the mixture distribution
% E[R] = sum_i (lambda_i / lambda_total) * E[S_i^2] / (2 * E[S_i])
%      = sum_i (lambda_i / lambda_total) * (1 + cs_i^2) / (2 * mu_i)
lambda_total = sum(lambda);
E_R = 0;
for i = 1:K
    p_i = lambda(i) / lambda_total;
    E_S_i = 1 / mu(i);
    E_S2_i = (1 + cs(i)^2) / mu(i)^2;
    E_R = E_R + p_i * E_S2_i / (2 * E_S_i);
end

% Compute response times using SETF formula
% SETF adds a non-preemptive penalty term to the FB formula
W = zeros(K, 1);

for k = 1:K
    % Mean service time for this class (job size x)
    x = 1 / mu(k);

    % For SETF formula (similar to FB but with residual):
    % rho_x = lambda * integral_0^x F_bar(t)dt
    rho_x = 0;
    for i = 1:K
        if cs(i) == 1  % Exponential case
            % integral_0^x exp(-mu_i*t)dt = (1 - exp(-mu_i*x)) / mu_i
            integral_Fbar = (1 - exp(-mu(i)*x)) / mu(i);
        else
            % For non-exponential, approximate using mean and variance
            % This is an approximation; exact requires knowing distribution
            integral_Fbar = min(x, 1/mu(i));  % Bounded approximation
        end
        rho_x = rho_x + lambda(i) * integral_Fbar;
    end

    % Numerator integral: lambda * integral_0^x t*F_bar(t)dt
    numerator = 0;
    for i = 1:K
        if cs(i) == 1  % Exponential case
            % integral_0^x t*exp(-mu_i*t)dt
            % = 1/mu_i^2 * (1 - exp(-mu_i*x)*(1 + mu_i*x))
            mu_i = mu(i);
            integral_tFbar = (1 - exp(-mu_i*x)*(1 + mu_i*x)) / mu_i^2;
        else
            % Approximation for non-exponential
            integral_tFbar = min(x^2/2, 1/mu(i)^2);
        end
        numerator = numerator + lambda(i) * integral_tFbar;
    end

    % SETF formula: E[T(x)]^SETF = E[T(x)]^FB + E[R] / (1 - rho_x)
    % where E[T(x)]^FB = numerator / (1-rho_x)^2 + x / (1-rho_x)
    if rho_x >= 1
        W(k) = Inf;
    else
        % FB waiting term
        fb_waiting_term = numerator / (1 - rho_x)^2;
        % FB service term (slowdown)
        fb_service_term = x / (1 - rho_x);
        % Non-preemptive penalty: residual service time adjusted
        np_penalty = E_R / (1 - rho_x);

        W(k) = fb_waiting_term + fb_service_term + np_penalty;
    end
end

% Compute rhohat = Q/(1+Q) to match qsys convention
Q = sum(lambda .* W);
rho_total = Q / (1 + Q);

end
