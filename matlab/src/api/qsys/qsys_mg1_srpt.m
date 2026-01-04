function [W, rho_total] = qsys_mg1_srpt(lambda, mu, cs)
% QSYS_MG1_SRPT Compute mean response time for M/G/1/SRPT queue
%
% [W, RHO] = QSYS_MG1_SRPT(LAMBDA, MU, CS) computes the mean response time
% for each job class in an M/G/1 queue with Shortest Remaining Processing
% Time (SRPT) scheduling.
%
% For SRPT, the mean response time for a job of size x is given by the
% Schrage-Miller formula (Section 4.2 of Wierman-Harchol-Balter 2003):
%
%   E[T(x)] = lambda*integral_0^x t*F_bar(t)dt / (1-rho(x))^2
%             + integral_0^x dt/(1-rho(t))
%
% where:
%   - rho(x) = lambda * integral_0^x t*f(t)dt (truncated load)
%   - F_bar(t) = 1 - F(t) (complementary CDF)
%
% For exponential service (M/M/1/SRPT), with multiple classes having
% different service rates, the analysis reduces to preemptive priority
% queueing since the memoryless property makes remaining time distribution
% equal to the original distribution.
%
% PARAMETERS:
%   lambda : Vector of arrival rates per class
%   mu     : Vector of service rates per class
%   cs     : Vector of coefficients of variation per class (cs=1 for exponential)
%
% RETURNS:
%   W   : Vector of mean response times per class (sorted by service time)
%   rho : Overall system utilization
%
% MULTICLASS CASE:
%   For multiple classes with different service times, SRPT effectively
%   gives preemptive priority to smaller jobs. Classes are internally
%   sorted by mean service time (ascending), and priority queue analysis
%   is applied.
%
% REFERENCES:
%   - L. E. Schrage and L. W. Miller, "The queue M/G/1 with the shortest
%     remaining processing time discipline", Operations Research, 14:670-684, 1966.
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
    error('qsys_mg1_srpt:InvalidInput', ...
        'lambda, mu, and cs must have the same length');
end

% Validate positive values
if any(lambda <= 0) || any(mu <= 0) || any(cs < 0)
    error('qsys_mg1_srpt:InvalidInput', ...
        'lambda and mu must be positive, cs must be non-negative');
end

K = length(lambda);

% Compute mean service times per class
mean_service = 1 ./ mu;

% Sort classes by mean service time (ascending) for SRPT priority
[sorted_service, sort_idx] = sort(mean_service, 'ascend');

% Reorder parameters according to sorted service times
lambda_sorted = lambda(sort_idx);
mu_sorted = mu(sort_idx);
cs_sorted = cs(sort_idx);

% Compute per-class utilizations (sorted order)
rho_i = lambda_sorted ./ mu_sorted;

% Overall utilization
rho_total = sum(rho_i);

% Stability check
if rho_total >= 1
    error('qsys_mg1_srpt:UnstableSystem', ...
        sprintf('System is unstable: utilization rho = %g >= 1', rho_total));
end

% For SRPT with general service distributions, we use numerical integration.
% For the exponential case (cs = 1), SRPT reduces to preemptive priority.
% We check if all cs values are approximately 1 (exponential case).
is_exponential = all(abs(cs - 1) < 1e-10);

if is_exponential || K == 1
    % Exponential case or single class: use preemptive priority formula
    W_sorted = qsys_mg1_srpt_exp(lambda_sorted, mu_sorted);
else
    % General case: use numerical integration for each class
    W_sorted = qsys_mg1_srpt_general(lambda_sorted, mu_sorted, cs_sorted);
end

% Restore original class ordering
[~, unsort_idx] = sort(sort_idx);
W = W_sorted(unsort_idx);

% Compute rhohat = Q/(1+Q) to match qsys convention
Q = sum(lambda .* W);
rho_total = Q / (1 + Q);

end


function W = qsys_mg1_srpt_exp(lambda, mu)
% SRPT for exponential service - uses preemptive priority formula
% Because of memoryless property, SRPT behaves like preemptive priority
% where priority is based on expected service time.
%
% For preemptive resume priority, class k only sees residual service time
% from classes 1 to k (higher priority classes plus itself), since it can
% preempt all lower priority classes.

K = length(lambda);
rho_i = lambda ./ mu;

% Compute response times using preemptive priority formula
W = zeros(K, 1);

for k = 1:K
    % Cumulative utilization up to class k-1 (higher priority classes)
    if k == 1
        rho_prev = 0;
    else
        rho_prev = sum(rho_i(1:k-1));
    end

    % Cumulative utilization up to class k
    rho_curr = sum(rho_i(1:k));

    % Cumulative mean residual service time from classes 1 to k
    % E[R_k] = sum_{i=1}^k lambda_i / mu_i^2
    % Class k only sees residual from classes that can't be preempted (1 to k)
    E_R_k = sum(lambda(1:k) ./ (mu(1:k).^2));

    % Mean waiting time for class k (preemptive priority)
    % W_q_k = E[R_k] / ((1 - rho_prev) * (1 - rho_curr))
    W_q = E_R_k / ((1 - rho_prev) * (1 - rho_curr));

    % Response time = waiting time + service time
    W(k) = W_q + 1 / mu(k);
end

end


function W = qsys_mg1_srpt_general(lambda, mu, cs)
% SRPT for general service distributions
% Uses numerical integration of the Schrage-Miller formula

K = length(lambda);
W = zeros(K, 1);

% Build mixture distribution parameters
% Aggregate arrival rate and mixture probabilities
lambda_total = sum(lambda);
p = lambda / lambda_total;  % mixture probabilities

% Compute overall moments
% E[S] = sum_i p_i * E[S_i]
% E[S^2] = sum_i p_i * E[S_i^2]
mean_S = sum(p ./ mu);
var_S = sum(p .* ((1 + cs.^2) ./ mu.^2)) - mean_S^2;
E_S2 = sum(p .* ((1 + cs.^2) ./ mu.^2));

% Overall load
rho = lambda_total * mean_S;

% For each class, compute response time at the class's mean service time
% This is an approximation - exact SRPT analysis requires the full
% job size distribution, but we can use the class-conditional formula.

for k = 1:K
    x = 1 / mu(k);  % Mean service time for class k

    % Compute truncated load rho(x) and integrals numerically
    % For the mixture distribution

    % Truncated load: rho(x) = lambda_total * integral_0^x t*f(t)dt
    % where f(t) = sum_i p_i * f_i(t)

    % For simplicity, approximate using the cumulative contribution
    % from classes with smaller mean service times
    smaller_idx = (1 ./ mu) <= x;
    rho_x = sum(lambda(smaller_idx) ./ mu(smaller_idx));

    % Approximate the integral terms
    % For the waiting time component
    if rho_x >= 1
        W(k) = Inf;
    else
        % Use approximation based on truncated workload
        % Numerator integral: lambda * integral_0^x t*F_bar(t)dt
        % For exponential mixture, approximate
        numerator = 0;
        for i = 1:K
            if 1/mu(i) <= x
                % Contribution from class i
                numerator = numerator + lambda(i) / (mu(i)^2);
            end
        end

        % Second integral: integral_0^x dt/(1-rho(t))
        % Approximate using quadrature
        second_term = 0;
        n_steps = 1000;
        dt = x / n_steps;
        for step = 1:n_steps
            t = (step - 0.5) * dt;
            % Compute rho(t) at this point
            rho_t = 0;
            for i = 1:K
                % For exponential: integral_0^t s*mu*exp(-mu*s)ds
                % = 1/mu * (1 - exp(-mu*t)*(1 + mu*t))
                rho_t = rho_t + lambda(i) / mu(i) * (1 - exp(-mu(i)*t) * (1 + mu(i)*t));
            end
            if rho_t < 1
                second_term = second_term + dt / (1 - rho_t);
            end
        end

        % Schrage-Miller formula
        W(k) = numerator / (1 - rho_x)^2 + second_term;
    end
end

end
