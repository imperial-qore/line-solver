function [W, rho] = qsys_mg1_prio(lambda, mu, cs)
% W=QSYS_MG1_PRIO(LAMBDA,MU,CS)
%
% Analyzes an M/G/1 queueing system with non-preemptive (Head-of-Line) priorities.
%
% This function computes per-class mean waiting times for multiple priority classes
% using the Pollaczek-Khinchine formula extended for priorities.
%
% For K priority classes (class 1 = highest priority), the waiting time for class k is:
%   W_k = (B_0 / (1 - sum_{i=1}^{k-1} rho_i)) * (1 / (1 - sum_{i=1}^{k} rho_i))
%
% Where:
%   - rho_i = lambda_i / mu_i (utilization of class i)
%   - B_0 = sum_{i=1}^{K} lambda_i * (1/mu_i^2) * (1 + cs_i^2)
%   - cs_i = coefficient of variation for class i service time
%
% PARAMETERS:
%   lambda : Vector of arrival rates per priority class (class 1 = highest)
%   mu     : Vector of service rates per priority class
%   cs     : Vector of coefficients of variation per priority class
%
% RETURNS:
%   W      : Vector of mean waiting times per priority class
%   rho    : Overall system utilization (scalar)
%
% EXAMPLE:
%   % Two priority classes
%   lambda = [0.2; 0.3];  % Arrival rates
%   mu = [1.0; 1.0];      % Service rates
%   cs = [0.5; 1.0];      % Coefficients of variation
%   [W, rho] = qsys_mg1_prio(lambda, mu, cs);
%
% REFERENCE:
%   Kleinrock, L., "Queueing Systems, Volume I: Theory", Wiley, 1975, Section 3.5

% Ensure inputs are column vectors
lambda = reshape(lambda, [], 1);
mu = reshape(mu, [], 1);
cs = reshape(cs, [], 1);

% Validate input lengths
if ~isequal(length(lambda), length(mu), length(cs))
    error('qsys_mg1_prio:InvalidInput', ...
        'lambda, mu, and cs must have the same length');
end

% Validate positive values
if any(lambda <= 0) || any(mu <= 0) || any(cs <= 0)
    error('qsys_mg1_prio:InvalidInput', ...
        'lambda, mu, and cs must all be positive');
end

% Compute per-class utilizations
rho_i = lambda ./ mu;

% Overall utilization
rho = sum(rho_i);

% Stability check
if rho >= 1
    error('qsys_mg1_prio:UnstableSystem', ...
        sprintf('System is unstable: utilization rho = %g >= 1', rho));
end

% Compute mean second moment of service time (used in waiting time formula)
% B_0 = sum_i lambda_i * E[S_i^2] / 2
% where E[S_i^2] = (1 + cs_i^2) / mu_i^2
B_0 = sum(lambda .* (1 + cs.^2) ./ (mu.^2)) / 2;

% Compute per-class waiting times (in queue, not including service)
K = length(lambda);
W_q = zeros(K, 1);

for k = 1:K
    % Cumulative sum of utilizations for higher priority classes (1 to k-1)
    if k == 1
        rho_prev = 0;
    else
        rho_prev = sum(rho_i(1:k-1));
    end

    % Cumulative sum including current class (1 to k)
    rho_curr = sum(rho_i(1:k));

    % Waiting time (in queue) for class k
    % W_q_k = B_0 / ((1 - rho_prev) * (1 - rho_curr))
    W_q(k) = B_0 / ((1 - rho_prev) * (1 - rho_curr));
end

% Response time = waiting time + service time
% W_k = W_q_k + 1/mu_k
W = W_q + 1 ./ mu;

% Compute rhohat = Q/(1+Q) to match qsys_mg1 convention
% Q is the mean number of customers in the system (Little's law)
% Q = sum_k lambda_k * W_k
Q = sum(lambda .* W);
rho = Q / (1 + Q);

end
