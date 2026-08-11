function [W_bar, W_bar_n] = map_m1ps_sojourn(C, D, mu, x, varargin)
% MAP_M1PS_SOJOURN Compute sojourn time distribution in MAP/M/1-PS queue
%
% W_bar = MAP_M1PS_SOJOURN(C, D, mu, x) computes the complementary
% distribution function of the sojourn time in a MAP/M/1 processor-sharing
% queue at the points specified in x.
%
% [W_bar, W_bar_n] = MAP_M1PS_SOJOURN(C, D, mu, x) also returns the
% conditional complementary distributions W_bar_n{i} for customers finding
% i-1 customers in the system on arrival.
%
% [...] = MAP_M1PS_SOJOURN(C, D, mu, x, 'Param', Value) specifies optional
% parameters:
%   'Epsilon'       - Truncation parameter for queue length (default: 1e-11)
%   'EpsilonPrime'  - Truncation parameter for uniformization (default: 1e-10)
%   'Verbose'       - Display computation progress (default: false)
%
% Input:
%   C   - M x M matrix governing MAP transitions without arrivals
%   D   - M x M matrix governing MAP transitions with arrivals
%   mu  - Service rate (scalar, mu > 0)
%   x   - Vector of time points at which to evaluate W_bar(x) = Pr[W > x]
%
% Output:
%   W_bar    - Vector of same size as x, containing Pr[W > x]
%   W_bar_n  - Cell array of conditional distributions (optional)
%
% The processor-sharing (PS) discipline shares the server equally among
% all customers. When n customers are present, each receives service at
% rate 1/n.
%
% The algorithm implements Theorem 1 from:
%   Masuyama, H., & Takine, T. (2003). Sojourn time distribution in a
%   MAP/M/1 processor-sharing queue. Operations Research Letters, 31(6),
%   406-412.
%
% Example:
%   % M/M/1-PS queue with lambda=0.8, mu=1
%   lambda = 0.8; mu = 1;
%   C = -lambda; D = lambda;  % Poisson arrivals
%   x = linspace(0, 10, 100);
%   W_bar = map_m1ps_sojourn(C, D, mu, x);
%   plot(x, W_bar);
%   xlabel('x'); ylabel('Pr[W > x]');
%   title('Sojourn time distribution for M/M/1-PS');
%
% See also: MAP_M1PS_H_RECURSIVE, MAP_COMPUTE_R

% Copyright (c) 2012-2025, Imperial College London
% All rights reserved.

%% Parse input arguments
p = inputParser;
addRequired(p, 'C', @(x) ismatrix(x) && size(x,1) == size(x,2));
addRequired(p, 'D', @(x) ismatrix(x) && size(x,1) == size(x,2));
addRequired(p, 'mu', @(x) isscalar(x) && x > 0);
addRequired(p, 'x', @(x) isnumeric(x) && all(x >= 0));
addParameter(p, 'Epsilon', 1e-11, @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'EpsilonPrime', 1e-10, @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'Verbose', false, @islogical);

parse(p, C, D, mu, x, varargin{:});
epsilon = p.Results.Epsilon;
epsilon_prime = p.Results.EpsilonPrime;
verbose = p.Results.Verbose;

%% Validate inputs
M = size(C, 1);
if size(D, 1) ~= M || size(D, 2) ~= M
    error('MAP_M1PS_SOJOURN:DimensionMismatch', 'C and D must have the same size');
end

I = eye(M);
e = ones(M, 1);

%% Compute MAP parameters
% Stationary probability vector pi of the underlying Markov chain
% Solve: pi*(C + D) = 0, pi*e = 1
Q = C + D;
A = [Q; e'];
b = [zeros(M, 1); 1];
pi = (A \ b)';

% Mean arrival rate
lambda = pi * D * e;

if verbose
    fprintf('MAP/M/1-PS Sojourn Time Computation\n');
    fprintf('  M (states): %d\n', M);
    fprintf('  lambda (arrival rate): %.4f\n', lambda);
    fprintf('  mu (service rate): %.4f\n', mu);
    fprintf('  rho (utilization): %.4f\n', lambda/mu);
end

% Check stability condition
rho = lambda / mu;
if rho >= 1
    error('MAP_M1PS_SOJOURN:Unstable', ...
        'System is unstable (rho = %.4f >= 1)', rho);
end

%% Compute R matrix
if verbose
    fprintf('Computing R matrix...\n');
end
R = map_compute_R(C, D, mu);

%% Determine truncation point N(epsilon) for queue length
% Find minimum N such that: (1/lambda) * sum_{n=0}^N pi_0 * R^n * D * e > 1 - epsilon
% where pi_0 = pi * (I - R)
pi_0 = pi * (I - R);

cumsum_prob = 0;
N_epsilon = 0;
for n = 0:1000
    cumsum_prob = cumsum_prob + (1/lambda) * pi_0 * (R^n) * D * e;
    if cumsum_prob > 1 - epsilon
        N_epsilon = n;
        break;
    end
end

if N_epsilon == 0
    N_epsilon = 100;  % Default fallback
    warning('MAP_M1PS_SOJOURN:TruncationDefault', ...
        'Could not determine N(epsilon), using default N=%d', N_epsilon);
end

if verbose
    fprintf('  N(epsilon): %d (queue length truncation)\n', N_epsilon);
end

%% Compute uniformization parameter theta
theta = max(abs(diag(C)));

%% Initialize output
x = x(:)';  % Ensure row vector
num_points = length(x);
W_bar = zeros(1, num_points);

if nargout > 1
    W_bar_n = cell(N_epsilon + 1, 1);
    for n = 0:N_epsilon
        W_bar_n{n+1} = zeros(1, num_points);
    end
end

%% Compute sojourn time distribution for each x
for idx = 1:num_points
    x_val = x(idx);

    % Determine truncation points L and R for uniformization
    % Find L and R such that: sum_{k=L}^R Poisson(theta+mu, x_val) > 1 - epsilon_prime
    mean_val = (theta + mu) * x_val;

    % Use Poisson quantiles to find L and R
    if mean_val > 0
        L = max(0, floor(mean_val - 10*sqrt(mean_val)));
        R = ceil(mean_val + 10*sqrt(mean_val));

        % Refine to meet epsilon_prime requirement
        pmf = poisspdf(L:R, mean_val);
        cumsum_pmf = sum(pmf);
        while cumsum_pmf < 1 - epsilon_prime && R < 10000
            R = R + 10;
            pmf = poisspdf(L:R, mean_val);
            cumsum_pmf = sum(pmf);
        end
    else
        L = 0;
        R = 0;
    end

    K_max = R;

    if verbose && idx == 1
        fprintf('  K(epsilon_prime): %d (uniformization truncation)\n', K_max);
    end

    % Compute h_{n,k} for n=0,...,N_epsilon and k=0,...,K_max
    if verbose && idx == 1
        fprintf('Computing h_{n,k} recursion...\n');
    end
    h = map_m1ps_h_recursive(C, D, mu, N_epsilon, K_max);

    % Compute W_bar(x) using equation (8)
    % W_bar(x) = (1/lambda) * sum_{n=0}^{N} pi_0 * R^n * D * sum_{k=L}^{R} ...
    %            [(theta+mu)^k * x^k / k!] * exp(-(theta+mu)*x) * h_{n,k}

    theta_plus_mu = theta + mu;
    exp_factor = exp(-theta_plus_mu * x_val);

    for n = 0:N_epsilon
        % Compute pi_0 * R^n * D
        if n == 0
            weight = pi_0 * D;
        else
            weight = pi_0 * (R^n) * D;
        end

        % Compute sum over k
        sum_k = zeros(M, 1);
        for k = L:K_max
            poisson_term = ((theta_plus_mu * x_val)^k / factorial(k)) * exp_factor;
            sum_k = sum_k + poisson_term * h{n+1, k+1};
        end

        W_bar(idx) = W_bar(idx) + (1/lambda) * weight * sum_k;

        % Store conditional distribution if requested
        if nargout > 1
            W_bar_n{n+1}(idx) = sum(sum_k) / M;  % Average over states
        end
    end
end

if verbose
    fprintf('Computation complete.\n');
end

end
