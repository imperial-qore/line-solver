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
% Stationary probability vector pi of the underlying Markov chain,
% pi*(C+D) = 0 with pi*e = 1. This is the LEFT null vector of Q: stacking Q
% itself rather than Q' solves Q*z = 0, whose solution is e/M for any
% generator, so the phase law came back uniform and lambda with it (0.85
% against the true 0.785714 on a two-phase MAP, which map_lambda already
% reported correctly).
Q = C + D;
pi = map_prob({C, D});
pi = pi(:)';

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

    % Determine truncation points K_lo and K_max for uniformization: the
    % smallest window with sum_{k=K_lo}^{K_max} Poisson(theta+mu, x_val)
    % > 1 - epsilon_prime.
    %
    % The upper index MUST NOT be called R: R is the rate matrix computed above
    % and is read again below as pi_0*R^n*D. Shadowing it made R^n a scalar
    % power, so W_bar(0) came back as 0.2 instead of 1 and the tail grew to
    % 1e+200 -- a complementary distribution outside [0,1] at every point.
    mean_val = (theta + mu) * x_val;

    if mean_val > 0
        K_lo = max(0, floor(mean_val - 10*sqrt(mean_val)));
        K_max = ceil(mean_val + 10*sqrt(mean_val));

        % Refine to meet epsilon_prime requirement
        cumsum_pmf = sum(poisspdf(K_lo:K_max, mean_val));
        while cumsum_pmf < 1 - epsilon_prime && K_max < 10000
            K_max = K_max + 10;
            cumsum_pmf = sum(poisspdf(K_lo:K_max, mean_val));
        end
    else
        K_lo = 0;
        K_max = 0;
    end

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

    Rpow = I;
    for n = 0:N_epsilon
        % weight = pi_0 * R^n * D, accumulated so R^n costs one product
        weight = pi_0 * Rpow * D;

        % Compute sum over k
        sum_k = zeros(M, 1);
        for k = K_lo:K_max
            poisson_term = exp(k*log(theta_plus_mu * x_val) - gammaln(k+1)) * exp_factor;
            if x_val == 0
                poisson_term = double(k == 0);
            end
            sum_k = sum_k + poisson_term * h{n+1, k+1};
        end

        W_bar(idx) = W_bar(idx) + (1/lambda) * weight * sum_k;

        % Conditional CCDF given the arrival finds n customers: the phase law at
        % such an arrival epoch is weight/(weight*e), so that
        % sum_n P(N=n)*W_bar_n = W_bar with P(N=n) = weight*e/lambda.
        if nargout > 1
            wsum = weight * ones(M, 1);
            if wsum > 0
                W_bar_n{n+1}(idx) = (weight * sum_k) / wsum;
            else
                W_bar_n{n+1}(idx) = 0;
            end
        end

        Rpow = Rpow * R;
    end
end

if verbose
    fprintf('Computation complete.\n');
end

end
