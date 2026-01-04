function [W_bar, W_bar_n] = map_m1ps_cdfrespt(C, D, mu, x, varargin)
% MAP_M1PS_CDFRESPT Compute sojourn time distribution in MAP/M/1-PS queue
%
% W_bar = MAP_M1PS_CDFRESPT(C, D, mu, x) computes the complementary
% distribution function of the sojourn time in a MAP/M/1 processor-sharing
% queue at the points specified in x.
%
% [W_bar, W_bar_n] = MAP_M1PS_CDFRESPT(C, D, mu, x) also returns the
% conditional complementary distributions W_bar_n{i} for customers finding
% i-1 customers in the system on arrival.
%
% [...] = MAP_M1PS_CDFRESPT(C, D, mu, x, 'Param', Value) specifies optional
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
%   W_bar = map_m1ps_cdfrespt(C, D, mu, x);
%   plot(x, W_bar);
%   xlabel('x'); ylabel('Pr[W > x]');
%   title('Sojourn time distribution for M/M/1-PS');

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
    error('MAP_M1PS_CDFRESPT:DimensionMismatch', 'C and D must have the same size');
end

I = eye(M);
e = ones(M, 1);

%% Compute MAP parameters
% Stationary probability vector pi of the underlying Markov chain
% Solve: pi*(C + D) = 0, pi*e = 1
Q = C + D;
% Replace last row of Q with normalization constraint
A = Q';
A(end,:) = e';
b = zeros(M, 1);
b(end) = 1;
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
    error('MAP_M1PS_CDFRESPT:Unstable', ...
        'System is unstable (rho = %.4f >= 1)', rho);
end

%% Compute R matrix
if verbose
    fprintf('Computing R matrix...\n');
end
R = compute_R_matrix(C, D, mu);

if verbose
    fprintf('  R = %.6f\n', R);
end

%% Determine truncation point N(epsilon) for queue length
% Find minimum N such that: (1/lambda) * sum_{n=0}^N pi_0 * R^n * D * e > 1 - epsilon
% where pi_0 = pi * (I - R)
pi_0 = pi * (I - R);

% Compute spectral radius of R to estimate required N
rho_R = max(abs(eig(R)));

if rho_R >= 1
    error('MAP_M1PS_CDFRESPT:UnstableR', 'R matrix has spectral radius >= 1');
end

% Estimate N based on spectral radius:
% We need rho_R^N < epsilon * (1 - rho_R)
% N > log(epsilon * (1 - rho_R)) / log(rho_R)
if rho_R > 0
    N_estimate = ceil(log(epsilon * (1 - rho_R)) / log(rho_R));
else
    N_estimate = 1;
end

% Cap at reasonable maximum to avoid excessive computation
N_max = 10000;
N_epsilon = min(max(N_estimate, 10), N_max);

if verbose
    fprintf('  Spectral radius of R: %.6f\n', rho_R);
    fprintf('  Estimated N: %d, using N(epsilon): %d (capped at %d)\n', N_estimate, N_epsilon, N_max);
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

    % Determine truncation points L and R_upper for uniformization
    % Find L and R_upper such that: sum_{k=L}^R_upper Poisson(theta+mu, x_val) > 1 - epsilon_prime
    mean_val = (theta + mu) * x_val;

    % Use Poisson quantiles to find L and R_upper
    if mean_val > 0
        L = max(0, floor(mean_val - 10*sqrt(mean_val)));
        R_upper = ceil(mean_val + 10*sqrt(mean_val));

        % Refine to meet epsilon_prime requirement
        pmf = poisspdf(L:R_upper, mean_val);
        cumsum_pmf = sum(pmf);
        while cumsum_pmf < 1 - epsilon_prime && R_upper < 10000
            R_upper = R_upper + 10;
            pmf = poisspdf(L:R_upper, mean_val);
            cumsum_pmf = sum(pmf);
        end
    else
        L = 0;
        R_upper = 0;
    end

    K_max = R_upper;

    if verbose && idx == 1
        fprintf('  K(epsilon_prime): %d (uniformization truncation)\n', K_max);
    end

    % Compute h_{n,k} for n=0,...,N_epsilon and k=0,...,K_max
    if verbose && idx == 1
        fprintf('Computing h_{n,k} recursion...\n');
    end
    h = compute_h_recursive(C, D, mu, N_epsilon, K_max, theta);

    % Compute W_bar(x) using equation (8)
    % W_bar(x) = (1/lambda) * sum_{n=0}^{N} pi_0 * R^n * D * sum_{k=L}^{R} ...
    %            [(theta+mu)^k * x^k / k!] * exp(-(theta+mu)*x) * h_{n,k}

    theta_plus_mu = theta + mu;

    for n = 0:N_epsilon
        % Compute pi_0 * R^n * D
        if n == 0
            weight = pi_0 * D;
            R_power_n = 1;
        else
            R_power_n = R^n;
            weight = pi_0 * R_power_n * D;
        end

        % Compute sum over k
        sum_k = zeros(M, 1);
        for k = L:K_max
            % Use poisspdf to avoid numerical overflow with factorial
            poisson_term = poisspdf(k, theta_plus_mu * x_val);
            sum_k = sum_k + poisson_term * h{n+1, k+1};
        end

        term_n = (1/lambda) * weight * sum_k;
        W_bar(idx) = W_bar(idx) + term_n;

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

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

function R = compute_R_matrix(C, D, mu)
% Compute rate matrix R for MAP/M/1 queue
% Solves: D + R(C - mu*I) + mu*R^2 = O

M = size(C, 1);
I = eye(M);

% For scalar case (M=1), use quadratic formula
if M == 1
    % mu*R^2 + (C - mu)*R + D = 0
    a = mu;
    b = C - mu;
    c = D;
    discriminant = b^2 - 4*a*c;

    if discriminant < 0
        error('compute_R_matrix:NoRealSolution', 'No real solution for R');
    end

    % Two solutions
    R1 = (-b - sqrt(discriminant)) / (2*a);
    R2 = (-b + sqrt(discriminant)) / (2*a);

    % Choose minimal nonnegative solution
    if R1 >= 0 && R1 < 1
        R = R1;
    elseif R2 >= 0 && R2 < 1
        R = R2;
    else
        error('compute_R_matrix:NoValidSolution', 'No valid solution in [0,1) for R');
    end
else
    % For matrix case, use fixed-point iteration
    % From D + R*(C - mu*I) + mu*R^2 = 0, we can derive:
    % R*(mu*I - C) = D + mu*R^2
    % R = (D + mu*R^2) * inv(mu*I - C)
    % which gives the fixed-point iteration: R_new = (D + mu*R_old^2) / (mu*I - C)

    % Initial approximation
    R = -D / (C - mu*I);
    max_iter = 5000;
    tol = 1e-10;

    % Precompute (mu*I - C) for efficiency
    muI_minus_C = mu*I - C;

    for iter = 1:max_iter
        R_old = R;
        % Fixed-point iteration
        R = (D + mu*(R*R)) / muI_minus_C;

        % Check convergence
        if norm(R - R_old, inf) < tol
            break;
        end
    end

    if iter == max_iter
        warning('compute_R_matrix:SlowConvergence', ...
            'R computation did not converge within %d iterations', max_iter);
    end
end

% Ensure non-negativity (numerical errors might produce small negative values)
R(R < 0) = 0;

% Verify solution
check = D + R*(C - mu*I) + mu*(R*R);
if norm(check, inf) > 1e-4
    warning('compute_R_matrix:LargeResidual', ...
        'Large residual in R computation: ||residual|| = %e', norm(check, inf));
end

end

function h = compute_h_recursive(C, D, mu, N, K, theta)
% Recursive computation of h_{n,k} coefficients
%
% The h_{n,k} vectors satisfy the recursion (Theorem 1):
%   h_{n,0} = e (vector of ones), for n = 0, 1, ...
%   h_{n,k+1} = 1/(theta+mu) * [n*mu/(n+1) * h_{n-1,k} + (theta*I + C) * h_{n,k}
%                                + D * h_{n+1,k}]
%   where h_{-1,k} = 0 for all k

M = size(C, 1);
I = eye(M);
e = ones(M, 1);

% Pre-compute constant matrices
theta_plus_mu = theta + mu;
theta_I_plus_C = theta * I + C;

% Initialize cell array to store h_{n,k}
% h{n+1, k+1} stores h_{n,k} (shift indices by 1)
h = cell(N+1, K+1);

% Base case: h_{n,0} = e for all n
for n = 0:N
    h{n+1, 1} = e;
end

% Recursive computation: fill column by column (increasing k)
for k = 0:K-1
    for n = 0:N
        % Compute h_{n,k+1} using the recursion formula
        % h_{n,k+1} = 1/(theta+mu) * [n*mu/(n+1) * h_{n-1,k}
        %                              + (theta*I + C) * h_{n,k}
        %                              + D * h_{n+1,k}]

        term1 = zeros(M, 1);
        term2 = theta_I_plus_C * h{n+1, k+1};
        term3 = zeros(M, 1);

        % First term: n*mu/(n+1) * h_{n-1,k}
        if n > 0
            term1 = (n * mu / (n + 1)) * h{n, k+1};
        end

        % Third term: D * h_{n+1,k}
        if n < N
            term3 = D * h{n+2, k+1};
        end

        h{n+1, k+2} = (term1 + term2 + term3) / theta_plus_mu;
    end
end

end
