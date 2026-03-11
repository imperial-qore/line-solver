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
theta_plus_mu = theta + mu;

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

%% Determine global K_max across all x points
K_max_global = 0;
L_per_point = zeros(1, num_points);
K_per_point = zeros(1, num_points);
for idx = 1:num_points
    x_val = x(idx);
    mean_val = theta_plus_mu * x_val;
    if mean_val > 0
        L_per_point(idx) = max(0, floor(mean_val - 10*sqrt(mean_val)));
        R_upper = ceil(mean_val + 10*sqrt(mean_val));
        pmf = poisspdf(L_per_point(idx):R_upper, mean_val);
        cumsum_pmf = sum(pmf);
        while cumsum_pmf < 1 - epsilon_prime && R_upper < 10000
            R_upper = R_upper + 10;
            pmf = poisspdf(L_per_point(idx):R_upper, mean_val);
            cumsum_pmf = sum(pmf);
        end
        K_per_point(idx) = R_upper;
    end
    K_max_global = max(K_max_global, K_per_point(idx));
end

if verbose
    fprintf('  K_max (global uniformization truncation): %d\n', K_max_global);
    fprintf('Computing h_{n,k} recursion...\n');
end

%% Precompute h_{n,k} using 3D matrix (M x (N+1) x (K+1)) for fast indexing
h = compute_h_matrix(C, D, mu, M, N_epsilon, K_max_global, theta);

%% Precompute pi_0 * R^n * D weights with early truncation
% Use iterative R^n multiplication and stop when weight is negligible
weight_tol = epsilon * 1e-2;
weights = zeros(N_epsilon + 1, M);  % weights(n+1,:) = pi_0 * R^n * D (row vector)
R_power = I;  % R^0 = I
N_actual = N_epsilon;
for n = 0:N_epsilon
    w = pi_0 * R_power * D;
    weights(n+1, :) = w;
    if n > 0 && norm(w, inf) < weight_tol
        N_actual = n;
        break;
    end
    R_power = R_power * R;
end

if verbose
    fprintf('  N_actual (early truncation): %d / %d\n', N_actual, N_epsilon);
end

%% Compute sojourn time distribution for each x
for idx = 1:num_points
    x_val = x(idx);
    L = L_per_point(idx);
    K_max = K_per_point(idx);

    % Precompute Poisson PMF values for this x point
    poisson_vals = poisspdf(L:K_max, theta_plus_mu * x_val);  % vector

    for n = 0:N_actual
        weight = weights(n+1, :);  % 1 x M row vector

        % Vectorized sum over k: sum_k(m) = sum_j poisson(j) * h(m, n+1, L+j+1)
        h_slice = squeeze(h(:, n+1, L+1:K_max+1));  % M x (K_max-L+1)
        if M == 1
            h_slice = h_slice(:)';  % Ensure row vector for M=1
        end
        sum_k = h_slice * poisson_vals(:);  % M x 1

        term_n = (1/lambda) * weight * sum_k;
        W_bar(idx) = W_bar(idx) + term_n;

        % Store conditional distribution if requested
        if nargout > 1
            W_bar_n{n+1}(idx) = sum(sum_k) / M;
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

function h = compute_h_matrix(C, D, mu, M, N, K, theta)
% Compute h_{n,k} coefficients as a 3D matrix h(M, N+1, K+1)
%
% The h_{n,k} vectors satisfy the recursion (Theorem 1):
%   h_{n,0} = e (vector of ones), for n = 0, 1, ...
%   h_{n,k+1} = 1/(theta+mu) * [n*mu/(n+1) * h_{n-1,k} + (theta*I + C) * h_{n,k}
%                                + D * h_{n+1,k}]
%   where h_{-1,k} = 0 for all k

e = ones(M, 1);
theta_plus_mu = theta + mu;
theta_I_plus_C = theta * eye(M) + C;
inv_tpm = 1 / theta_plus_mu;

% Store h as M x (N+1) x (K+1) matrix
h = zeros(M, N+1, K+1);

% Base case: h_{n,0} = e for all n
for n = 0:N
    h(:, n+1, 1) = e;
end

% Precompute n*mu/(n+1) coefficients
nmu_coeff = zeros(N+1, 1);
for n = 1:N
    nmu_coeff(n+1) = n * mu / (n + 1);
end

% Recursive computation: fill column by column (increasing k)
for k = 0:K-1
    % Vectorized over n: extract column k for all n values
    h_col_k = h(:, :, k+1);  % M x (N+1)

    % Compute (theta*I + C) * h_{n,k} for all n at once
    term2_all = theta_I_plus_C * h_col_k;  % M x (N+1)

    % Compute D * h_{n+1,k} for n=0..N-1
    term3_all = D * h_col_k(:, 2:end);  % M x N (for n=0..N-1)

    for n = 0:N
        term2 = term2_all(:, n+1);

        % First term: n*mu/(n+1) * h_{n-1,k}
        if n > 0
            term1 = nmu_coeff(n+1) * h_col_k(:, n);
        else
            term1 = zeros(M, 1);
        end

        % Third term: D * h_{n+1,k}
        if n < N
            term3 = term3_all(:, n+1);
        else
            term3 = zeros(M, 1);
        end

        h(:, n+1, k+2) = inv_tpm * (term1 + term2 + term3);
    end
end

end
