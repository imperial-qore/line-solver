function h = map_m1ps_h_recursive(C, D, mu, N, K)
% MAP_M1PS_H_RECURSIVE Recursive computation of h_{n,k} coefficients
%
% h = MAP_M1PS_H_RECURSIVE(C, D, mu, N, K) computes the h_{n,k} vectors
% used in the sojourn time distribution for a MAP/M/1-PS queue.
%
% The h_{n,k} vectors satisfy the recursion (Theorem 1 in paper):
%   h_{n,0} = e (vector of ones), for n = 0, 1, ...
%   h_{n,k+1} = 1/(theta+mu) * [n*mu/(n+1) * h_{n-1,k} + (theta*I + C) * h_{n,k}
%                                + D * h_{n+1,k}]
%   where h_{-1,k} = 0 for all k
%
% Input:
%   C   - M x M matrix governing MAP transitions without arrivals
%   D   - M x M matrix governing MAP transitions with arrivals
%   mu  - Service rate (scalar)
%   N   - Maximum value of n to compute (determines rows)
%   K   - Maximum value of k to compute (determines columns)
%
% Output:
%   h   - Cell array of size (N+1) x (K+1)
%         h{n+1, k+1} contains the M x 1 vector h_{n,k}
%         (indices shifted by 1 for MATLAB 1-based indexing)
%
% Reference:
%   Masuyama, H., & Takine, T. (2003). Sojourn time distribution in a
%   MAP/M/1 processor-sharing queue. Operations Research Letters, 31(6),
%   406-412.
%
% See also: MAP_M1PS_SOJOURN, MAP_COMPUTE_R

% Copyright (c) 2012-2025, Imperial College London
% All rights reserved.

M = size(C, 1);
I = eye(M);
e = ones(M, 1);

% Compute theta (uniformization parameter)
theta = max(abs(diag(C)));

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
