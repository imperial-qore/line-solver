function [pi, L] = qbd_finite_solve(A_minus, A_0, A_plus, N)
% QBD_FINITE_SOLVE Stationary distribution of a finite-capacity QBD process.
%
%   [PI, L] = QBD_FINITE_SOLVE(A_MINUS, A_0, A_PLUS, N)
%
%   Solves a homogeneous finite-capacity QBD with buffer size N.
%   The generator has the block-tridiagonal structure:
%
%       Q = [ B_0   A_plus                          ]
%           [ A_minus A_0   A_plus                   ]
%           [         A_minus A_0   A_plus            ]
%           [                  ...   ...   ...        ]
%           [                       A_minus  B_N     ]
%
%   where B_0 = A_0 + A_minus (absorption of downward transitions at 0)
%   and   B_N = A_0 + A_plus  (absorption of upward transitions at N).
%
%   Uses block LDU elimination (forward reduction + back substitution)
%   to solve pi*Q = 0, sum(pi) = 1 in O(N*m^3) time where m is the
%   phase dimension.
%
%   Inputs:
%     A_minus  - m x m sub-diagonal block (downward transitions)
%     A_0      - m x m diagonal block (internal transitions, negative diagonal)
%     A_plus   - m x m super-diagonal block (upward transitions)
%     N        - buffer capacity (levels 0, 1, ..., N)
%
%   Outputs:
%     pi       - cell array of length N+1, pi{k+1} is the 1 x m
%                stationary probability vector at level k
%     L        - mean number in system (sum over k of k * sum(pi{k+1}))
%
% Copyright (c) 2026, Imperial College London. BSD 3-Clause License.

m = size(A_0, 1);

% Boundary blocks
B_0 = A_0 + A_minus;   % level 0: no downward transitions
B_N = A_0 + A_plus;    % level N: no upward transitions

% Forward reduction (block Gaussian elimination)
% Eliminate sub-diagonal: transform Q into upper block-bidiagonal form.
% F{k} stores the modified diagonal block at level k.
% G{k} stores the super-diagonal block at level k.
F = cell(N + 1, 1);
G = cell(N + 1, 1);

F{1} = B_0;
G{1} = A_plus;

for k = 1:N-1
    % Multiplier: A_minus * inv(F{k})
    M = A_minus / F{k};
    F{k + 1} = A_0 - M * A_plus;
    G{k + 1} = A_plus;
end

% Last level
if N > 0
    M = A_minus / F{N};
    F{N + 1} = B_N - M * A_plus;
else
    F{1} = A_0; % single level, no transitions possible
end

% Solve pi * [upper bidiagonal] = 0 with normalization.
% Last level: pi_N is in the left null space of F{N+1}.
% Replace first column of F{N+1}' with normalization placeholder,
% then back-substitute.

% Find pi_N from the nullspace of F{N+1}
% pi_N * F{N+1} = 0 => F{N+1}' * pi_N' = 0
[~, ~, V] = svd(F{N + 1}.');
pi_N = V(:, end).';
pi_N = pi_N / sum(pi_N); % temporary normalization

% Back substitution: pi_{k} * F{k} + pi_{k+1} * A_minus = 0
% => pi_{k} = -pi_{k+1} * A_minus * inv(F{k})
pi = cell(N + 1, 1);
pi{N + 1} = pi_N;

for k = N:-1:1
    pi{k} = -pi{k + 1} * A_minus / F{k};
end

% Normalize so that sum of all probabilities = 1
total = 0;
for k = 1:(N + 1)
    total = total + sum(pi{k});
end
for k = 1:(N + 1)
    pi{k} = pi{k} / total;
end

% Mean number in system
L = 0;
for k = 1:(N + 1)
    L = L + (k - 1) * sum(pi{k});
end
end
