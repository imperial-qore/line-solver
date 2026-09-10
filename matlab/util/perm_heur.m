function p_est = perm_heur(A)
% PERM_HEUR  Heuristic approximation to the permanent of a positive matrix.
%
% P_EST = PERM_HEUR(A) computes an approximate permanent of matrix A using
% a heuristic method based on Sinkhorn scaling and mean-field approximation.
%
% The algorithm:
% 1. Uses Sinkhorn scaling to make A approximately doubly stochastic
% 2. Applies a mean-field approximation with van der Waerden bounds
% 3. Combines with a Gurvits-like capacity bound
% 4. Scales the result back to the original matrix scale
%
% Input:
%   A - Strictly positive matrix (all elements must be > 0; a zero entry is
%       refused by PERM_REQUIRE_SUPPORT rather than floored)
%
% Output:
%   p_est - Approximate permanent value
%
% Note: This is a heuristic approximation suitable for large matrices where
% exact computation is computationally prohibitive. The approximation quality
% depends on the structure of the input matrix.

% Validate input
if any(A(:) < 0)
    error('perm_heur:InvalidInput', 'Matrix must be non-negative.');
end

n = size(A, 1);
tol = 1e-10;
maxiter = 1000;

% A zero used to be replaced by 1e-15 here. That is not invertible: it changes
% the permanent by n!*eps, which is O(1) by n=18. Refuse instead.
perm_require_support(A, mfilename);
B = A;

% Sinkhorn scaling to make matrix approximately doubly stochastic
r = ones(n, 1);
c = ones(n, 1);

converged = false;
for iter = 1:maxiter
    r = 1 ./ (B * c);
    c = 1 ./ (B' * r);
    if max(abs(r .* (B * c) - 1)) < tol
        converged = true;
        break;
    end
end
if ~converged
    line_error(mfilename, ['The Sinkhorn scaling did not converge to a doubly ' ...
        'stochastic matrix in %d sweeps (margin error %g against a tolerance of ' ...
        '%g). The estimate below assumes convergence, so no value is returned. ' ...
        'The usual cause is a matrix without total support.'], ...
        maxiter, max(abs(r .* (B * c) - 1)), tol);
end
B = diag(r) * B * diag(c);

% Compute approximate permanent of doubly stochastic matrix B
% Mean-field approximation: product of row sums / n^n * n!
row_prod = prod(sum(B, 2));
p_meanfield = factorial(n) * (row_prod / n^n);

% Gurvits-like capacity bound (optional refinement)
cap = exp(sum(log(sum(B, 2))) / n);
p_gurvits = factorial(n) * (cap / n)^n;

% Combine (simple average)
p_est = 0.5 * (p_meanfield + p_gurvits);

% Undo scaling
scale_factor = prod(1 ./ r) * prod(1 ./ c);
p_est = p_est * scale_factor;
end
