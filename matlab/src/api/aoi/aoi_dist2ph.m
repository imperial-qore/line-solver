function [alpha, T] = aoi_dist2ph(proc)
%AOI_DIST2PH Convert LINE process representation to PH format for aoi-fluid
%
% [alpha, T] = AOI_DIST2PH(proc)
%
% Converts LINE's {D0, D1} MAP representation to (alpha, T) PH format
% where alpha is the initial probability vector and T is the sub-generator matrix.
%
% In PH representation:
%   - alpha: Row vector of initial probabilities (sums to 1)
%   - T: Sub-generator matrix (negative diagonal, non-negative off-diagonal)
%   - Absorption rates: -T * ones(n,1)
%
% Parameters:
%   proc (cell): LINE process representation {D0, D1} where:
%       D0: Sub-generator matrix (like T in PH)
%       D1: Completion/transition matrix
%
% Returns:
%   alpha: Initial probability row vector
%   T: Sub-generator matrix
%
% See also: aoi_extract_params, solveBufferless, solveSingleBuffer

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(proc) || ~iscell(proc) || length(proc) < 2
    error('aoi_dist2ph:InvalidInput', 'proc must be a cell array {D0, D1}');
end

D0 = proc{1};
D1 = proc{2};

% Check for NaN (disabled class)
if any(isnan(D0(:))) || any(isnan(D1(:)))
    error('aoi_dist2ph:DisabledClass', 'Process contains NaN - class may be disabled');
end

% Get dimensions
n = size(D0, 1);

if n ~= size(D0, 2) || n ~= size(D1, 1) || n ~= size(D1, 2)
    error('aoi_dist2ph:DimensionMismatch', 'D0 and D1 must be square matrices of the same size');
end

% T is the sub-generator (D0 in MAP notation)
T = D0;

% Compute stationary distribution of the underlying CTMC
% The full generator is Q = D0 + D1
Q = D0 + D1;

% Check if Q is a proper generator (rows should sum to ~0)
rowSums = sum(Q, 2);
if max(abs(rowSums)) > 1e-10
    % Q is not a proper generator, normalize
    for i = 1:n
        Q(i, i) = Q(i, i) - rowSums(i);
    end
end

% Compute stationary distribution pi of Q
% pi * Q = 0 and sum(pi) = 1
% Solve [Q'; ones(1,n)] * [pi'; 1] = [zeros(n,1); 1]
A = [Q'; ones(1, n)];
b = [zeros(n, 1); 1];
pi = (A \ b)';

% Ensure pi is non-negative and normalized
pi = max(pi, 0);
pi = pi / sum(pi);

% Initial probability alpha is the stationary distribution weighted by
% the probability of starting in each phase after a completion
% For a renewal process: alpha = pi * D1 / (pi * D1 * ones)
% This gives the probability of starting in each phase
completionRates = D1 * ones(n, 1);
numerator = pi .* completionRates';

if sum(numerator) > 0
    alpha = numerator / sum(numerator);
else
    % Fallback: use stationary distribution
    alpha = pi;
end

% Ensure alpha is a row vector
alpha = alpha(:)';

% Normalize to ensure sum = 1
alpha = alpha / sum(alpha);

% Validate T is a proper sub-generator
% - Diagonal elements should be non-positive
% - Off-diagonal elements should be non-negative
% - Row sums should be non-positive (not zero, since it's a sub-generator)
for i = 1:n
    if T(i, i) > 0
        error('aoi_dist2ph:InvalidSubGenerator', 'T(%d,%d) = %g > 0, expected non-positive', i, i, T(i, i));
    end
end

end
