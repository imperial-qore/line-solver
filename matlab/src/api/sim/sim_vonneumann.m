function result = sim_vonneumann(x, alpha)
% OA_VONNEUMANN Von Neumann ratio test for randomness of a sequence.
%
% RESULT = OA_VONNEUMANN(X) applies the two-sided randomness test of von
% Neumann (1941) to the sequence X, at the default significance level 0.05.
%
% RESULT = OA_VONNEUMANN(X, ALPHA) uses significance level ALPHA.
%
% The statistic is the ratio of the mean square successive difference to the
% variance,
%   ratio = sum_{i=1}^{b-1} (x_{i+1}-x_i)^2 / sum_{i=1}^{b} (x_i - xbar)^2,
% with b = numel(X). Under the null hypothesis that X is i.i.d. normal the
% ratio has mean 2 and variance 4(b-2)/((b-1)(b+1)), and (ratio-2)/sd is
% asymptotically standard normal, so the two-sided p-value is
% 2(1 - Phi(|z|)). Serial correlation of either sign moves the ratio away
% from 2: positive correlation shrinks the successive differences and pushes
% the ratio below 2, negative correlation pushes it above.
%
% The null mean and variance above were confirmed by Monte Carlo over
% b = 10, 16, 24, 32, 50 to within 0.3%.
%
% Returns a struct with fields:
%   ratio    - The von Neumann ratio
%   zscore   - Standardized statistic (ratio-2)/sd
%   pvalue   - Two-sided p-value
%   reject   - true when pvalue < ALPHA, i.e. randomness is rejected
%   nobs     - Number of observations b
%   analyzer - Identifier string
%
% Examples:
%   sim_vonneumann(randn(1,50)).reject      % false, i.i.d. input
%   sim_vonneumann(cumsum(randn(1,50))).reject  % true, random walk
%
% Reference: J. von Neumann, "Distribution of the Ratio of the Mean Square
% Successive Difference to the Variance", Ann. Math. Statist. 12(4), 1941;
% L. C. Young, "Randomness in Ordered Sequences", Ann. Math. Statist. 12, 1941.
%
% See also OA_SHAPIROWILK, OA_FQUEST, OA_FIRQUEST
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(alpha)
    alpha = 0.05;
end
if ~isscalar(alpha) || ~isreal(alpha) || alpha <= 0 || alpha >= 1
    line_error(mfilename, 'alpha must be a real scalar in (0,1)');
end

x = x(:);
b = numel(x);
if b < 3
    line_error(mfilename, 'at least 3 observations are required, got %d', b);
end
if any(~isfinite(x))
    line_error(mfilename, 'the sequence must be finite');
end

den = sum((x - mean(x)).^2);
if den <= 0
    line_error(mfilename, 'the sequence is constant, the ratio is undefined');
end

ratio = sum(diff(x).^2) / den;
sd = sqrt(4 * (b - 2) / ((b - 1) * (b + 1)));
zscore = (ratio - 2) / sd;
pvalue = 2 * (1 - sim_normcdf(abs(zscore)));

result = struct('ratio', ratio, 'zscore', zscore, 'pvalue', pvalue, ...
    'reject', pvalue < alpha, 'nobs', b, 'analyzer', 'sim_vonneumann');
end
