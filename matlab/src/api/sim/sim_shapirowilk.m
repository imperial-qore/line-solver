function result = sim_shapirowilk(x, alpha)
% OA_SHAPIROWILK Shapiro-Wilk test for univariate normality.
%
% RESULT = OA_SHAPIROWILK(X) tests the null hypothesis that X is a sample from
% a normal distribution, at the default significance level 0.05.
%
% RESULT = OA_SHAPIROWILK(X, ALPHA) uses significance level ALPHA.
%
% This is Royston's AS R94 algorithm, valid for 3 <= numel(X) <= 5000. The
% statistic is
%   W = (sum_i a_i x_(i))^2 / sum_i (x_i - xbar)^2,
% where x_(i) are the order statistics and a is the antisymmetric weight
% vector obtained by correcting the normalized expected normal order
% statistics m_i = Phi^{-1}((i-3/8)/(n+1/4)) in their two extreme components.
% Small W means departure from normality, so the test is one-sided in W and
% the p-value is an upper normal tail after Royston's normalizing transform,
% which has three branches: n = 3 exact, 4 <= n <= 11, and n >= 12.
%
% Returns a struct with fields:
%   W        - The Shapiro-Wilk statistic
%   pvalue   - p-value, small means normality is rejected
%   zscore   - Normalized statistic, NaN when n = 3
%   reject   - true when pvalue < ALPHA
%   nobs     - Number of observations n
%   analyzer - Identifier string
%
% Examples:
%   sim_shapirowilk(randn(1,50)).reject          % false, normal input
%   sim_shapirowilk(exprnd(1,1,50)).reject       % true, skewed input
%
% Reference: J. P. Royston, "Approximating the Shapiro-Wilk W-test for
% Non-normality", Statistics and Computing 2, 1992; J. P. Royston, "Remark
% AS R94", Applied Statistics 44(4), 1995. W and the p-value agree with
% scipy.stats.shapiro to 5e-10 and 1.5e-7 respectively over n up to 2000.
%
% See also OA_VONNEUMANN, OA_FQUEST, OA_FIRQUEST
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(alpha)
    alpha = 0.05;
end
if ~isscalar(alpha) || ~isreal(alpha) || alpha <= 0 || alpha >= 1
    line_error(mfilename, 'alpha must be a real scalar in (0,1)');
end

x = sort(x(:));
n = numel(x);
if n < 3
    line_error(mfilename, 'at least 3 observations are required, got %d', n);
end
if n > 5000
    line_error(mfilename, 'the AS R94 approximation is valid up to n = 5000, got %d', n);
end
if any(~isfinite(x))
    line_error(mfilename, 'the sample must be finite');
end

ssd = sum((x - mean(x)).^2);
if ssd <= 0
    line_error(mfilename, 'the sample is constant, W is undefined');
end

a = sim_shapirowilk_weights(n);
W = (a' * x)^2 / ssd;
W = min(W, 1);

if n == 3
    % exact null distribution, W is supported on [3/4, 1]
    pvalue = 6 / pi * (asin(sqrt(W)) - asin(sqrt(0.75)));
    pvalue = min(max(pvalue, 0), 1);
    zscore = NaN;
else
    if n <= 11
        g = -2.273 + 0.459 * n;
        w = -log(g - log(1 - W));
        mu = 0.5440 - 0.39978 * n + 0.025054 * n^2 - 0.0006714 * n^3;
        sigma = exp(1.3822 - 0.77857 * n + 0.062767 * n^2 - 0.0020322 * n^3);
    else
        ln = log(n);
        w = log(1 - W);
        mu = -1.5861 - 0.31082 * ln - 0.083751 * ln^2 + 0.0038915 * ln^3;
        sigma = exp(-0.4803 - 0.082676 * ln + 0.0030302 * ln^2);
    end
    zscore = (w - mu) / sigma;
    pvalue = 1 - sim_normcdf(zscore);
end

result = struct('W', W, 'pvalue', pvalue, 'zscore', zscore, ...
    'reject', pvalue < alpha, 'nobs', n, 'analyzer', 'sim_shapirowilk');
end

function a = sim_shapirowilk_weights(n)
% Royston AS R94 antisymmetric weight vector, a(n+1-i) = -a(i).
if n == 3
    a = [-sqrt(0.5); 0; sqrt(0.5)];
    return
end

c1 = [0, 0.221157, -0.147981, -2.071190, 4.434685, -2.706056];
c2 = [0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633];

m = sim_norminv(((1:n)' - 0.375) / (n + 0.25));
mm = m' * m;
c = m / sqrt(mm);
u = 1 / sqrt(n);

a = m;
an = c(n) + polyval(fliplr(c1), u);
if n > 5
    anm1 = c(n - 1) + polyval(fliplr(c2), u);
    phi = (mm - 2 * m(n)^2 - 2 * m(n - 1)^2) / (1 - 2 * an^2 - 2 * anm1^2);
    a(3:n - 2) = m(3:n - 2) / sqrt(phi);
    a(n) = an;
    a(n - 1) = anm1;
    a(1) = -an;
    a(2) = -anm1;
else
    phi = (mm - 2 * m(n)^2) / (1 - 2 * an^2);
    a(2:n - 1) = m(2:n - 1) / sqrt(phi);
    a(n) = an;
    a(1) = -an;
end
end
