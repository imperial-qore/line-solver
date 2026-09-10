function stats = sim_sts_quantile_areas(Y, b, m, p, weight)
% OA_STS_QUANTILE_AREAS Standardized time series areas of the batched quantile process.
%
% STATS = OA_STS_QUANTILE_AREAS(Y, B, M, P) splits the B*M observations in Y
% into B nonoverlapping batches of size M and returns the signed standardized
% time series (STS) areas of the quantile-estimation process, the batched
% quantile estimators, and the three variance-parameter estimators that
% OA_FQUEST and OA_FIRQUEST build confidence intervals from.
%
% STATS = OA_STS_QUANTILE_AREAS(Y, B, M, P, WEIGHT) uses the constant STS
% weight function WEIGHT instead of the default sqrt(12). The requirement on a
% weight function w is that int_0^1 w(t)B(t)dt be standard normal for a
% standard Brownian bridge B; for a constant w = c that variance is c^2/12, so
% c = sqrt(12) is the normalizing choice.
%
% With yhat_p(j,m) the empirical P-quantile of batch J and yhat_p(j,k) the
% empirical P-quantile of its first K observations, the STS process of batch J
% is
%   T_{j,m}(k/m) = (k/sqrt(m)) (yhat_p(j,m) - yhat_p(j,k)),
% its signed area is
%   A_p(w;j,m) = m^{-1} sum_{k=1}^{m} w(k/m) T_{j,m}(k/m),
% and the three variance-parameter estimators of
% sigma_p^2 = lim n Var(ytilde_p(n)) are
%   A_p(w;b,m) = b^{-1} sum_j A_p(w;j,m)^2                    (STS area)
%   N_p(b,m)   = (b-1)^{-1} m sum_j (yhat_p(j,m)-ytilde_p(n))^2  (NBQ)
%   V_p(w;b,m) = [b A_p(w;b,m) + (b-1) N_p(b,m)] / (2b-1)      (combined)
% where ytilde_p(n) is the full-sample empirical P-quantile over all n = B*M
% observations. The first two have limiting chi-square laws on B and B-1
% degrees of freedom and are asymptotically independent, so the combined
% estimator carries 2B-1 degrees of freedom and is about sqrt(2) less variable
% than either component.
%
% The prefix quantiles yhat_p(j,k) are exact order statistics, obtained from a
% Fenwick tree over the within-batch ranks that is advanced across all B
% batches simultaneously, so the cost is O(B*M log M) with the M loop carrying
% only vectorized statements.
%
% Returns a struct with fields:
%   areas    - B x 1 signed STS areas A_p(w;j,m)
%   bqe      - B x 1 batched quantile estimators yhat_p(j,m)
%   quantile - Full-sample empirical P-quantile ytilde_p(n), n = B*M
%   Ap       - Batched STS area estimator A_p(w;b,m)
%   Np       - NBQ variance-parameter estimator N_p(b,m)
%   Vp       - Combined variance-parameter estimator V_p(w;b,m)
%   b        - Batch count
%   m        - Batch size
%   n        - Number of observations used, B*M
%   analyzer - Identifier string
%
% Examples:
%   s = sim_sts_quantile_areas(exprnd(1, 32000, 1), 32, 1000, 0.9);
%   s.Vp    % estimates p(1-p)/f(y_p)^2 = 9 for i.i.d. Exp(1) at p = 0.9
%
% Reference: C. Alexopoulos, D. Goldsman, A. Lolos, K. D. Dingec, J. R. Wilson,
% "Steady-State Quantile Estimation Using Standardized Time Series", 2020/2023;
% A. Lolos et al., Proc. Winter Simulation Conference, 2023, theorems 1-3.
%
% See also OA_FQUEST, OA_FIRQUEST
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5 || isempty(weight)
    weight = sqrt(12);
end
if ~isscalar(b) || b ~= floor(b) || b < 1
    line_error(mfilename, 'the batch count b must be a positive integer');
end
if ~isscalar(m) || m ~= floor(m) || m < 1
    line_error(mfilename, 'the batch size m must be a positive integer');
end
if ~isscalar(p) || ~isreal(p) || p <= 0 || p >= 1
    line_error(mfilename, 'p must be a real scalar in (0,1)');
end
if ~isscalar(weight) || ~isreal(weight) || weight == 0
    line_error(mfilename, 'weight must be a nonzero real scalar');
end

Y = Y(:);
n = b * m;
if numel(Y) ~= n
    line_error(mfilename, 'Y must hold exactly b*m = %d observations, got %d', ...
        n, numel(Y));
end
if any(~isfinite(Y))
    line_error(mfilename, 'the sample path must be finite');
end

Ym = reshape(Y, m, b);
[sorted, ord] = sort(Ym, 1);

colOff = (0:b - 1) * m;
rnk = zeros(m, b);
rnk(ord + colOff) = repmat((1:m)', 1, b);

bqe = sorted(ceil(m * p), :);

% Fenwick tree per batch over the within-batch ranks 1..m, advanced in
% lockstep across batches so the k loop stays vectorized.
F = zeros(m, b);
LOG = floor(log2(m));
acc = zeros(1, b);

for k = 1:m
    pos = rnk(k, :);
    while true
        act = pos <= m;
        if ~any(act)
            break
        end
        pa = pos(act);
        lin = pa + colOff(act);
        F(lin) = F(lin) + 1;
        pos(act) = pa + (pa - bitand(pa, pa - 1));
    end

    L = ceil(p * k);
    pos = zeros(1, b);
    rem = L * ones(1, b);
    step = 2^LOG;
    while step >= 1
        cand = pos + step;
        oki = find(cand <= m);
        if ~isempty(oki)
            fv = reshape(F(cand(oki) + colOff(oki)), 1, []);
            mv = fv < rem(oki);
            sel = oki(mv);
            if ~isempty(sel)
                rem(sel) = rem(sel) - fv(mv);
                pos(sel) = cand(sel);
            end
        end
        step = step / 2;
    end

    qk = sorted(pos + 1 + colOff);
    acc = acc + k * (bqe - qk);
end

areas = weight * acc / (m * sqrt(m));
allSorted = sort(Y);
quantile = allSorted(ceil(n * p));

Ap = mean(areas.^2);
if b >= 2
    Np = m * sum((bqe - quantile).^2) / (b - 1);
    Vp = (b * Ap + (b - 1) * Np) / (2 * b - 1);
else
    % a single batch carries no between-batch degrees of freedom; areas and bqe
    % stay valid and OA_FIRQUEST pools them across replications instead
    Np = NaN;
    Vp = NaN;
end

stats = struct('areas', areas(:), 'bqe', bqe(:), 'quantile', quantile, ...
    'Ap', Ap, 'Np', Np, 'Vp', Vp, 'b', b, 'm', m, 'n', n, ...
    'analyzer', 'sim_sts_quantile_areas');
end
