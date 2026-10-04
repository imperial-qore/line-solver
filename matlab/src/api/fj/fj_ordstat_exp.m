function m = fj_ordstat_exp(ri, k)
% M = FJ_ORDSTAT_EXP(RI, K)
%
% Mean of the K-th smallest of N independent EXPONENTIAL branch completion
% times whose means are RI, i.e. the instant a k-of-n (quorum) join fires.
% K = N is the ordinary AND-join, the maximum, and K = 1 the minimum.
%
% With lambda_i = 1/RI(i) and m = N-K stragglers allowed,
%
%   E[X_(K)] = sum_{j=m+1..N} (-1)^(j-m-1) C(j-1,m) e_j,
%   e_j      = sum_{|S|=j} 1 / sum_{i in S} lambda_i
%
% the inclusion-exclusion identity for the order statistics of independent
% exponentials. At m = 0 it collapses to sum_j (-1)^(j-1) e_j, the classical
% expression for the maximum, TERM BY TERM: a full join therefore evaluates
% exactly as it did before this function existed.
%
% The sum has 2^N terms and its signs alternate, so it is evaluated exactly
% only while the branch count is small. Beyond MAXEXACT branches a genuine
% quorum (K < N) is evaluated by FJ_QUORUM_MOMENTS instead, whose
% Poisson-binomial recurrence adds no cancellation; a full join keeps the
% exact path at every N so that no existing result moves.
%
% Reference:
% A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
% ACM Computing Surveys 47(2), Article 17, 2014, Sec. 3 (Eq. 18-19).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

MAXEXACT = 15;

ri = ri(:)';
ri = ri(~isnan(ri) & ~isinf(ri));
n = numel(ri);
if n == 0
    m = 0;
    return
end
if k < 1 || k > n
    line_error(mfilename, 'k must satisfy 1 <= k <= n. Got k=%d, n=%d.', k, n);
end
% A branch of zero mean completes instantly: it never delays the join and it
% counts toward the quorum at once. Removing it here keeps the reciprocal
% below finite, which an exact arithmetic requires and IEEE only tolerates.
nzero = sum(ri <= 0);
if nzero > 0
    k = k - nzero;
    if k <= 0
        m = 0;
        return
    end
    ri = ri(ri > 0);
    n = numel(ri);
end
if n == 1
    m = ri(1);
    return
end

if k < n && n > MAXEXACT
    % Branch times are taken as exponential, so the variance is the square
    % of the mean.
    m = fj_quorum_moments(ri, ri.^2, k);
    return
end

lambdai = 1 ./ ri;
nstrag = n - k;
m = 0;
for j = (nstrag+1):n
    ej = sum(1 ./ sum(nchoosek(lambdai, j), 2));
    m = m + (-1)^(j-nstrag-1) * nchoosek(j-1, nstrag) * ej;
end
end
