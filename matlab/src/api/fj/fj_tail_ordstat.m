%{ @file fj_tail_ordstat.m
 %  @brief Tail latency of a k-of-n (quorum) fork-join request
 %
 %  @author LINE Development Team
%}

%{
 % @brief Tail latency of a k-of-n (quorum) fork-join request
 %
 % @details
 % Approximates the p-th percentile of the response time of a request that
 % forks into N parallel tasks and joins on the K-th of them, from the mean
 % and variance of the per-branch task response times alone. K = N is the
 % ordinary AND-join and reproduces FJ_TAIL_FORKTAIL exactly; K = 1 is the
 % first completion.
 %
 % Each branch is treated as the same black box FJ_TAIL_FORKTAIL uses: its
 % task response time is fitted by a generalized exponential law
 %
 %   F_i(x) = (1 - exp(-x/beta_i))^alpha_i
 %
 % matched on the branch mean and variance by GE_FIT. The request completes
 % once K of the N branches have, so its law is the K-th ORDER STATISTIC of
 % independent, not identically distributed branch times,
 %
 %   F_X(x) = P(at least K of the N branches are done by x),
 %
 % which is the upper tail of a Poisson-binomial with success probabilities
 % F_i(x). It is evaluated by the standard convolution recurrence, which adds
 % no cancellation, and inverted by bisection. With homogeneous branches the
 % recurrence collapses to the regularized incomplete beta function
 % I_{F(x)}(K, N-K+1), used directly.
 %
 % At K = N both routes reduce TERM BY TERM to prod_i F_i(x), the product of
 % the branch CDFs that FJ_TAIL_FORKTAIL inverts, so a full join evaluates
 % exactly as it did before this function existed.
 %
 % BRANCH INDEPENDENCE is assumed, as in ForkTail: the branches of one
 % request are positively correlated through their shared arrival instant, so
 % the true quorum percentile is somewhat larger than this one. The same
 % heavy-traffic caveat applies, see FJ_TAIL_FORKTAIL.
 %
 % @par Syntax:
 % @code
 % [xp, alpha, beta] = fj_tail_ordstat(ET, VT, K, p, kreq)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>ET<td>Mean task response time: a scalar (homogeneous branches) or a vector with one entry per branch
 % <tr><td>VT<td>Variance of the task response time, same shape as ET
 % <tr><td>K<td>Number of branches; used only when ET is a scalar (default 1)
 % <tr><td>p<td>Percentile, either a fraction in (0,1) or a percentage in (0,100) (default 99)
 % <tr><td>kreq<td>Quorum: the join fires on the kreq-th branch (default: every branch)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>xp<td>Predicted p-th percentile of the request response time
 % <tr><td>alpha<td>Fitted shape parameter(s) of the generalized exponential
 % <tr><td>beta<td>Fitted scale parameter(s)
 % </table>
 %
 % @par References:
 % M. Nguyen, S. Alesawi, N. Li, H. Che, H. Jiang, "ForkTail: A Black-Box
 % Fork-Join Tail Latency Prediction Model for User-Facing Datacenter
 % Workloads", ACM HPDC 2018, pp. 206-217, for the branch law.
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems", ACM
 % Computing Surveys 47(2), Article 17, 2014, Sec. 3, for the quorum.
%}
function [xp, alpha, beta] = fj_tail_ordstat(ET, VT, K, p, kreq)

if nargin < 3 || isempty(K)
    K = 1;
end
if nargin < 4 || isempty(p)
    p = 99;
end
if p > 1
    p = p / 100;
end
if p <= 0 || p >= 1
    line_error(mfilename, 'The percentile must lie strictly between 0 and 1 (or 0 and 100).');
end

ET = ET(:)';
VT = VT(:)';
if numel(VT) ~= numel(ET)
    line_error(mfilename, 'ET and VT must have the same number of entries.');
end
if any(ET <= 0) || any(VT <= 0)
    line_error(mfilename, 'The task response time mean and variance must be positive.');
end

nbranch = numel(ET);
if nbranch == 1
    nsib = max(1, round(K(1)));
else
    nsib = nbranch;
end
if nargin < 5 || isempty(kreq)
    kreq = nsib;
end
kreq = round(kreq);
if kreq < 1 || kreq > nsib
    line_error(mfilename, 'The quorum must satisfy 1 <= kreq <= %d. Got %d.', nsib, kreq);
end

alpha = zeros(1, nbranch);
beta = zeros(1, nbranch);
for i = 1:nbranch
    [alpha(i), beta(i)] = ge_fit(ET(i), VT(i));
end

% The AND-join percentile is an upper bound for every quorum, and the SINGLE
% branch percentile a lower one, so the two bracket the root without a search.
xhi = fj_tail_forktail(ET, VT, nsib, p);
if kreq == nsib
    xp = xhi;
    return
end
xlo = min(-beta .* log(1 - p .^ (1 ./ alpha)));

if nbranch == 1
    % Homogeneous: the count done by x is Binomial(nsib, F(x)), so the
    % quorum CDF is the regularized incomplete beta of its upper tail.
    cdfk = @(x) betainc(gecdf(x, alpha, beta), kreq, nsib - kreq + 1);
else
    cdfk = @(x) poissbin_upper(gecdf(x, alpha, beta), kreq);
end

residual = @(x) cdfk(x) - p;
if residual(xlo) > 0
    % A single branch already meets the percentile; the quorum is met earlier.
    while residual(xlo) > 0 && xlo > eps
        xlo = xlo / 2;
    end
end
xp = fzero(residual, [xlo, xhi]);
end

function u = gecdf(x, alpha, beta)
% U = GECDF(X, ALPHA, BETA)
% The generalized-exponential CDF (1-exp(-x/beta))^alpha, per branch, clamped
% into [0,1] so that a rounding excursion cannot leave the beta/binomial
% routines out of domain.
u = exp(alpha .* log1p(-exp(-x ./ beta)));
u(~isfinite(u)) = 0;
u = min(1, max(0, u));
end

function q = poissbin_upper(u, kreq)
% Q = POISSBIN_UPPER(U, KREQ)
% P(at least KREQ successes) for independent Bernoulli trials of success
% probabilities U, by the convolution recurrence over the trials. All the
% terms are non-negative, so no cancellation is introduced.
n = numel(u);
pmf = zeros(1, n+1);
pmf(1) = 1;                 % pmf(j+1) = P(j successes so far)
for i = 1:n
    pmf(2:i+1) = pmf(2:i+1)*(1-u(i)) + pmf(1:i)*u(i);
    pmf(1) = pmf(1)*(1-u(i));
end
q = sum(pmf(kreq+1:end));
end
