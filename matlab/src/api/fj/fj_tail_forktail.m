%{ @file fj_tail_forktail.m
 %  @brief ForkTail black-box tail-latency approximation for fork-join requests
 %
 %  @author LINE Development Team
%}

%{
 % @brief ForkTail black-box tail-latency approximation for fork-join requests
 %
 % @details
 % Approximates the p-th percentile of the response time of a request that
 % forks into K parallel tasks and joins on the last of them, from the mean
 % and variance of the per-branch task response times alone. Each branch is
 % treated as a black box: its task response time is fitted by a generalized
 % exponential law
 %
 %   F_T(x) = (1 - exp(-x/beta))^alpha
 %
 % whose two parameters are matched on the branch mean and variance,
 %
 %   E[T] = beta*(psi(alpha+1) - psi(1)),  V[T] = beta^2*(psi'(1) - psi'(alpha+1)),
 %
 % and the request response time is the maximum over the branches, taken as
 % the product of the branch CDFs (exact only for independent branches):
 %
 %   F_X(x) = prod_i (1 - exp(-x/beta_i))^alpha_i,  x_p = F_X^{-1}(p).
 %
 % In the homogeneous case this inverts in closed form,
 % x_p = -beta*log(1 - p^(1/(K*alpha))).
 %
 % When the fanout itself is random, i.e. a request spawns K_i tasks with
 % probability P_i (a service whose requests touch different numbers of
 % shards), the request law is the mixture
 %
 %   F_X(x) = sum_i P_i * (1 - exp(-x/beta))^(K_i*alpha),
 %
 % which is inverted numerically.
 %
 % The approximation rests on the central limit theorem for G/G/m queues in
 % heavy traffic, so it is a HIGH-LOAD result: the reference reports errors
 % within 20% and 15% at 80% and 90% utilization respectively, and makes no
 % claim at low load, where the tail is dominated by the service law rather
 % than by queueing and the branch dependence is strongest. Use
 % fj_is_homogeneous plus the FJ_codes route when the model is in the
 % homogeneous MAP/PH/1 class, which is more accurate there; ForkTail covers
 % the heterogeneous branches and mixed service laws that route rejects.
 %
 % @par Syntax:
 % @code
 % [xp, alpha, beta] = fj_tail_forktail(ET, VT, K, p)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>ET<td>Mean task response time: a scalar (homogeneous branches) or a vector with one entry per branch
 % <tr><td>VT<td>Variance of the task response time, same shape as ET
 % <tr><td>K<td>Number of branches, or a vector of distinct fanouts when the fanout is random; ignored when ET is a vector (default 1)
 % <tr><td>p<td>Percentile, either a fraction in (0,1) or a percentage in (0,100) (default 99)
 % <tr><td>P<td>Probabilities of the fanouts in K, required when K is a vector of more than one entry
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
 % Workloads", ACM HPDC 2018, pp. 206-217.
%}
function [xp, alpha, beta] = fj_tail_forktail(ET, VT, K, p, P)

if nargin < 3 || isempty(K)
    K = 1;
end
if nargin < 4 || isempty(p)
    p = 99;
end
if nargin < 5
    P = [];
end
K = K(:)';
P = P(:)';
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
alpha = zeros(1, nbranch);
beta = zeros(1, nbranch);
for i = 1:nbranch
    [alpha(i), beta(i)] = ge_fit(ET(i), VT(i));
end

if nbranch == 1 && numel(K) > 1
    % Random fanout: mix the homogeneous request laws over the fanout
    % distribution and invert numerically
    if numel(P) ~= numel(K)
        line_error(mfilename, 'A vector of fanouts K needs a probability vector P of the same length.');
    end
    if abs(sum(P) - 1) > GlobalConstants.CoarseTol || any(P < 0)
        line_error(mfilename, 'The fanout probabilities P must be non-negative and sum to one.');
    end
    mixcdf = @(x) sum(P .* (1 - exp(-x/beta)).^(K*alpha));
    xlo = -beta * log(1 - p^(1/(min(K)*alpha)));
    xhi = -beta * log(1 - p^(1/(max(K)*alpha)));
    xp = fzero(@(x) mixcdf(x) - p, [xlo, xhi]);
    return
end

if nbranch == 1
    % Homogeneous: every branch shares (alpha,beta), so the product of K
    % identical CDFs raises the shape to K*alpha and inverts in closed form
    xp = -beta * log(1 - p^(1/(K*alpha)));
    return
end

% Inhomogeneous: solve prod_i (1-exp(-x/beta_i))^alpha_i = p. F_X is bounded
% above by the CDF of any single branch, so the request percentile is at
% least the largest branch percentile; expand from there until the root is
% bracketed.
logp = log(p);
residual = @(x) sum(alpha .* log1p(-exp(-x ./ beta))) - logp;
xlo = max(-beta .* log(1 - p.^(1 ./ alpha)));
xhi = xlo;
while residual(xhi) < 0
    xhi = 2 * xhi;
    if ~isfinite(xhi)
        line_error(mfilename, 'Could not bracket the ForkTail percentile.');
    end
end
if residual(xlo) > 0
    xlo = xlo / 2;
    while residual(xlo) > 0 && xlo > eps
        xlo = xlo / 2;
    end
end
xp = fzero(residual, [xlo, xhi]);
end

function [alpha, beta] = ge_fit(ET, VT)
% [ALPHA, BETA] = GE_FIT(ET, VT)
% Match a generalized exponential law on a mean and a variance. The squared
% coefficient of variation depends on the shape alone and decreases
% monotonically in it, so the shape is recovered by a scalar root-find on a
% logarithmic scale and the scale then follows in closed form. SCV = 1 is the
% exponential case alpha = 1, kept exact.
scv = VT / ET^2;
if abs(scv - 1) < GlobalConstants.FineTol
    alpha = 1;
else
    scvOf = @(a) (psi(1,1) - psi(1,a+1)) / (psi(0,a+1) - psi(0,1))^2;
    residual = @(la) scvOf(exp(la)) - scv;
    lo = -30; hi = 30;
    while residual(lo) < 0 && lo > -700
        lo = lo - 30;   % smaller shape -> larger SCV
    end
    while residual(hi) > 0 && hi < 700
        hi = hi + 30;   % larger shape -> smaller SCV
    end
    alpha = exp(fzero(residual, [lo, hi]));
end
beta = ET / (psi(0,alpha+1) - psi(0,1));
end
