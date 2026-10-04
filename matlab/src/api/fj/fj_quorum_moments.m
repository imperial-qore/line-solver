function [m, v] = fj_quorum_moments(branchMeans, branchVars, k)
% FJ_QUORUM_MOMENTS Mean and variance of a k-of-n (quorum) join completion time.
%
% [M, V] = FJ_QUORUM_MOMENTS(BRANCHMEANS, BRANCHVARS, K) returns the first two
% moments of the K-th smallest of N independent branch completion times, each
% branch being given by its mean and variance.
%
% Each branch is expanded into a discrete step CDF by a two-moment fit. At each
% time the number of completed branches is Poisson-binomial, so its distribution
% is built by the recurrence q_j <- q_{j-1}*F_i + q_j*(1-F_i) over branches i, and
% the K-th order statistic is the upper tail sum_{j>=K} q_j. For K = N this
% reduces to the product of the branch CDFs, i.e. the ordinary AND-join, and for
% K = 1 to 1 - prod(1 - F_i), i.e. the minimum.
%
% This is the same quantity as the inclusion-exclusion identity
%   F_(k)(t) = sum_{i=k..n} (-1)^(i-k) * C(i-1, k-1) * e_i(F_1(t), ..., F_n(t))
% used by LQNS, but every term here is a probability in [0,1] and none is
% subtracted, so it avoids the catastrophic cancellation the alternating binomial
% sum incurs as N grows.
%
% Follows the formulation of Omari, Franks, Woodside and Pan, as implemented in
% LQNS 6.x (randomvar.cc). Mirrors jline.api.fj.FJ_quorum in the Java runtime.

MAX_BRANCHES = 512;

n = numel(branchMeans);
if n ~= numel(branchVars)
    line_error(mfilename, 'branchMeans and branchVars must have the same number of elements.');
end
if n == 0
    m = 0; v = 0; return;
end
if k < 1 || k > n
    line_error(mfilename, sprintf('k must satisfy 1 <= k <= n. Got k=%d, n=%d.', k, n));
end
if n > MAX_BRANCHES
    % Evaluation is cubic in the branch count; fail immediately rather than hang.
    line_error(mfilename, sprintf(['quorum join has %d branches, above the supported maximum ' ...
        'of %d; evaluation is cubic in the branch count.'], n, MAX_BRANCHES));
end

% Fit each branch to a three-point discrete distribution matching its moments.
branches = cell(1, n);
for i = 1:n
    branches{i} = local_threePointFit(branchMeans(i), branchVars(i));
    if isempty(branches{i}.t)
        % A branch with no mass completes instantly, so that it neither delays the
        % join nor suppresses the completion counts below.
        branches{i} = struct('t', 0, 'A', 1);
    end
end

% Evaluate the order statistic pointwise on the union of the branch grids.
grid = [];
for i = 1:n
    grid = [grid, branches{i}.t]; %#ok<AGROW>
end
grid = unique(grid);

ngrid = numel(grid);

% Tabulate each branch along the merged grid by a single monotone walk, so the
% evaluation below is linear rather than quadratic in the grid size.
fv = zeros(n, ngrid);
for i = 1:n
    bt = branches{i}.t; bA = branches{i}.A;
    p = 1; cur = 0;
    for g = 1:ngrid
        while p <= numel(bt) && bt(p) <= grid(g)
            cur = bA(p);
            p = p + 1;
        end
        fv(i, g) = cur;
    end
end

total = struct('t', grid, 'A', zeros(1, ngrid));
for g = 1:ngrid
    q = zeros(1, n + 1);
    q(1) = 1;   % q(j+1) holds P(exactly j branches complete)
    for i = 1:n
        f = fv(i, g);
        for j = min(i, n):-1:1
            q(j + 1) = q(j) * f + q(j + 1) * (1 - f);
        end
        q(1) = q(1) * (1 - f);
    end
    total.A(g) = sum(q((k + 1):(n + 1)));
end

m = local_mean(total);
v = local_variance(total, m);
end

function v = local_at(f, x)
% Value of the step function at time x.
v = 0;
for i = 1:numel(f.t)
    if f.t(i) > x
        break;
    end
    v = f.A(i);
end
end

function f = local_threePointFit(mu, var)
% Two-moment fit to a three-point discrete distribution.
if mu < 0 || var < 0
    line_error(mfilename, 'mean and variance must be non-negative.');
end
if mu == 0
    f = struct('t', [], 'A', []); return;
end
sd = sqrt(var);
if sd == 0
    f = struct('t', mu, 'A', 1); return;
end
if mu > sd
    t1 = mu - sd;
else
    t1 = 0;
end
t2 = mu;
if sd >= mu
    t3 = mu + 2 * var / mu;
else
    t3 = mu + 2 * sd;
end
delta = t1^2 * (t3 - t2) + t2^2 * (t1 - t3) + t3^2 * (t2 - t1);
if delta == 0
    % The abscissae are not distinct, so the fit is undetermined.
    f = struct('t', mu, 'A', 1); return;
end
temp = var + mu^2;
a1 = (temp * (t3 - t2) + t2^2 * (mu - t3) + t3^2 * (t2 - mu)) / delta;
a3 = (t1^2 * (mu - t2) + t2^2 * (t1 - mu) + temp * (t2 - t1)) / delta;
f = struct('t', [t1, t2, t3], 'A', [a1, 1 - a3, 1]);
end

function m = local_mean(f)
m = 0; prev = 0;
for i = 1:numel(f.t)
    m = m + (f.A(i) - prev) * f.t(i);
    prev = f.A(i);
end
end

function v = local_variance(f, m)
v = 0; prev = 0;
for i = 1:numel(f.t)
    v = v + (f.A(i) - prev) * (f.t(i) - m)^2;
    prev = f.A(i);
end
v = max(v, 0);
end

