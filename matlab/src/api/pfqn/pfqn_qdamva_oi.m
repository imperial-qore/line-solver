function [X, Qoi, Qli, Qdelay, it] = pfqn_qdamva_oi(Z, N, mu, Dli, options)
% [X, QOI, QLI, QDELAY, IT] = PFQN_QDAMVA_OI(Z, N, MU, DLI, OPTIONS)
%
% Queue-dependent approximate mean value analysis (QD-AMVA) of a closed
% product-form network composed of an aggregated infinite-server (delay) node,
% any number of load-independent (LI) single-server queues, and any number of
% order-independent (OI) stations with empty swap graph. This is the
% approximate counterpart of PFQN_MVAOI: it returns the same quantities from a
% fixed point of M*R equations whose cost is INDEPENDENT of the population,
% whereas PFQN_MVAOI is exact but recurs over a lattice of size
% prod_r nchoosek(N_r+K+1, K+1).
%
% The algorithm specializes to OI stations the QD-AMVA method of Casale, Perez
% and Wang, "QD-AMVA: evaluating systems with queue-dependent service
% requirements" (Performance Evaluation 91:80-98, 2015). Three steps take the
% exact OI mean-value expression to the fixed point solved here.
%
% 1) EXACT OI MEAN-VALUE FORM. Let Phi_i be the balanced-fairness balance
%    function of OI station i, the function tabulated by PFQN_NCOI, and let
%      mu_{i,r}(n) = Phi_i(n - e_r) / Phi_i(n)
%    be the class-r departure rate of station i in state n. Matching the BCMP
%    queue-dependent template pi(n) = G^{-1} prod_i C_i(n_i) F_i(n_i) with the
%    multinomial C_i(n) = |n|!/prod_r n_r! gives the queue-dependent demand
%    D_{i,r}(n) = F_i(n)/F_i(n - e_r) = n_r / (|n| mu_{i,r}(n)), and the exact
%    state-dependent MVA identity (eq. 4 of the reference)
%      x_{i,r}(N) = T_r sum_n |n| D_{i,r}(n) pi(n - e_r | N - e_r)
%    collapses, under n = m + e_r, to the OI arrival theorem
%      x_{i,r}(N) = T_r * m_r / mu_{i,r}(m),   m = e_r + x_i(N - e_r),
%    the |m| of the multinomial factor cancelling the 1/|n| of D. The OI
%    residence time is therefore Little's law inside the station: the jobs of
%    class r found on arrival, the tagged one included, over the rate at which
%    the station clears class r.
%
% 2) CLASS SPLIT OF THE RATE. The balanced-fairness recursion is equivalent to
%      sum_{r: n_r > 0} mu_{i,r}(n) = mu_i(n),
%    which holds exactly but is a single equation in R unknowns, so the split
%    of the station rate among the classes must be closed. Two solvable
%    extremes pin the closure down:
%      mu_i(n) = sum_r g_r(n_r)   (separable)  =>  Phi_i factorizes over the
%                                                  classes and mu_{i,r} = g_r(n_r);
%      mu_i(n) = M(|n|)           (shared)     =>  mu_{i,r} = (n_r/|n|) M(|n|).
%    The Aumann-Shapley decomposition of the rate along the ray from the empty
%    state to n,
%      mu_i(n) - mu_i(0) = sum_r n_r * int_0^1 (d mu_i / d n_r)(t n) dt,
%    reproduces BOTH exactly, is an exact additive decomposition of mu_i for
%    any differentiable rate, and is adopted here as the closure
%      mu_{i,r}(m) = mu_i(m) * w_r / sum_s w_s,
%      w_r = m_r * int_0^1 (d mu_i / d n_r)(t m) dt.
%    The integral is evaluated by Gauss-Legendre quadrature, so the cost is
%    independent of the population. When every w_r vanishes, which happens only
%    for a rate that is constant in n, the split falls back to w_r = m_r, the
%    exact shared-rate answer. Verified against exact lattice values, the
%    shares are correct to machine precision on separable rates, on
%    total-population rates including the multiserver rate min(c, |n|), and at
%    every state with a single job per class.
%
% 3) SCHWEITZER INTERPOLATION (eq. 7 of the reference), applied componentwise
%    to the vector argument,
%      x_{i,r}(N - e_r) = delta_r x_{i,r},  delta_r = (N_r - 1)/N_r,
%      x_{i,s}(N - e_r) = x_{i,s},          s ~= r,
%    so that m_r = 1 + delta_r x_{i,r} and m_s = x_{i,s}. The reference also
%    offers the class-independent form of its eq. (8), adopted there to secure
%    uniqueness of the fixed point; on OI stations that form loses the per-class
%    arrival-instant correction and was measured to be three to ten times less
%    accurate, so eq. (7) is used instead.
%
% The system, solved by successive substitution as in eq. (11) of the reference,
% is
%   W_{0,r} = Z_r
%   W_{j,r} = Dli(j,r) (1 + x_j - x_{j,r}/N_r)                     LI queue j
%   W_{i,r} = m_r / mu_{i,r}(m),  m = e_r + x_i(N - e_r)           OI station i
%   T_r     = N_r / (W_{0,r} + sum_j W_{j,r} + sum_i W_{i,r})
%   x_{k,r} = T_r W_{k,r}
% where x_k = sum_r x_{k,r}. A sweep costs O((J + K) R) residence times and
% O(K R^2 q) rate evaluations with q quadrature nodes, independently of the
% population. The number of sweeps is nearly so: 19 to 83 at the default
% tolerance from |N| = 9 to |N| = 45000 on a saturating three-class station.
% It does grow once an OI station saturates, because the fixed point is then
% neutrally stable -- with mu_i(n) -> sum_r a_r and the shares -> a_r/sum_s a_s,
% every x_i satisfies T_r = a_r -- and the trailing digits creep while the
% throughputs are already settled. That is why the tolerance is relative and
% defaults to 1e-6: tightening it to 1e-10 costs 16182 sweeps at |N| = 4500 and
% moves the throughputs by 3e-9, far below the 1e-2 error of the approximation
% itself. With every mu_i constant the OI equations reduce exactly to
% Bard-Schweitzer AMVA, since then mu_{i,r}(m) = mu_i m_r/|m| and
% W_{i,r} = |m|/mu_i.
%
% ACCURACY. Against PFQN_MVAOI on a battery of separable, shared, multiserver,
% mixed and multi-station OI models with two and three classes, the worst
% relative error on X and on the per-station queue-lengths was 8.6e-2, reached
% on a rate mixing a per-class and a shared term, with typical errors between
% 1e-3 and 6e-2. The residual is dominated by the Schweitzer step, which alone
% accounts for 1.5e-2 on a single-class queue-dependent station.
%
% SMOOTHNESS. Theorem 1 of the reference requires the queue-dependent factors to
% be differentiable on the real box, so MU IS ASSUMED TO BE SUPPLIED ALREADY
% SMOOTHED BY THE CALLER: each mu{i} must accept a real, non-negative occupancy
% vector and return a strictly positive rate on it. A rate written through the
% support indicator 1{n_r > 0}, as accepted by PFQN_NCOI and PFQN_MVAOI, is a
% step function of n and must be replaced by a smooth surrogate before being
% passed here, for instance by substituting n_r/(1 + n_r) or 1 - exp(-n_r) for
% the indicator. The same smoothed handle should then be passed to PFQN_MVAOI
% when the two are compared, since it defines the reference model.
%
% Parameters:
%   Z   - (1 x R) think-time demand vector of the aggregated delay node.
%   N   - (1 x R) closed population vector, finite.
%   mu  - cell array {mu_1,...,mu_K} of function handles; mu_i(n) returns the OI
%         total service rate of station i at the REAL occupancy vector n
%         (1 x R). A bare function handle is accepted as the single-station
%         shorthand. May be empty when there is no OI station.
%   Dli - (J x R) per-class demand matrix of the LI single-server queues
%         (D_{j,r} = V_{j,r}/rate_{j,r}); empty or omitted when J = 0.
%   options - solver options (optional). Recognized fields:
%         options.tol      successive-substitution tolerance (default 1e-6) on
%                          the RELATIVE increment max |x^{n+1} - x^n| / |x^{n+1}|
%                          over all (k,r). The criterion must be relative: the
%                          queue-lengths scale with the population, so a fixed
%                          absolute threshold demands more precision than the
%                          iteration can deliver at large N;
%         options.iter_max maximum number of sweeps (default 1000). Reaching it
%                          raises a warning naming the residual, never a silent
%                          unconverged return;
%         options.nodes    Gauss-Legendre nodes for the Aumann-Shapley integral
%                          (default 16);
%         options.dmu      cell array of gradient handles, dmu{i}(n) returning
%                          the (1 x R) gradient of mu_i at n. When omitted the
%                          gradient is taken by central differences;
%         options.init     (1+K+J) x R initial queue-length guess, rows ordered
%                          as [delay; OI stations; LI queues]. Defaults to the
%                          population spread evenly over the queueing stations.
%
% Returns:
%   X      - (1 x R) per-class throughput.
%   Qoi    - (K x R) per-class mean queue-length at each OI station (row i).
%   Qli    - (J x R) per-class mean queue-length at each LI queue (row j).
%   Qdelay - (1 x R) per-class mean queue-length at the delay node (X.*Z).
%   it     - number of successive-substitution sweeps performed.
%
% Example (smoothed two-class OI station, one LI queue):
%   sm   = @(v) v ./ (1 + v);                       % smoothed indicator
%   rate = @(n) 1.0*sm(n(1)) + 1.5*sm(n(2));        % smoothed OI rate
%   [X, Qoi, Qli] = pfqn_qdamva_oi([1 0.5], [8 6], {rate}, [0.3 0.2]);
%
% See also PFQN_MVAOI, PFQN_MVAOI_MARG, PFQN_NCOI, PFQN_OI_INSVC, PFQN_BS.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4
    Dli = [];
end
if nargin < 5
    options = struct();
end
if nargin < 3
    mu = {};
end
if isa(mu, 'function_handle')
    mu = {mu};
end
if isempty(mu)
    mu = {};
end
if ~iscell(mu)
    line_error(mfilename, 'mu must be a cell of OI rate handles.');
end
K = numel(mu);
for i = 1:K
    if ~isa(mu{i}, 'function_handle')
        line_error(mfilename, 'each mu{i} must be a function handle mu_i(n).');
    end
end

R = numel(N);
N = round(N(:)');
Z = Z(:)';
if numel(Z) ~= R
    line_error(mfilename, 'Z and N must have the same number of classes.');
end
if any(~isfinite(N))
    line_error(mfilename, 'pfqn_qdamva_oi requires finite (closed) populations.');
end
if isempty(Dli)
    Dli = zeros(0, R);
end
if size(Dli, 2) ~= R
    line_error(mfilename, 'Dli must have one column per class.');
end
J = size(Dli, 1);

tol = 1e-6;
if isfield(options, 'tol') && ~isempty(options.tol)
    tol = options.tol;
end
iter_max = 1000;
if isfield(options, 'iter_max') && ~isempty(options.iter_max)
    iter_max = options.iter_max;
end
nodes = 16;
if isfield(options, 'nodes') && ~isempty(options.nodes)
    nodes = round(options.nodes);
end
dmu = cell(1, K);
if isfield(options, 'dmu') && ~isempty(options.dmu)
    dmu = options.dmu;
    if isa(dmu, 'function_handle')
        dmu = {dmu};
    end
    if numel(dmu) ~= K
        line_error(mfilename, 'options.dmu must have one gradient handle per OI station.');
    end
end

X = zeros(1, R);
Qoi = zeros(K, R);
Qli = zeros(J, R);
Qdelay = zeros(1, R);
it = 0;
if sum(N) == 0
    return
end

% Schweitzer per-class ratios of eq. (7). A class with a single job sees the
% station empty of its own class on arrival, delta_r = 0.
deltar = zeros(1, R);
nz = N > 0;
deltar(nz) = (N(nz) - 1) ./ N(nz);

[qt, qw] = gauss_legendre01(nodes);

% Queue-length state: row 1 is the delay node, rows 2:K+1 the OI stations, rows
% K+2:K+J+1 the LI queues.
nq = K + J;
if isfield(options, 'init') && ~isempty(options.init)
    x = options.init;
    if ~isequal(size(x), [1 + nq, R])
        line_error(mfilename, 'options.init must be (1+K+J) x R.');
    end
else
    x = zeros(1 + nq, R);
    if nq > 0
        for r = 1:R
            x(2:end, r) = N(r) / nq;
        end
    else
        x(1, :) = N;
    end
end

W = zeros(1 + nq, R);
W(1, :) = Z;                       % the delay node is population-independent
for it = 1:iter_max
    xtot = sum(x, 2).';            % 1 x (1+nq) aggregate queue-lengths

    % OI residence times W_{i,r} = m_r / mu_{i,r}(m), the class rate being split
    % by the Aumann-Shapley shares of mu_i at the arrival-instant occupancy m.
    for i = 1:K
        xi = x(1 + i, :);
        for r = 1:R
            if N(r) == 0
                W(1 + i, r) = 0;
                continue
            end
            m = xi;
            m(r) = 1 + deltar(r) * xi(r);
            rate = mu{i}(m);
            if ~(rate > 0) || ~isfinite(rate)
                line_error(mfilename, sprintf(['mu{%d} returned a non-positive or non-finite rate at a ' ...
                    'real occupancy vector; it must be a smooth, strictly positive extension.'], i));
            end
            f = aumann_shapley(mu{i}, dmu{i}, m, R, qt, qw);
            W(1 + i, r) = m(r) / (rate * f(r));
        end
    end

    % LI residence times: the standard arrival-theorem term under eq. (7).
    for j = 1:J
        for r = 1:R
            if N(r) == 0
                W(1 + K + j, r) = 0;
                continue
            end
            W(1 + K + j, r) = Dli(j, r) * (1 + xtot(1 + K + j) - x(1 + K + j, r) / N(r));
        end
    end

    % Throughputs by population conservation, then the updated queue-lengths.
    xnew = zeros(1 + nq, R);
    for r = 1:R
        if N(r) == 0
            continue
        end
        den = sum(W(:, r));
        if ~(den > 0)
            line_error(mfilename, sprintf('class %d has no positive residence time at any station.', r));
        end
        X(r) = N(r) / den;
        xnew(:, r) = X(r) * W(:, r);
    end

    % Relative criterion. An absolute one is unusable here: the queue-lengths
    % scale with the population, so a fixed absolute threshold silently demands
    % more precision than the iteration can deliver at large N.
    err = max(abs(xnew(:) - x(:)) ./ max(abs(xnew(:)), eps));
    x = xnew;
    if err < tol
        break
    end
end
if err >= tol
    line_warning(mfilename, ['successive substitution stopped at iter_max = %d with relative residual ' ...
        '%g > tol = %g; the returned means are not converged.'], iter_max, err, tol);
end

Qdelay = x(1, :);
Qoi = x(2:1 + K, :);
Qli = x(2 + K:end, :);
end

% =========================================================================
function f = aumann_shapley(murate, dmurate, m, R, qt, qw)
% Aumann-Shapley shares of murate(m) among the R classes: the exact additive
% decomposition murate(m) - murate(0) = sum_r m_r int_0^1 d murate/d n_r (t m) dt.
% Exact for separable and for total-population rates, hence exact whenever the
% balanced-fairness split itself is available in closed form.
w = zeros(1, R);
for q = 1:numel(qt)
    p = qt(q) * m;
    if isempty(dmurate)
        for r = 1:R
            if m(r) <= 0
                continue
            end
            h = 1e-6 * max(1, abs(p(r)));
            pp = p; pp(r) = p(r) + h;
            pm = p; pm(r) = max(0, p(r) - h);
            w(r) = w(r) + qw(q) * (murate(pp) - murate(pm)) / (pp(r) - pm(r));
        end
    else
        g = dmurate(p);
        w = w + qw(q) * g(:).';
    end
end
w = max(m .* w, 0);
w(m <= 0) = 0;
s = sum(w);
if ~(s > 0)
    % Rate constant in n: the shared-rate split m_r/|m| is the exact answer.
    w = m;
    w(m <= 0) = 0;
    s = sum(w);
end
f = w / s;
end

% =========================================================================
function [t, w] = gauss_legendre01(n)
% Gauss-Legendre nodes and weights on [0,1] by the Golub-Welsch algorithm.
n = max(2, n);
k = 1:n-1;
b = k ./ sqrt(4 * k.^2 - 1);
T = diag(b, 1) + diag(b, -1);
[V, D] = eig(T);
[xs, ord] = sort(diag(D));
ws = 2 * V(1, ord).^2;
t = 0.5 * (xs(:).' + 1);
w = 0.5 * ws(:).';
end
