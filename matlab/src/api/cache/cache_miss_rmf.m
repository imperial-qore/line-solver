%{ @file cache_miss_rmf.m
 %  @brief Computes miss rates using the refined mean field (RMF) method
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes cache miss rates for RANDOM(m) replacement via refined mean field
 %
 % @details
 % This function computes global, per-user, and per-item miss rates for
 % multi-list caches with RANDOM(m) replacement using the DDPP mean field
 % approximation with 1/N correction (refined mean field). The cache
 % occupancy is computed from the aggregate request stream; per-user miss
 % rates weight the per-item miss probabilities by each user's rates.
 %
 % Reference:
 %   N. Gast, "Expected Values Estimated via Mean-Field Approximation are
 %   1/N-Accurate", Proc. ACM Meas. Anal. Comput. Syst., 2017.
 %
 % @par Syntax:
 % @code
 % [M, MU, MI, pi0] = cache_miss_rmf(gamma, m, lambda)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item access factors (used for sizing only)
 % <tr><td>m<td>Cache capacity vector
 % <tr><td>lambda<td>Arrival rates per user per item per list
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>M<td>Global miss rate
 % <tr><td>MU<td>Per-user miss rate
 % <tr><td>MI<td>Per-item miss rate
 % <tr><td>pi0<td>Per-item miss probability
 % </table>
%}
function [M,MU,MI,pi0,tout,pi0_t,MU_t,xtraj,xss,Wcov] = cache_miss_rmf(gamma, m, lambda, tspan, x0init, accost) %#ok<INUSL>
% XSS is the converged occupancy vector (model_dim x 1, item-major through
% RMF_INDEX) and WCOV the stationary covariance of the SAME process under the
% linear noise approximation, i.e. the solution of the Lyapunov equation
% F'(pi) W + W F'(pi)' + Q(pi) = 0 that RMF_EXPANSION_STEADY_STATE already
% solves on the way to the 1/N mean correction. It was discarded there; it is
% the second moment of the cache occupancy and is what SolverFLD's moment
% closure reports for a cache node. WCOV is empty when the refinement did not
% run (a non-linear access graph has no 1/N path) or when its reduced linear
% system was singular.
%
% Optional TSPAN = [t0, t1] and initial occupancy X0INIT request the
% transient mean-field trajectory: integrating the same drift F(x) that
% RMF_FIXED_POINT drives to steady state, over the finite window. Returns
% the time grid TOUT, the per-item list-0 occupancy PI0_T (n_items x nt)
% and the per-user miss-rate trajectory MU_T (u x nt). Omitting TSPAN
% preserves the original steady-state-only contract.
%
% Optional ACCOST is the per-(user,item) access graph (each an (h+1)x(h+1)
% matrix): row 1 is miss admission (col 1 = reject, col 1+l = admit to list l),
% row 1+i is a hit in list i (col 1+b = promote to list b>=i). A non-linear
% graph modulates admission/promotion per item; the drift then follows the
% general RANDOM(m) dynamics (RMF_DRIFT_GRAPH) with the plain mean-field fixed
% point, while the standard linear chain keeps the 1/N-refined path unchanged.
if nargin < 4, tspan = []; end
if nargin < 5, x0init = []; end
if nargin < 6, accost = []; end
u = size(lambda,1);
n_items = size(lambda,2);
h = length(m);
m = m(:)';

% aggregate per-item request rates over users
lam_i = zeros(1, n_items);
for v = 1:u
    row = lambda(v,:,1);
    row(~isfinite(row)) = 0;
    lam_i = lam_i + row;
end
p = lam_i / sum(lam_i);

model_dim = n_items * (h + 1);

% Build initial state: first m(1) items in list 1, next m(2) in
% list 2, etc.; remaining items outside cache (list 0)
x0 = zeros(model_dim, 1);
obj_idx = 0;
for k = 1:h
    for jj = 1:m(k)
        obj_idx = obj_idx + 1;
        if obj_idx <= n_items
            x0(rmf_index(obj_idx, k, n_items)) = 1.0;
        end
    end
end
for i = (obj_idx + 1):n_items
    x0(rmf_index(i, 0, n_items)) = 1.0;
end

% see _kb/09-ldes-and-cache.md (mean-field cache fixed points and TTL edge cases)
% A non-default access graph (accost) modulates admission (row 1) and promotion
% (row 1+i) per item; the general drift honours it with the plain mean-field
% fixed point, while the linear chain keeps the 1/N-refined path.
G = rmf_build_item_graphs(accost, lambda, n_items, h);
Wcov = [];
if ~isempty(G)
    xss = rmf_fixed_point_graph(x0, p, G, m, n_items, h, model_dim);
else
    xss = rmf_fixed_point(x0, p, m, n_items, h, model_dim);
    % The refinement's reduced linear system is expected to be singular for
    % small/non-hyperbolic fixed points; silence only that specific warning here
    % since the non-finite result is detected and rejected below.
    ws1 = warning('off', 'MATLAB:singularMatrix');
    ws2 = warning('off', 'MATLAB:nearlySingularMatrix');
    restoreWarn = onCleanup(@() warning([ws1 ws2]));
    try
        [pi_mf, V] = rmf_expansion_steady_state(x0, p, m, n_items, h, model_dim);
        xref = pi_mf + V / n_items;
        if all(isfinite(xref))
            xss = xref;
        end
    catch
        % Fall back to plain mean field
    end
    try
        W = rmf_lna_covariance(xss, p, m, n_items, h, model_dim);
        if all(all(isfinite(W)))
            Wcov = W;
        end
    catch
        % no stationary covariance at this fixed point; the mean stands
    end
    clear restoreWarn
end

% Per-item miss probability (occupancy of list 0), clipped to [0,1]
pi0 = zeros(n_items, 1);
for i = 1:n_items
    pi0(i) = max(0, min(1, xss(rmf_index(i, 0, n_items))));
end

MI = lam_i(:) .* pi0;
MU = zeros(1, u);
for v = 1:u
    row = lambda(v,:,1);
    row(~isfinite(row)) = 0;
    MU(v) = row * pi0;
end
M = sum(MI);

% Transient mean-field trajectory (optional): integrate the drift over the
% requested window from the supplied (or default) initial occupancy. XTRAJ is
% the full DDPP occupancy trajectory (model_dim x nt), used to carry the cache
% mean occupancy across environment switches.
tout = []; pi0_t = []; MU_t = []; xtraj = [];
if ~isempty(tspan)
    if isempty(x0init)
        x0init = x0;
    else
        x0init = x0init(:);
    end
    ode_func = @(t, x) rmf_drift(x, p, m, n_items, h, model_dim);
    odeopt = odeset('AbsTol', 1e-10, 'RelTol', 1e-8);
    [tout, xtraj] = ode15s(ode_func, tspan, x0init, odeopt);
    xtraj = xtraj'; % model_dim x nt
    nt = numel(tout);
    pi0_t = zeros(n_items, nt);
    for i = 1:n_items
        pi0_t(i,:) = max(0, min(1, xtraj(rmf_index(i, 0, n_items), :)));
    end
    MU_t = zeros(u, nt);
    for v = 1:u
        row = lambda(v,:,1);
        row(~isfinite(row)) = 0;
        MU_t(v,:) = row * pi0_t;
    end
end
end

function idx = rmf_index(i, k, n_items)
% RMF_INDEX Map (item i, list k) to flat state index
%   i: item index (1-based)
%   k: list index (0 = outside cache, 1..h = cache lists)
%   n_items: total number of items
idx = i + k * n_items;
end

function hr = rmf_hit_rate(x, p, list_number, n_items)
% RMF_HIT_RATE Compute hit rate contribution from a specific list
%   hr = sum_i p(i) * x(index(i, list_number))
hr = 0.0;
for i = 1:n_items
    hr = hr + p(i) * x(rmf_index(i, list_number, n_items));
end
end

function dX = rmf_drift(x, p, m, n_items, h, model_dim)
% RMF_DRIFT Compute mean field drift F(x) for RANDOM(m) replacement
%
% dx[i,k]/dt = -p(i)*x[i,k] + hitRate[k]*x[i,k+1]/m(k)  (promotion from k to k+1)
% dx[i,k+1]/dt = p(i)*x[i,k] - hitRate[k]*x[i,k+1]/m(k)

hit_rates = zeros(1, h + 1);
for k = 0:h
    hit_rates(k + 1) = rmf_hit_rate(x, p, k, n_items);
end

dX = zeros(model_dim, 1);
for i = 1:n_items
    for k = 0:(h - 1)
        flow = p(i) * x(rmf_index(i, k, n_items)) ...
             - hit_rates(k + 1) * x(rmf_index(i, k + 1, n_items)) / m(k + 1);
        dX(rmf_index(i, k, n_items))     = dX(rmf_index(i, k, n_items))     - flow;
        dX(rmf_index(i, k + 1, n_items)) = dX(rmf_index(i, k + 1, n_items)) + flow;
    end
end
end

function Fp = rmf_jacobian(x, p, m, n_items, h, model_dim)
% RMF_JACOBIAN Compute Jacobian dF/dx at state x

hit_rates = zeros(1, h + 1);
for k = 0:h
    hit_rates(k + 1) = rmf_hit_rate(x, p, k, n_items);
end

Fp = zeros(model_dim, model_dim);
for i = 1:n_items
    for k = 0:(h - 1)
        ik  = rmf_index(i, k, n_items);
        ik1 = rmf_index(i, k + 1, n_items);

        % Direct rate terms
        Fp(ik,  ik)  = Fp(ik,  ik)  - p(i);
        Fp(ik1, ik)  = Fp(ik1, ik)  + p(i);
        Fp(ik,  ik1) = Fp(ik,  ik1) + hit_rates(k + 1) / m(k + 1);
        Fp(ik1, ik1) = Fp(ik1, ik1) - hit_rates(k + 1) / m(k + 1);

        % Indirect terms via hit rate dependence on x(j,k)
        for j = 1:n_items
            jk  = rmf_index(j, k, n_items);
            jk1 = rmf_index(j, k + 1, n_items);
            Fp(ik,  jk1) = Fp(ik,  jk1) - p(i) * x(ik) / m(k + 1);
            Fp(ik1, jk1) = Fp(ik1, jk1) + p(i) * x(ik) / m(k + 1);
            Fp(ik,  jk)  = Fp(ik,  jk)  + p(j) * x(ik1) / m(k + 1);
            Fp(ik1, jk)  = Fp(ik1, jk)  - p(j) * x(ik1) / m(k + 1);
        end
    end
end
end

function Fpp = rmf_hessian(~, p, m, n_items, h, model_dim)
% RMF_HESSIAN Compute Hessian d^2F/dx^2 (constant for this quadratic drift)

Fpp = zeros(model_dim, model_dim, model_dim);
for i = 1:n_items
    for k = 0:(h - 1)
        ik  = rmf_index(i, k, n_items);
        ik1 = rmf_index(i, k + 1, n_items);
        for j = 1:n_items
            if j ~= i
                jk  = rmf_index(j, k, n_items);
                jk1 = rmf_index(j, k + 1, n_items);
                % d^2 F[ik] / (d x[jk] d x[ik1])
                Fpp(ik,  jk,  ik1) = Fpp(ik,  jk,  ik1) + p(j) / m(k + 1);
                Fpp(ik,  ik1, jk)  = Fpp(ik,  ik1, jk)  + p(j) / m(k + 1);
                % d^2 F[ik] / (d x[jk1] d x[ik])
                Fpp(ik,  jk1, ik)  = Fpp(ik,  jk1, ik)  - p(i) / m(k + 1);
                Fpp(ik,  ik,  jk1) = Fpp(ik,  ik,  jk1) - p(i) / m(k + 1);
                % Symmetric for ik1
                Fpp(ik1, jk,  ik1) = Fpp(ik1, jk,  ik1) - p(j) / m(k + 1);
                Fpp(ik1, ik1, jk)  = Fpp(ik1, ik1, jk)  - p(j) / m(k + 1);
                Fpp(ik1, jk1, ik)  = Fpp(ik1, jk1, ik)  + p(i) / m(k + 1);
                Fpp(ik1, ik,  jk1) = Fpp(ik1, ik,  jk1) + p(i) / m(k + 1);
            end
        end
    end
end
end

function Q = rmf_noise_matrix(x, p, m, n_items, h, model_dim)
% RMF_NOISE_MATRIX Compute noise intensity matrix Q(x) for the DDPP
%
% Q[a,b] = sum_ell ell[a]*ell[b]*beta_ell(x)
% Each transition swaps items i and j across lists k and k+1.

Q = zeros(model_dim, model_dim);
signs = [-1, 1, 1, -1];
for i = 1:n_items
    for k = 0:(h - 1)
        for j = 1:n_items
            rate = p(i) * x(rmf_index(i, k, n_items)) ...
                 * x(rmf_index(j, k + 1, n_items)) / m(k + 1);
            indices = [rmf_index(i, k, n_items), ...
                       rmf_index(j, k, n_items), ...
                       rmf_index(i, k + 1, n_items), ...
                       rmf_index(j, k + 1, n_items)];
            for ia = 1:4
                for ib = 1:4
                    Q(indices(ia), indices(ib)) = Q(indices(ia), indices(ib)) ...
                        + rate * signs(ia) * signs(ib);
                end
            end
        end
    end
end
end

function pi = rmf_fixed_point(x0, p, m, n_items, h, model_dim)
% RMF_FIXED_POINT Compute mean field fixed point by ODE integration
%
% Integrates dx/dt = F(x) until steady state using ode15s.

tmax = 10000;
ode_func = @(t, x) rmf_drift(x, p, m, n_items, h, model_dim);
odeopt = odeset('AbsTol', 1e-10, 'RelTol', 1e-8);
[~, xvec] = ode15s(ode_func, [0, tmax], x0, odeopt);
%[~, xvec] = lsoda_solve(ode_func, [0, tmax], x0, odeopt);
pi = xvec(end, :)';
end

function [C, Cinv, rk] = rmf_dimension_reduction(Fp, n_items, h, model_dim)
% RMF_DIMENSION_REDUCTION Compute change-of-basis for singular Jacobian
%
% The Jacobian is singular because item populations are conserved
% (sum over lists for each item = 1). Returns matrices to project
% onto the non-singular subspace.

rk = rank(Fp);

C = zeros(model_dim, model_dim);
d = 0;
for l_idx = 0:h
    for i = 1:(n_items - 1)
        d = d + 1;
        C(d, rmf_index(i, l_idx, n_items)) = 1.0;
    end
end

[U, ~, ~] = svd(Fp);
C((rk + 1):model_dim, :) = U(:, (rk + 1):model_dim)';
Cinv = inv(C);
end

function W = rmf_lna_covariance(x, p, m, n_items, h, model_dim)
% W = RMF_LNA_COVARIANCE(X, P, M, N_ITEMS, H, MODEL_DIM)
%
% Stationary covariance of the RANDOM(m) occupancy process under the linear
% noise approximation: the solution of F'(x) W + W F'(x)' + Q(x) = 0 with the
% same Jacobian RMF_JACOBIAN and the same noise matrix RMF_NOISE_MATRIX that
% the 1/N correction is built from, so the mean and the covariance linearise
% about the identical drift.
%
% THE SUBSPACE IS THE POINT. The Jacobian is singular twice over, because the
% cache conserves two things: each item is in exactly one list
% (sum_k x(i,k) = 1) and each list holds exactly its capacity
% (sum_i x(i,k) = m(k)). The fluctuation therefore lives on the span of the
% jump directions, and every jump is a SWAP,
%   l(i,j,k) = (e_i - e_j) tensor (e_{k+1} - e_k),
% so that span is the tensor product of the zero-sum item space with the
% zero-sum list space -- the double-centred subspace, of dimension
% (n_items-1)*h. Restricting to an orthonormal basis of it is an EXACT
% reduction, and it is what makes the answer respect both conservation laws:
% the covariance of a deterministic total must be zero, and on n=5, m=2 the
% reduction used by RMF_EXPANSION_STEADY_STATE (which drops the last item's
% rows and pads with SVD null vectors) returns sum(W(:)) = +0.47 for the miss
% indicators where the exact chain gives 0.
%
% Parameters:
%   x         - occupancy at which to linearise (model_dim x 1)
%   p         - per-item request probabilities (1 x n_items)
%   m         - list capacities (1 x h)
%   n_items   - number of items
%   h         - number of lists
%   model_dim - n_items*(h+1)
%
% Returns:
%   W - model_dim x model_dim covariance, item-major through RMF_INDEX
%
% See also CACHE_MISS_RMF, RMF_NOISE_MATRIX, FLUID_LYAPUNOV.

Fp = rmf_jacobian(x, p, m, n_items, h, model_dim);
Q = rmf_noise_matrix(x, p, m, n_items, h, model_dim);
Q = (Q + Q') / 2;

Ui = local_centered_basis(n_items);   % n_items x (n_items-1)
Ul = local_centered_basis(h + 1);     % (h+1) x h
V = kron(Ul, Ui);                     % item-major flat index i + k*n_items
if isempty(V)
    W = zeros(model_dim);
    return
end

Ar = V' * Fp * V;
Qr = V' * Q * V;
Qr = (Qr + Qr') / 2;

% the linear noise approximation has a stationary covariance only at an
% exponentially stable fixed point
if max(real(eig(Ar))) >= -sqrt(eps)
    line_error(mfilename, ['The cache fluid fixed point is not exponentially stable on the ' ...
        'reachable subspace, so the occupancy process has no stationary covariance.']);
end

Wr = lyap(Ar, Qr);
Wr = (Wr + Wr') / 2;
W = V * Wr * V';
W = (W + W') / 2;
end

function U = local_centered_basis(n)
% Orthonormal basis of {u in R^n : sum(u) = 0}, n-by-(n-1).
if n <= 1
    U = zeros(n, 0);
    return
end
[U, ~] = qr(eye(n) - ones(n)/n, 0);
U = U(:, 1:(n-1));
end

function [pi, V] = rmf_expansion_steady_state(x0, p, m, n_items, h, model_dim)
% RMF_EXPANSION_STEADY_STATE Compute refined mean field steady-state expansion
%
% Computes the mean field fixed point pi and the 1/N correction V
% using the Lyapunov equation approach with dimension reduction.
%
% The refined approximation for a system of N items is:
%   E[X] ~ pi + V/N + O(1/N^2)

pi = rmf_fixed_point(x0, p, m, n_items, h, model_dim);

Fp  = rmf_jacobian(pi, p, m, n_items, h, model_dim);
Fpp = rmf_hessian(pi, p, m, n_items, h, model_dim);
Q   = rmf_noise_matrix(pi, p, m, n_items, h, model_dim);

% Dimension reduction: project onto non-singular subspace
[C, Cinv, rk] = rmf_dimension_reduction(Fp, n_items, h, model_dim);

Fp_r = (C * Fp * Cinv);
Fp_r = Fp_r(1:rk, 1:rk);

% Reduce Hessian: Fpp_r(a,b,c) = sum_{i,j,k} C(a,i)*Fpp(i,j,k)*Cinv(j,b)*Cinv(k,c)
% First contraction: tmp1(a,j,k) = sum_i C(a,i)*Fpp(i,j,k)
tmp1 = zeros(model_dim, model_dim, model_dim);
for a = 1:rk
    for j = 1:model_dim
        for k = 1:model_dim
            tmp1(a, j, k) = C(a, :) * Fpp(:, j, k);
        end
    end
end
% Second contraction: tmp2(a,b,k) = sum_j tmp1(a,j,k)*Cinv(j,b)
tmp2 = zeros(rk, rk, model_dim);
for a = 1:rk
    for b = 1:rk
        for k = 1:model_dim
            tmp2(a, b, k) = tmp1(a, :, k) * Cinv(:, b);
        end
    end
end
% Third contraction: Fpp_r(a,b,c) = sum_k tmp2(a,b,k)*Cinv(k,c)
Fpp_r = zeros(rk, rk, rk);
for a = 1:rk
    for b = 1:rk
        Fpp_r(a, b, :) = reshape(tmp2(a, b, :), 1, []) * Cinv(:, 1:rk);
    end
end

Q_r = C * Q * C';
Q_r = Q_r(1:rk, 1:rk);

% Solve Lyapunov equation: Fp_r * W_r + W_r * Fp_r' + Q_r = 0
W_r = lyap(Fp_r, Q_r);

% First-order correction: V_r = -Fp_r \ (C_r / 2)
% where C_r = sum_{b,c} Fpp_r(:,b,c) * W_r(b,c)
C_r = zeros(rk, 1);
for a = 1:rk
    for b = 1:rk
        for c = 1:rk
            C_r(a) = C_r(a) + Fpp_r(a, b, c) * W_r(b, c);
        end
    end
end
V_r = -Fp_r \ (C_r / 2.0);

% Expand back to full dimension. The dropped coordinates span the null
% directions of the Jacobian, i.e. the per-item conservation sum_k x(i,k) = 1,
% which carries no fluctuation, so the expansion of W is exact rather than a
% truncation.
V = Cinv(:, 1:rk) * V_r;
W = Cinv(:, 1:rk) * W_r * Cinv(:, 1:rk)';
W = (W + W') / 2;
end

function g = rmf_linear_graph(h)
% RMF_LINEAR_GRAPH Standard linear chain: miss->list1, hit in list a->list a+1,
% self-loop on the top list. (h+1)x(h+1), 1-based (col/row 1 = out/reject).
g = zeros(h+1, h+1);
g(1, 2) = 1;
for a = 1:(h-1)
    g(a+1, a+2) = 1;
end
g(h+1, h+1) = 1;
end

function G = rmf_build_item_graphs(accost, lambda, n, h)
% RMF_BUILD_ITEM_GRAPHS Per-item (h+1)x(h+1) access graph aggregated over users
% by request rate. Returns {} when accost is absent or the standard linear
% chain (so the caller keeps the refined linear path). ACCOST is cell{v,k}.
G = {};
if isempty(accost)
    return;
end
lin = rmf_linear_graph(h);
u = size(accost, 1);
Gc = cell(1, n);
isLinear = true;
for k = 1:n
    num = zeros(h+1, h+1); den = 0;
    for v = 1:u
        wv = sum(lambda(v, k, 1));
        if ~isfinite(wv), wv = 0; end
        gvk = accost{v, k};
        if isempty(gvk), continue; end
        num = num + wv * gvk;
        den = den + wv;
    end
    if den > 0
        gk = num / den;
    else
        gk = accost{1, k};
        if isempty(gk), gk = lin; end
    end
    for a = 1:(h+1)
        srow = sum(gk(a, :));
        if srow > 0, gk(a, :) = gk(a, :) / srow; end
    end
    Gc{k} = gk;
    if ~all(all(abs(gk - lin) < 1e-9))
        isLinear = false;
    end
end
if ~isLinear
    G = Gc;
end
end

function dX = rmf_drift_graph(x, p, G, m, n, h, model_dim)
% RMF_DRIFT_GRAPH General RANDOM(m) mean-field drift honouring per-item access
% graph G{k} (h+1)x(h+1). Miss admission weighted by row 1, hit promotion by
% row 1+i; a uniformly random occupant of the target list is displaced (evicted
% out on a miss, swapped to the source list on a hit), matching the exact RR
% sample path (State.afterEventCache). Reduces to the linear eq-8 drift when G
% is the standard chain.
x = max(0, min(1, x));
A = zeros(h+1, h+1);   % A(s+1, i): insertion/promotion into list i from source s
for s = 0:h
    for j = 1:n
        xjs = x(rmf_index(j, s, n));
        if xjs == 0, continue; end
        gj = G{j};
        for i = 1:h
            A(s+1, i) = A(s+1, i) + p(j) * xjs * gj(s+1, i+1);
        end
    end
end
dX = zeros(model_dim, 1);
for k = 1:n
    outk = x(rmf_index(k, 0, n));
    gk = G{k};
    for i = 1:h
        xki = x(rmf_index(k, i, n));
        infl = p(k) * outk * gk(1, i+1);
        for s = 1:(i-1)
            infl = infl + p(k) * x(rmf_index(k, s, n)) * gk(s+1, i+1);
        end
        for b = (i+1):h
            infl = infl + A(i+1, b) * x(rmf_index(k, b, n)) / m(b);
        end
        outfl = p(k) * xki * (1 - gk(i+1, i+1));
        disp = 0;
        for s = 0:(i-1)
            disp = disp + A(s+1, i);
        end
        outfl = outfl + disp * xki / m(i);
        dX(rmf_index(k, i, n)) = dX(rmf_index(k, i, n)) + infl - outfl;
    end
    acc = 0;
    for i = 1:h
        acc = acc + dX(rmf_index(k, i, n));
    end
    dX(rmf_index(k, 0, n)) = -acc;
end
end

function pj = rmf_fixed_point_graph(x0, p, G, m, n, h, model_dim)
% RMF_FIXED_POINT_GRAPH Plain mean-field fixed point of the general drift.
tmax = 20000;
ode_func = @(t, x) rmf_drift_graph(x, p, G, m, n, h, model_dim);
odeopt = odeset('AbsTol', 1e-10, 'RelTol', 1e-8);
[~, xvec] = ode15s(ode_func, [0, tmax], x0, odeopt);
pj = xvec(end, :)';
end
