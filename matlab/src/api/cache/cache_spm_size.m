%{ @file cache_spm_size.m
 %  @brief Ray (WKB) asymptotic expansion of the cost-capped cache normalizing constant
 %
 %  @author LINE Development Team
%}

%{
 % @brief Approximates the cost-constrained cache normalizing constant E(m,k) by the ray expansion
 %
 % @details
 % Evaluates the geometrical-optics (WKB) asymptotic expansion of the cost-capped
 % normalizing constant computed exactly by cache_erec(gamma,m,sigma,k), and
 % returns it in the SAME normalization, so that the two are interchangeable.
 % This is the item-size extension of retrieval_rayint, which carries the
 % size-free expansion; call that one when there are no storage costs.
 %
 % Writing E = prod_j m_j! * H, the size-free recursion
 %   E(m,n) = E(m,n-1) + sum_j gamma_{n,j} m_j E(m-1_j,n-1),  E(0,0)=1
 % relaxes to H ~ exp(phi/eps) with n = y/eps, m_j = x_j/eps, whose eikonal
 % e^{phi_y} = 1 + sum_j gamma_j(y) e^{-phi_j} carries the ray constants
 % xi_j = e^{-phi_j}. With per-item storage costs sigma_i and per-list cost caps
 % k_j the recursion gains the cost coordinate,
 %   E(m,k) = E_i(m,k) + sum_j m_j gamma_ij E_i(m-1_j, k-sigma_i 1_j),
 % so the shift 1_j becomes e_j(y) = (1_j, s(y) 1_j) in the enlarged space
 % X = (x,kappa) and the eikonal picks up the size tilt
 %   e^{phi_y} = 1 + sum_j gamma_j(y) e^{-phi_{x_j} - s(y) phi_{kappa_j}},
 % with the second family of ray constants zeta_j = e^{-phi_{kappa_j}}. The rays
 % integrate to the discrete saddle point of the product generating function
 %   sum_{m,k} H(m,k) prod_j z_j^{m_j} w_j^{k_j}
 %       = prod_i ( 1 + sum_j gamma_ij z_j w_j^{sigma_i} ),
 % namely, with D_i = 1 + sum_j gamma_ij xi_j zeta_j^{sigma_i} and
 % Psi = sum_i log D_i,
 %   m_j = sum_i gamma_ij xi_j zeta_j^{sigma_i} / D_i,
 %   k_j = sum_i sigma_i gamma_ij xi_j zeta_j^{sigma_i} / D_i,
 %   log H(m,k) ~ Psi - sum_j m_j log xi_j - sum_j k_j log zeta_j
 %                - (d/2) log(2 pi) - (1/2) log det grad^2 Psi,
 % where d is the number of saddle coordinates and, with
 % pi_ij = gamma_ij xi_j zeta_j^{sigma_i} / D_i and
 % Q^i_{jl} = delta_{jl} pi_ij - pi_ij pi_il,
 %   grad^2 Psi = sum_i [1; sigma_i] [1; sigma_i]' (x) Q^i .
 % Setting zeta_j = 1 recovers the size-free expansion exactly.
 %
 % CAPS ARE CUMULATIVE. cache_erec sums over the states of cost AT MOST k_j, so
 % this function does the same by default ('atmost'). The shadow price
 % eta_j = log zeta_j <= 0 obeys complementary slackness: a list whose
 % unconstrained mean cost already meets its cap is SLACK, keeps zeta_j = 1, and
 % drops out of the saddle, which then degenerates continuously to the size-free
 % expansion; a list whose cap BINDS sits at eta_j < 0, and the states below the
 % boundary decay geometrically with ratio zeta_j, contributing the amplitude
 % factor 1/(1-zeta_j). Pass 'exact' to obtain instead the constant resolving the
 % cost exactly at k_j, which is the raw Laplace formula above with no such factor.
 %
 % SIZE DIVERSITY IS REQUIRED. The Hessian integrand [1;sigma_i][1;sigma_i]' (x) Q^i
 % has rank h, not 2h, so grad^2 Psi is nonsingular only if the sizes actually vary.
 % This is not an artefact: with a single item size the cost of list j is
 % sigma*m_j identically and the cap carries no information. That case is detected
 % and answered exactly (slack, or zero when the cap cannot be met) rather than
 % passed to a singular saddle. If the sizes share a common divisor the cost lives
 % on a sublattice; the sizes and caps are divided through by their gcd, which is
 % an exact reduction and removes the corresponding lattice factor.
 %
 % OCCUPANCY. out.pij is the saddle occupancy
 % pi_il = gamma_il zeta_l^{sigma_i} xi_l / D_i and out.K its per-list cost. These
 % are EXACT-COST quantities: the saddle conditions are sum_i pi_ij = m_j and
 % sum_i sigma_i pi_ij = k_j, so out.K equals the cap exactly on every binding
 % list. Under cumulative caps the true mean cost is strictly below the cap; use
 % cache_cost and cache_prob_erec for that.
 %
 % ACCURACY. The expansion is O(1/n) at fixed occupancy. With a well-separated cap
 % the observed error in log E is around 1e-2 at n = 200 and halves at each
 % doubling of n. It degrades as a binding zeta_j approaches 1, i.e. in the
 % transition between the binding and slack regimes, where the geometric
 % resummation 1/(1-zeta_j) is no longer sharp; out.zeta and out.binding report
 % where the saddle sits and a warning is raised inside that region.
 %
 % @par Syntax:
 % @code
 % E = cache_spm_size(gamma, m, sigma, k)
 % E = cache_spm_size(gamma, m, sigma, k, costmode)
 % [E, logE, out] = cache_spm_size(...)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities (access factors), n x h
 % <tr><td>m<td>Cache list capacities, 1 x h
 % <tr><td>sigma<td>Item storage costs (sizes), 1 x n, positive integers
 % <tr><td>k<td>Per-list storage cost caps, 1 x h, non-negative integers
 % <tr><td>costmode<td>(Optional) 'atmost' (default, matches cache_erec) or 'exact'
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>E<td>Normalizing constant, same normalization as cache_erec (Inf on overflow; use logE)
 % <tr><td>logE<td>Natural logarithm of E, safe for large n
 % <tr><td>out<td>Struct with xi, zeta, binding, pij, K, phi, logdetSigma, span, method, relerrEst, iter
 % </table>
 %
 % @par References:
 % G. Casale, N. Gast, "Performance Analysis Methods for List-Based Caches With
 % Non-Uniform Access", IEEE/ACM Trans. Networking 29(2), 2021 (product form and
 % the cost-capped recursion, Sec. IX); G. Casale, "Accelerating Performance
 % Inference over Closed Systems by Asymptotic Methods", ACM SIGMETRICS, 2017
 % (the asymptotic expansion of the normalizing constant integral).
 %
 % @see cache_erec, cache_cost, cache_prob_erec, retrieval_rayint, cache_spm
%}
function [E, logE, out] = cache_spm_size(gamma, m, sigma, k, costmode)

if nargin < 4
    line_error(mfilename, ['Four arguments are required: the access factors, the list ' ...
        'capacities, the item sizes and the cost caps. Use retrieval_rayint for the ' ...
        'size-free expansion.']);
end
if ~ismatrix(gamma) || isempty(gamma)
    line_error(mfilename, 'The access factors must be a non-empty n x h matrix.');
end
[n0, h0] = size(gamma);
m = m(:).';
if numel(m) ~= h0
    line_error(mfilename, 'The capacity vector must have one entry per cache list (%d given, %d expected).', numel(m), h0);
end
if any(m < 0) || any(abs(m - round(m)) > 0)
    line_error(mfilename, 'List capacities must be non-negative integers.');
end
if isempty(sigma) || isempty(k)
    line_error(mfilename, ['The item sizes and the cost caps are both required. Use ' ...
        'retrieval_rayint for the size-free expansion.']);
else
    sigma = sigma(:).'; k = k(:).';
    if numel(sigma) ~= n0
        line_error(mfilename, 'The item size vector must have one entry per item.');
    end
    if numel(k) ~= h0
        line_error(mfilename, 'The cost cap vector must have one entry per cache list.');
    end
    if any(sigma <= 0) || any(abs(sigma - round(sigma)) > 0)
        line_error(mfilename, 'Item sizes must be positive integers.');
    end
    if any(abs(k - round(k)) > 0)
        line_error(mfilename, 'Storage cost caps must be integers.');
    end
end
capped = true;      % cleared below when a single item size makes the cap uninformative
if nargin < 5 || isempty(costmode)
    costmode = 'atmost';
end
costmode = lower(costmode);
if ~any(strcmp(costmode, {'atmost','exact'}))
    line_error(mfilename, 'The cost mode must be ''atmost'' or ''exact'' (''%s'' given).', costmode);
end

out = struct('xi', zeros(1,h0), 'zeta', ones(1,h0), 'binding', false(1,h0), ...
             'pij', [], 'K', zeros(1,h0), 'phi', NaN, 'logdetSigma', NaN, ...
             'span', 1, 'method', '', 'relerrEst', NaN, 'iter', 0);
out.pij = [ones(n0,1), zeros(n0,h0)];

% --- boundaries, matching cache_erec ---
if sum(m) > n0 || (capped && any(k < 0))
    E = 0; logE = -Inf; out.method = 'boundary'; return
end
if sum(m) == 0
    E = 1; logE = 0; out.method = 'boundary';
    if capped && strcmp(costmode,'exact') && any(k > 0)
        E = 0; logE = -Inf;
    end
    return
end

% --- items that can never be cached and lists of zero capacity drop out ---
alive = sum(gamma, 2) > 0;
live  = m > 0;
G = gamma(alive, live);
mk = m(live);
hk = numel(mk);
n  = size(G, 1);
sg = sigma(alive);
kk = k(live);
if sum(mk) > n
    E = 0; logE = -Inf; out.method = 'boundary'; return
end
if sum(mk) == n
    line_error(mfilename, ['The expansion requires sum(m) < n: at sum(m) = n the saddle point ' ...
        'escapes to infinity. Use cache_erec for a full cache.']);
end

% --- exact reductions on the cost lattice ---
if capped
    span = sg(1);
    for i = 2:n
        span = gcd(span, sg(i));
    end
    out.span = span;
    if strcmp(costmode, 'exact') && any(mod(kk, span) ~= 0)
        E = 0; logE = -Inf; out.method = 'lattice'; return   % unreachable off the sublattice
    end
    sg = sg / span;
    kk = floor(kk / span);
    % per-list feasibility: the m_j cheapest (dearest) reachable items bound the cost
    for j = 1:hk
        idx = find(G(:,j) > 0);
        if numel(idx) < mk(j)
            E = 0; logE = -Inf; out.method = 'boundary'; return
        end
        srt = sort(sg(idx));
        if sum(srt(1:mk(j))) > kk(j)
            E = 0; logE = -Inf; out.method = 'boundary'; return
        end
        if strcmp(costmode, 'exact') && sum(srt(end-mk(j)+1:end)) < kk(j)
            E = 0; logE = -Inf; out.method = 'boundary'; return
        end
    end
    % a single item size makes the cost of list j equal to sigma*m_j identically,
    % so the cap carries no information and the 2h saddle is singular (rank h)
    if all(sg == sg(1))
        cost = sg(1) * mk;
        if strcmp(costmode, 'exact')
            feasible = all(cost == kk);
        else
            feasible = all(cost <= kk);
        end
        if ~feasible
            E = 0; logE = -Inf; out.method = 'uniform-size'; return
        end
        capped = false;                    % fall through to the size-free expansion
        out.method = 'uniform-size';
    end
end

% --- saddle point ---
if capped
    [th, et, bind, it, P, D] = local_saddle_cost(G, mk, sg, kk, strcmp(costmode,'atmost'));
    % reshape, because at h = 1 find() returns 0x0 rather than a 1x0 row, and the
    % empty matrix product (0x1)*(1x0) is EMPTY rather than zero -- which silently
    % turns phi into [] instead of erroring. Sum elementwise for the same reason.
    ix  = reshape(find(bind), 1, []);
    dof = hk + numel(ix);
    phi = sum(log(D)) - sum(mk.*th) - sum(kk(ix).*et(ix));
    Sig = local_hessian_cost(P, sg, ix);
    logdet = local_logdet(Sig);
    logH = phi - (dof/2)*log(2*pi) - 0.5*logdet;
    if strcmp(costmode, 'atmost') && ~isempty(ix)
        logH = logH - sum(log1p(-exp(et(ix))));   % geometric resummation below the cap
    end
    if isempty(out.method)
        out.method = 'spm-size';
    end
else
    [th, it, P, D] = local_saddle(G, mk);
    et   = zeros(1, hk);
    bind = false(1, hk);
    phi  = sum(log(D)) - sum(mk.*th);
    Sig  = local_hessian_cost(P, zeros(1,n), []);
    logdet = local_logdet(Sig);
    logH = phi - (hk/2)*log(2*pi) - 0.5*logdet;
    if isempty(out.method)
        out.method = 'spm';
    end
end

logE = logH + sum(gammaln(m + 1));       % back to the cache_erec normalization
E    = exp(logE);

% --- ray quantities, reported on the original item and list indexing ---
xi = zeros(1,h0);   xi(live)   = exp(th);
ze = ones(1,h0);    ze(live)   = exp(et/out.span);
bd = false(1,h0);   bd(live)   = bind;
pij = zeros(n0, h0+1);
pij(alive, [false, live]) = P;
pij(:,1) = 1 - sum(pij(:,2:end), 2);
out.xi = xi;  out.zeta = ze;  out.binding = bd;  out.pij = pij;
out.K = sigma * pij(:,2:end);
out.phi = phi;  out.logdetSigma = logdet;  out.iter = it;
out.relerrEst = 0.14*(1/min(mk) + 1/(n - sum(mk)));
if min([mk, n - sum(mk)]) < 2
    line_warning(mfilename, ['The smallest occupancy is %d, so the expansion is only ' ...
        'qualitative here (estimated relative error %.0f%%); cache_erec is exact.\n'], ...
        min([mk, n - sum(mk)]), 100*out.relerrEst);
end
if strcmp(costmode,'atmost') && any(bind) && max(exp(et(bind))) > 0.8
    line_warning(mfilename, ['A binding cost cap has zeta = %.3f, i.e. it sits in the transition ' ...
        'between the binding and slack regimes where the geometric resummation below the cap ' ...
        'is not sharp; the reported relative error does not cover it.\n'], max(exp(et(bind))));
end
end

% ---------------------------------------------------------------------------

function [th, it, P, D] = local_saddle(G, tgt)
% Size-free saddle: solve sum_i gamma_ij xi_j / D_i = m_j by Newton on
% theta = log xi. The objective sum_i log(1+sum_l gamma_il e^{theta_l}) - m'*theta
% is strictly convex, so the root is unique and damped Newton converges globally.
[~, h] = size(G);
th = local_theta0(G, tgt);
for it = 1:200
    [P, D] = local_occupancy(G, zeros(1,size(G,1)), th, zeros(1,h));
    g = sum(P,1) - tgt;
    if norm(g, inf) <= 1e-12*max(1, norm(tgt, inf))
        break
    end
    H = local_hessian_cost(P, zeros(1,size(G,1)), []);
    d = local_newton_step(H, g);
    th = th + local_damp(d)*d;
end
[P, D] = local_occupancy(G, zeros(1,size(G,1)), th, zeros(1,h));
end

function [th, et, bind, it, P, D] = local_saddle_cost(G, tgtm, sg, tgtk, cumulative)
% Cost-constrained saddle. Minimises the convex dual
%   f(theta,eta) = sum_i log D_i - m'*theta - k'*eta
% over eta <= 0 when the caps are cumulative, so that complementary slackness
% selects the binding lists; over all of R^{2h} when the cost is resolved exactly.
[n, h] = size(G);
th = local_theta0(G, tgtm);
et = zeros(1, h);
bind = true(1, h);
for it = 1:200
    [P, D] = local_occupancy(G, sg, th, et);
    gth = sum(P,1) - tgtm;
    get = sg*P - tgtk;
    if cumulative
        bind = (et < 0) | (get > 0);    % at eta_j = 0 the cap binds when the mean cost exceeds it
    end
    ix = reshape(find(bind), 1, []);
    g  = [gth, get(ix)];
    if norm(g, inf) <= 1e-12*max(1, max(norm(tgtm,inf), norm(tgtk,inf)))
        break
    end
    H = local_hessian_cost(P, sg, ix);
    d = local_newton_step(H, g);
    step = local_damp(d);
    fcur = local_obj(G, sg, th, et, tgtm, tgtk, bind);
    % Backtrack until the dual decreases. The slack is essential, not cosmetic:
    % Newton reaches the floating-point floor of f in a handful of steps, and a
    % strict test then rejects every step and halves to zero without converging.
    ftol = 1e-12*(1 + abs(fcur));
    for ls = 1:40
        thn = th + step*d(1:h);
        etn = et;
        etn(ix) = et(ix) + step*d(h+1:end);
        if cumulative
            etn = min(etn, 0);
        end
        if local_obj(G, sg, thn, etn, tgtm, tgtk, bind) <= fcur + ftol
            break
        end
        step = step/2;
    end
    moved = max([abs(thn - th), abs(etn - et)]);
    th = thn;
    et = etn;
    if moved <= 1e-13                    % the iterate can no longer move: at the floor
        break
    end
end
[P, D] = local_occupancy(G, sg, th, et);
if cumulative
    bind = et < 0;
end
end

function th = local_theta0(G, tgt)
n = size(G,1);
slack = max(1 - sum(tgt)/n, 1e-9);
th = log(max(tgt, 1e-12) ./ max(sum(G,1)*slack, 1e-12));
end

function [P, D] = local_occupancy(G, sg, th, et)
% pi_ij = gamma_ij xi_j zeta_j^{sigma_i} / D_i, D_i = 1 + sum_j gamma_ij xi_j zeta_j^{sigma_i}
A = G .* exp(th + sg(:)*et);
D = 1 + sum(A, 2);
P = A ./ D;
end

function f = local_obj(G, sg, th, et, tgtm, tgtk, bind)
% Summed elementwise, not as a matrix product: at h = 1 with no binding list
% tgtk(bind) is 0x0 and (0x1)*(1x0) is EMPTY, which would make f [] and leave
% every line-search comparison false.
[~, D] = local_occupancy(G, sg, th, et);
f = sum(log(D)) - sum(tgtm.*th) - sum(tgtk(bind).*et(bind));
end

function H = local_hessian_cost(P, sg, ix)
% grad^2 Psi in (theta, eta), restricted to the free eta coordinates ix. Each
% block is sum_i w_i Q^i with Q^i = diag(pi_i) - pi_i pi_i' and w = 1, sigma, sigma^2.
Qw  = @(w) diag(sum(w(:).*P, 1)) - P.'*(w(:).*P);
Htt = Qw(ones(size(P,1),1));
if isempty(ix)
    H = Htt;
    return
end
Hte = Qw(sg(:));
Hee = Qw(sg(:).^2);
H = [Htt, Hte(:,ix); Hte(ix,:), Hee(ix,ix)];
end

function d = local_newton_step(H, g)
d = -(H\g.');
d = d.';
if ~all(isfinite(d))
    line_error(mfilename, ['The saddle-point Newton step is not finite. With item sizes this ' ...
        'is the rank-h degeneracy of the size-tilted Hessian: the sizes must genuinely vary ' ...
        'for the cost coordinate to carry information.']);
end
end

function step = local_damp(d)
step = 1;
while max(abs(step*d)) > 2               % keep the tilts within a factor e^2 per iteration
    step = step/2;
end
end

function ld = local_logdet(H)
[R, p] = chol((H + H.')/2);
if p ~= 0
    line_error(mfilename, ['The saddle-point Hessian is not positive definite; the ray map is ' ...
        'singular here. With item sizes this happens when the sizes do not vary over the items ' ...
        'the cache can hold, in which case the cost cap carries no information.']);
end
ld = 2*sum(log(diag(R)));
end
