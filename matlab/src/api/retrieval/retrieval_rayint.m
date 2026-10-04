%{ @file retrieval_rayint.m
 %  @brief Ray (WKB) asymptotic expansion of the list-based cache normalizing constant
 %
 %  @author LINE Development Team
%}

%{
 % @brief Approximates the cache normalizing constant E(m) by the ray/eikonal expansion
 %
 % @details
 % Evaluates the geometrical-optics (WKB) asymptotic expansion of the normalizing
 % constant computed exactly by cache_erec, and returns it in the SAME
 % normalization, so that the two are interchangeable:
 %
 %   E(m,n) = E(m,n-1) + sum_j gamma_{n,j} m_j E(m-1_j,n-1),   E(0,0)=1
 %
 % Writing E = prod_j m_j! * Et, the relaxation Et ~ H exp(phi/eps) with
 % n = y/eps and m_j = x_j/eps turns the recursion into the eikonal equation
 % e^{phi_y} = 1 + sum_j gamma_j(y) e^{-phi_j}, whose rays carry constants
 % xi_j = e^{-phi_j}. With S(v) = 1 + sum_j gamma_j(v) xi_j the solution is
 %
 %   x_j    = int_0^y gamma_j(v) xi_j / S(v) dv          (the saddle conditions)
 %   phi    = int_0^y log S(v) dv - sum_j x_j log xi_j
 %   H      = (2 pi)^{-h/2} sqrt(S(y)/S(0)) / sqrt(prod_j xi_j * det A)
 %   A_{ik} = d x_i / d xi_k
 %
 % and E ~ prod_j m_j! * eps^{h/2} H exp(phi/eps). Two evaluations are offered,
 % selected by the type of the first argument.
 %
 % DISCRETE (gamma given as an n x h matrix). The ray integrals are the sums
 % they discretize and the expansion collapses to the Laplace form
 %
 %   Et ~ (2 pi)^{-h/2} exp(sum_k log D_k - sum_j m_j log xi_j) / sqrt(det Sigma)
 %
 % with D_k = 1 + sum_j gamma_{k,j} xi_j, sum_k gamma_{k,j} xi_j / D_k = m_j and
 % Sigma = A * diag(xi) the Hessian in log xi. This is the default and is the
 % more accurate of the two; the sqrt(S(y)/S(0)) factor above is exactly the
 % Euler-Maclaurin term relating sum_k to int dv and is already accounted for.
 %
 % CONTINUUM (gamma given as a function handle v -> [numel(v) x h] on v in [0,1]).
 % The ray integrals are evaluated by composite Simpson quadrature on the profile
 % itself, which is the form written in the note; it costs roughly a factor two in
 % accuracy relative to the discrete form but does not need the n rows.
 %
 % ACCURACY. The relative error is O(1/n) at fixed occupancy but is governed by the
 % smallest occupancy rather than by n, tracking
 %   0.14 * ( 1/min_j m_j + 1/(n - sum_j m_j) ),
 % so a per cent needs every m_j and n - sum_j m_j above about 15 and a part in a
 % thousand needs them above about 150. The estimate is returned in OUT.relerrEst.
 % Lists with m_j = 0 contribute nothing and are dropped before the saddle is
 % solved, so they neither degrade the expansion nor make it singular.
 %
 % This covers the no-fetch (q=0) case, i.e. the same quantity as cache_erec. The
 % delayed-hit extension carrying the fetch coordinates is not implemented: its
 % eikonal is known but its amplitude has not been derived, and a partial
 % amplitude would be silently wrong rather than merely approximate.
 %
 % @par Syntax:
 % @code
 % E = retrieval_rayint(gamma, m)                 % gamma is n x h
 % E = retrieval_rayint(gfun, m, n)               % gfun is @(v) -> numel(v) x h
 % [E, logE, out] = retrieval_rayint(...)
 % [...] = retrieval_rayint(..., nquad)           % Simpson nodes, continuum form
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Access factors gamma(k,j), n x h; or a handle @(v) returning numel(v) x h for v in [0,1]
 % <tr><td>m<td>Cache list capacities, 1 x h
 % <tr><td>n<td>Number of items (required, and used only, with a function handle)
 % <tr><td>nquad<td>(Optional) composite Simpson nodes for the continuum form (default 4097, forced odd)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>E<td>Normalizing constant, same normalization as cache_erec (Inf on overflow; use logE)
 % <tr><td>logE<td>Natural logarithm of E, safe for large n
 % <tr><td>out<td>Struct with the ray quantities: xi, phi, logdetSigma, S0, Sy, method, relerrEst, iter
 % </table>
 %
 % @par References:
 % G. Casale, "Accelerating Performance Inference over Closed Systems by
 % Asymptotic Methods", ACM SIGMETRICS, 2017 (asymptotic expansion of the
 % normalizing constant integral); the ray/eikonal form used here follows the
 % geometrical-optics construction for the RR cache.
 %
 % @see cache_erec, retrieval_nc, retrieval_mva, retrieval_fpi
%}
function [E, logE, out] = retrieval_rayint(gamma, m, n, nquad)

if nargin < 2
    line_error(mfilename, 'Two arguments are required: the access factors and the list capacities.');
end
m = m(:).';
if any(m < 0) || any(abs(m - round(m)) > 0)
    line_error(mfilename, 'List capacities must be non-negative integers.');
end

isprofile = isa(gamma, 'function_handle');
if isprofile
    if nargin < 3 || isempty(n)
        line_error(mfilename, 'The number of items n must be given when the access factors are a function handle.');
    end
    if nargin < 4 || isempty(nquad)
        nquad = 4097;
    end
    nquad = max(5, 2*floor(nquad/2) + 1);   % composite Simpson needs an odd node count
    h = size(gamma(0), 2);
else
    if ~ismatrix(gamma) || isempty(gamma)
        line_error(mfilename, 'The access factors must be a non-empty n x h matrix.');
    end
    [n, h] = size(gamma);
    if nargin >= 3 && ~isempty(n) && nargin >= 3
        % n is ignored for the matrix form; its value is size(gamma,1)
    end
end
if numel(m) ~= h
    line_error(mfilename, 'The capacity vector must have one entry per cache list (%d given, %d expected).', numel(m), h);
end

out = struct('xi', [], 'phi', NaN, 'logdetSigma', NaN, 'S0', NaN, 'Sy', NaN, ...
             'method', '', 'relerrEst', NaN, 'iter', 0);

% --- boundaries, matching cache_erec ---
if sum(m) > n
    E = 0; logE = -Inf; out.method = 'boundary'; return
end
if sum(m) == 0
    E = 1; logE = 0; out.method = 'boundary'; return
end
if sum(m) == n
    line_error(mfilename, ['The expansion requires sum(m) < n: at sum(m) = n the saddle point ' ...
        'escapes to infinity. Use cache_erec for a full cache.']);
end

% --- lists with zero capacity contribute nothing: drop them ---
keep = m > 0;
mk = m(keep);
hk = numel(mk);
if isprofile
    gfun = @(v) subcols(gamma(v), keep);
else
    Gk = gamma(:, keep);
end

% --- saddle point, then the expansion ---
if isprofile
    v = linspace(0, 1, nquad).';
    w = ones(nquad,1); w(2:2:end-1) = 4; w(3:2:end-2) = 2; w = w/(3*(nquad-1));
    G = gfun(v);
    if size(G,1) ~= nquad
        line_error(mfilename, 'The access-factor handle must return one row per evaluation point.');
    end
    x = mk/n;
    [xi, it] = local_saddle(G, w, x);
    S  = 1 + G*xi(:);
    phi = w.'*log(S) - x(:).'*log(xi(:));
    Sig = local_hessian(G, w, xi);
    S0 = 1 + gfun(0)*xi(:);
    Sy = 1 + gfun(1)*xi(:);
    logdet = local_logdet(Sig);
    logEt = -(hk/2)*log(n) - (hk/2)*log(2*pi) + n*phi - 0.5*logdet + 0.5*log(Sy/S0);
    out.method = 'rayint';
    out.phi = phi;  out.S0 = S0;  out.Sy = Sy;
else
    w = ones(n,1);
    [xi, it] = local_saddle(Gk, w, mk);
    D   = 1 + Gk*xi(:);
    phi = sum(log(D)) - mk(:).'*log(xi(:));
    Sig = local_hessian(Gk, w, xi);
    logdet = local_logdet(Sig);
    logEt = -(hk/2)*log(2*pi) + phi - 0.5*logdet;
    out.method = 'saddle';
    out.phi = phi;  out.S0 = D(1);  out.Sy = D(end);
end

logE = logEt + sum(gammaln(m + 1));     % back to the cache_erec normalization
E = exp(logE);

xifull = zeros(1,h); xifull(keep) = xi;
out.xi = xifull;
out.logdetSigma = logdet;
out.iter = it;
out.relerrEst = 0.14*(1/min(mk) + 1/(n - sum(m)));
if min([mk, n - sum(m)]) < 2
    line_warning(mfilename, ['The smallest occupancy is %d, so the expansion is only ' ...
        'qualitative here (estimated relative error %.0f%%); cache_erec is exact.\n'], ...
        min([mk, n - sum(m)]), 100*out.relerrEst);
end
end

% ---------------------------------------------------------------------------

function G = subcols(G, keep)
G = G(:, keep);
end

function [xi, it] = local_saddle(G, w, tgt)
% Solve sum_k w_k gamma_{k,j} xi_j / (1 + sum_l gamma_{k,l} xi_l) = tgt_j by Newton
% on theta = log xi. The objective sum_k w_k log(1+sum_l gamma_{k,l} e^{theta_l})
% - tgt'*theta is strictly convex, so the root is unique and Newton is globally
% convergent under damping.
tgt = tgt(:);
gb  = (w.'*G).';
slack = max(1 - sum(tgt)/sum(w), 1e-9);
th  = log(max(tgt, 1e-12) ./ max(gb*slack, 1e-12));
it  = 0;
for it = 1:200
    xi = exp(th);
    a  = G .* xi(:).';
    S  = 1 + sum(a, 2);
    g  = (sum(w.*a./S, 1)).' - tgt;
    if norm(g, inf) <= 1e-12*max(1, norm(tgt, inf))
        break
    end
    H = local_hessian(G, w, xi);
    d = -(H\g);
    if ~all(isfinite(d))
        line_error(mfilename, 'The saddle-point Newton step is not finite; check that the access factors are positive.');
    end
    step = 1;
    while max(abs(step*d)) > 2      % keep xi within a factor e^2 per iteration
        step = step/2;
    end
    th = th + step*d;
end
xi = exp(th);
end

function H = local_hessian(G, w, xi)
% Hessian in theta = log xi; equals A*diag(xi) with A_{ik} = d x_i / d xi_k
a = G .* xi(:).';
S = 1 + sum(a, 2);
H = diag(sum(w.*a./S, 1)) - (a./S).'*(w.*(a./S));
end

function ld = local_logdet(H)
[R, p] = chol((H + H.')/2);
if p ~= 0
    line_error(mfilename, 'The saddle-point Hessian is not positive definite; the ray map is singular here.');
end
ld = 2*sum(log(diag(R)));
end
