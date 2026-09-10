function [L, W, Ca, Cd, lambda, rho, X, iter] = me_cqn(M, R, N, mu, Cs, P, c, refstat, insens, options)
%ME_CQN Maximum Entropy algorithm for Closed Queueing Networks
%
% Implements the two-stage ME algorithm from Kouvatsos (1994) "Entropy
% Maximisation and Queueing Network Models", Section 3.3, for closed
% multiclass networks of G/G/1 and G/G/inf queues.
%
% Stage 1 solves a pseudo-open network (no external arrivals) subject to
% job flow conservation and the population constraints sum_i L(i,r)=N(r),
% using the GE-type fixed point of the open algorithm (Section 3.2).
% Stage 2 builds the ME product-form solution (3.8) from the Lagrangian
% coefficients of Stage 1, computes the normalising constant Z(N) by a
% multiclass convolution, and iterates the flow (work rate) equations
% until the class throughputs implied by the closed ME solution agree
% with those used to parametrise the building blocks.
%
% INPUTS:
%   M       - Number of queues (stations)
%   R       - Number of job classes
%   N       - Class populations [1 x R]
%   mu      - Service rates [M x R matrix]
%   Cs      - Service scv [M x R matrix]
%   P       - Routing probability matrix [M x M x R], P(j,i,r) = p_ji,r
%   c       - (optional) Servers per queue [M x 1]; Inf marks an
%             infinite-server (IS) queue; finite values must be 1
%             (default: ones(M,1))
%   refstat - (optional) Reference station per class [1 x R] for visit
%             normalisation (default: first station visited by the class)
%   options - (optional) struct with fields:
%             .tol     - convergence tolerance (default: 1e-6)
%             .maxiter - maximum iterations (default: 1000)
%             .verbose - print iteration info (default: false)
%
% OUTPUTS:
%   L      - Mean queue lengths [M x R matrix] (sum_i L(i,r) = N(r))
%   W      - Mean response times [M x R matrix], W = L ./ lambda
%   Ca     - Arrival scv at each queue [M x R matrix] (pseudo-open)
%   Cd     - Departure scv at each queue [M x R matrix] (pseudo-open)
%   lambda - Class throughputs at each queue [M x R], lambda(i,r)=X(r)*v(i,r)
%   rho    - Utilizations [M x R matrix] from the closed ME solution
%            (mean jobs in service for IS queues)
%   X      - Class throughputs at the reference station [1 x R]
%   iter   - Total number of fixed-point iterations
%
% Reference:
%   D.D. Kouvatsos, "Entropy Maximisation and Queueing Network Models",
%   Annals of Operations Research, 48:63-126, 1994. Section 3.3.

if nargin < 7 || isempty(c)
    c = ones(M, 1);
end
if nargin < 8 || isempty(refstat)
    refstat = zeros(1, R);
end
if nargin < 9 || isempty(insens)
    insens = false(M, 1);
end
if nargin < 10
    options = struct();
end
insens = logical(insens(:));
if ~isfield(options, 'tol')
    options.tol = 1e-6;
end
if ~isfield(options, 'maxiter')
    options.maxiter = 1000;
end
if ~isfield(options, 'verbose')
    options.verbose = false;
end
c = c(:);
N = N(:)';

% Feedback correction (as in the open algorithm)
P_eff = P;
mu_eff = mu;
Cs_eff = Cs;
selfp = zeros(M, R);
for i = 1:M
    for r = 1:R
        pii = P(i, i, r);
        if pii > 0
            selfp(i, r) = pii;
            mu_eff(i, r) = mu(i, r) * (1 - pii);
            Cs_eff(i, r) = pii + (1 - pii) * Cs(i, r);
            P_eff(i, :, r) = P(i, :, r) / (1 - pii);
            P_eff(i, i, r) = 0;
        end
    end
end

% Visit ratios from the original routing (visit-inclusive), normalised at
% the reference station of each class
V = zeros(M, R);
for r = 1:R
    A = eye(M) - P(:, :, r)';
    if refstat(r) <= 0
        ref = find(mu(:, r) > 0, 1);
    else
        ref = refstat(r);
    end
    A(ref, :) = 0;
    A(ref, ref) = 1;
    b = zeros(M, 1);
    b(ref) = 1;
    V(:, r) = A \ b;
    V(abs(V(:, r)) < 1e-14, r) = 0;
    refstat(r) = ref;
end

% Stage 1: pseudo-open network, find X such that sum_i L(i,r) = N(r)
% Initialise throughputs at half the single-station capacity bound
X = zeros(1, R);
for r = 1:R
    capr = Inf;
    for i = 1:M
        if ~isinf(c(i)) && V(i, r) > 0 && mu(i, r) > 0
            capr = min(capr, mu(i, r) / V(i, r));
        end
    end
    if isinf(capr) % IS-only class
        capr = 1;
    end
    X(r) = 0.5 * capr / R;
end

Ca = ones(M, R);
iter = 0;
Lpo = zeros(M, R);
rho_po = zeros(M, R);
Cd = ones(M, R);
% Stage 1 only seeds the Stage 2 flow iteration, so a moderate iteration
% budget suffices; updates are damped and step-clamped to avoid limit
% cycles when the population target sits beyond the saturation cap
maxit1 = min(options.maxiter, 100);
for it1 = 1:maxit1
    iter = iter + 1;
    X = me_cqn_capacity_cap(X, V, mu, c, M, R);
    lambda = V .* repmat(X, M, 1);
    [Lpo, Ca, Cd, rho_po] = me_cqn_pseudoopen(M, R, lambda, mu, mu_eff, Cs_eff, P_eff, selfp, c, insens, Ca, options);
    Ltot = sum(Lpo, 1);
    err1 = 0;
    for r = 1:R
        if N(r) > 0 && Ltot(r) > 0
            err1 = max(err1, abs(Ltot(r) - N(r)) / N(r));
        end
    end
    if options.verbose
        fprintf('Stage 1 iteration %d: max population error = %e\n', int32(it1), err1);
    end
    if err1 < options.tol
        break;
    end
    Xold = X;
    for r = 1:R
        if Ltot(r) > 0
            fac = min(max((N(r) / Ltot(r))^0.5, 0.25), 4); % clamped step
            X(r) = 0.5 * X(r) + 0.5 * X(r) * fac;          % damped update
        end
    end
    % Stall guard: the stability cap can bind before the population
    % target is met (bottleneck saturation); stop when X no longer moves
    if max(abs(me_cqn_capacity_cap(X, V, mu, c, M, R) - Xold) ./ max(Xold, 1e-12)) < options.tol
        break;
    end
end

% Stage 2: closed ME solution by convolution, iterated on the flow
% (work rate) equations
sz = N + 1;
[Dec, PIdx] = me_cqn_lattice(N, R);
L = Lpo;
rho = rho_po;
for it2 = 1:options.maxiter
    iter = iter + 1;
    % Lagrangian coefficient functions f_i over the population lattice
    F = me_cqn_coefficients(M, R, PIdx, Dec, sz, Lpo, rho_po, lambda, mu_eff, Cs_eff, Ca, c, selfp);
    % Convolution and per-station marginals
    [L, U] = me_cqn_convolve(M, R, N, PIdx, Dec, sz, F, c);
    % Utilization split by pseudo-open per-class load; implied throughputs
    % from the work rate theorem, averaged with visit weights
    rho = zeros(M, R);
    Xhat = zeros(1, R);
    for r = 1:R
        num = 0;
        den = 0;
        for i = 1:M
            if lambda(i, r) > 0
                if isinf(c(i))
                    rho(i, r) = L(i, r);
                    % IS work rate: lambda_eff = L*mu_eff, revisits add 1/(1-p)
                    num = num + L(i, r) * mu_eff(i, r) / (1 - selfp(i, r));
                else
                    rho_i = sum(rho_po(i, :));
                    if rho_i > 0
                        rho(i, r) = U(i) * rho_po(i, r) / rho_i;
                    end
                    num = num + rho(i, r) * mu(i, r);
                end
                den = den + V(i, r);
            end
        end
        if den > 0
            Xhat(r) = num / den;
        end
    end
    err2 = 0;
    for r = 1:R
        if X(r) > 0
            err2 = max(err2, abs(Xhat(r) - X(r)) / X(r));
        end
    end
    if options.verbose
        fprintf('Stage 2 iteration %d: max flow error = %e\n', int32(it2), err2);
    end
    if err2 < options.tol
        break;
    end
    % Damped throughput update with stability cap, then refresh the
    % pseudo-open decomposition at the new flows
    Xold = X;
    X = 0.5 * X + 0.5 * Xhat;
    X = me_cqn_capacity_cap(X, V, mu, c, M, R);
    % Stall guard: when the stability cap binds, X stops moving even
    % though the residual flow error stays above tolerance
    if max(abs(X - Xold) ./ max(Xold, 1e-12)) < options.tol
        break;
    end
    lambda = V .* repmat(X, M, 1);
    [Lpo, Ca, Cd, rho_po] = me_cqn_pseudoopen(M, R, lambda, mu, mu_eff, Cs_eff, P_eff, selfp, c, insens, Ca, options);
end

if it2 == options.maxiter && err2 >= options.tol
    warning('me_cqn:noconverge', 'Did not converge within %d iterations (flow error=%e)', int32(options.maxiter), err2);
end

% Response times by Little's law on the visit-inclusive throughputs
lambda = V .* repmat(X, M, 1);
W = zeros(M, R);
for i = 1:M
    for r = 1:R
        if lambda(i, r) > 0
            W(i, r) = L(i, r) / lambda(i, r);
        end
    end
end

end

%% ------------------------------------------------------------------------
function X = me_cqn_capacity_cap(X, V, mu, c, M, R)
% Scales the class throughputs uniformly so that every single-server
% queue in the pseudo-open network remains stable
maxrho = 0;
for i = 1:M
    if ~isinf(c(i))
        rho_i = 0;
        for r = 1:R
            if V(i, r) > 0 && mu(i, r) > 0
                rho_i = rho_i + X(r) * V(i, r) / mu(i, r);
            end
        end
        maxrho = max(maxrho, rho_i);
    end
end
if maxrho >= 0.999
    X = X * (0.999 / maxrho);
end
end

function [L, Ca, Cd, rho] = me_cqn_pseudoopen(M, R, lambda, mu, mu_eff, Cs_eff, P_eff, selfp, c, insens, Ca, options)
% GE-type fixed point of the open algorithm (Section 3.2) on the
% pseudo-open network: no external arrivals, flows given by lambda.
% The flow scvs are computed on the class-composed (aggregate) streams
% and disaggregated per class by thinning, following the class
% composition and disaggregation principle of the closed ME algorithm
% (Kouvatsos 1994, Section 3.3 discussion); this keeps the closed
% solution consistent with the single-class one when classes are
% statistically identical.
lambda_eff = lambda .* (1 - selfp);
rho = zeros(M, R);
for i = 1:M
    for r = 1:R
        if mu(i, r) > 0
            if isinf(c(i))
                rho(i, r) = lambda_eff(i, r) / mu_eff(i, r);
            else
                rho(i, r) = lambda(i, r) / mu(i, r);
            end
        end
    end
end
% Class composition per station: aggregate flow, service process moments
% and flow-weighted aggregate routing
lam_a = sum(lambda_eff, 2)';
mu_a = zeros(1, M);
Cs_a = ones(1, M);
for i = 1:M
    if lam_a(i) > 0
        ES = 0;
        ES2 = 0;
        for u = 1:R
            if lambda_eff(i, u) > 0 && mu_eff(i, u) > 0
                wu = lambda_eff(i, u) / lam_a(i);
                ES = ES + wu / mu_eff(i, u);
                ES2 = ES2 + wu * (Cs_eff(i, u) + 1) / mu_eff(i, u)^2;
            end
        end
        if ES > 0
            mu_a(i) = 1 / ES;
            Cs_a(i) = ES2 / ES^2 - 1;
        end
    end
end
Pa = zeros(M, M);
for j = 1:M
    if lam_a(j) > 0
        for i = 1:M
            num = 0;
            for r = 1:R
                if lambda_eff(j, r) > 0
                    num = num + lambda_eff(j, r) * P_eff(j, i, r);
                end
            end
            Pa(j, i) = num / lam_a(j);
        end
    end
end
% Fixed point on the aggregate arrival scvs
Ca_a = ones(1, M);
for i = 1:M
    if lam_a(i) > 0 && any(lambda_eff(i, :) > 0)
        wr = find(lambda_eff(i, :) > 0, 1);
        Ca_a(i) = 1 + (Ca(i, wr) - 1) * lam_a(i) / max(lambda_eff(i, wr), realmin); % warm start
    end
end
Cd_a = ones(1, M);
L_a = zeros(1, M);
for it = 1:options.maxiter
    Ca_old = Ca_a;
    for i = 1:M
        if lam_a(i) <= 0
            continue;
        end
        rho_i = sum(rho(i, :));
        if isinf(c(i))
            % GE/GE/inf: L = lambda/mu, departures inherit the arrival scv
            L_a(i) = lam_a(i) / mu_a(i);
            Cd_a(i) = Ca_a(i);
        elseif rho_i < 1
            if insens(i)
                % Insensitive disciplines (PS, LCFS-PR): product-form mql
                L_a(i) = rho_i / (1 - rho_i);
            else
                % Single-class GE/GE/1 mql, eq. (3.6)
                L_a(i) = rho_i * (Ca_a(i) + 1) / 2 + rho_i^2 * (Ca_a(i) + Cs_a(i)) / (2 * (1 - rho_i));
            end
            Cd_a(i) = 2 * L_a(i) * (1 - rho_i) + Ca_a(i) * (1 - 2 * rho_i);
        end
    end
    % GE-type merging, eq. (3.7) with lambda_o = 0, on aggregate flows
    for i = 1:M
        if lam_a(i) > 0
            sum_inv = 0;
            for j = 1:M
                if Pa(j, i) > 0 && lam_a(j) > 0
                    Cdji = 1 + Pa(j, i) * (Cd_a(j) - 1);
                    sum_inv = sum_inv + (lam_a(j) * Pa(j, i) / lam_a(i)) / (Cdji + 1);
                end
            end
            if sum_inv > 0
                Ca_a(i) = -1 + 1 / sum_inv;
            end
        end
    end
    delta = max(abs(Ca_a - Ca_old));
    if delta < options.tol
        break;
    end
end
% Disaggregation: per-class arrival scvs by thinning of the composed
% stream, then per-class mean queue lengths (Section 3.1.1)
L = zeros(M, R);
Cd = ones(M, R);
for i = 1:M
    rho_i = sum(rho(i, :));
    for r = 1:R
        if lambda_eff(i, r) > 0
            pr = lambda_eff(i, r) / lam_a(i);
            Ca(i, r) = 1 + pr * (Ca_a(i) - 1);
            Cd(i, r) = 1 + pr * (Cd_a(i) - 1);
        end
    end
    if isinf(c(i))
        for r = 1:R
            if lambda_eff(i, r) > 0 && mu_eff(i, r) > 0
                L(i, r) = lambda_eff(i, r) / mu_eff(i, r);
            end
        end
    elseif rho_i < 1
        if insens(i)
            % Insensitive disciplines (PS, LCFS-PR): product-form mql
            for r = 1:R
                if lambda_eff(i, r) > 0 && mu_eff(i, r) > 0
                    L(i, r) = rho(i, r) / (1 - rho_i);
                end
            end
        else
            resid = 0;
            for u = 1:R
                if lambda_eff(i, u) > 0 && mu_eff(i, u) > 0
                    resid = resid + lambda_eff(i, u) * (Cs_eff(i, u) + Ca(i, u)) / mu_eff(i, u)^2;
                end
            end
            for r = 1:R
                if lambda_eff(i, r) > 0 && mu_eff(i, r) > 0
                    L(i, r) = rho(i, r) * (Ca(i, r) + 1) / 2 + lambda_eff(i, r) * resid / (2 * (1 - rho_i));
                end
            end
        end
    end
end
end

function [Dec, PIdx] = me_cqn_lattice(N, R)
% Enumerates the population lattice {0..N(1)} x ... x {0..N(R)} in mixed
% radix order; Dec(p,:) is the population vector of linear index p
sz = N + 1;
PIdx = prod(sz);
Dec = zeros(PIdx, R);
for p = 1:PIdx
    q = p - 1;
    for r = 1:R
        Dec(p, r) = mod(q, sz(r));
        q = floor(q / sz(r));
    end
end
end

function F = me_cqn_coefficients(M, R, PIdx, Dec, sz, Lpo, rho_po, lambda, mu_eff, Cs_eff, Ca, c, selfp)
% Auxiliary functions f_i(n) of the ME solution (3.8): the right-hand
% sides of (3.2) and (3.4) with the (1-rho) factor removed, evaluated
% from the Stage 1 Lagrangian coefficients. Each f_i is rescaled by its
% maximum for numerical stability (per-station constants cancel in the
% marginal probabilities).
F = zeros(PIdx, M);
lambda_eff = lambda .* (1 - selfp);
lgamma = gammaln(1:(sum(sz - 1) + 2)); % lgamma(k) = log((k-1)!)
for i = 1:M
    if isinf(c(i))
        % GE/GE/inf: f(n) = prod_r prod_{k=1}^{n_r} g_r(k), with g_r(j)
        % from the exact ME solution of the GE/GE/inf queue
        logg = cell(1, R);
        for r = 1:R
            logg{r} = zeros(1, sz(r) - 1);
            for j = 1:(sz(r) - 1)
                if lambda_eff(i, r) > 0 && mu_eff(i, r) > 0
                    gj = (lambda_eff(i, r) * (1 + Cs_eff(i, r)) + (j - 1) * mu_eff(i, r) * (Ca(i, r) - 1)) ...
                        / (j * mu_eff(i, r) * (Ca(i, r) + Cs_eff(i, r)));
                    logg{r}(j) = log(max(gj, 0));
                else
                    logg{r}(j) = -Inf;
                end
            end
        end
        for p = 1:PIdx
            n = Dec(p, :);
            val = 0;
            for r = 1:R
                for j = 1:n(r)
                    val = val + logg{r}(j);
                end
            end
            F(p, i) = exp(val);
        end
        F(1, i) = 1;
    else
        % GE/GE/1 (non-priority): f(n) = ((|n|-1)!/prod_r n_r!)
        %   * sum_r n_r*(g_r*x_r)*x_r^(n_r-1)*prod_{s~=r} x_s^{n_s}
        % with x_r = (L_r-rho_r)/L and g_r*x_r = rho_r*rho/((1-rho)*L)
        rho_i = sum(rho_po(i, :));
        Li = sum(Lpo(i, :));
        x = zeros(1, R);
        gx = zeros(1, R);
        if Li > 0 && rho_i < 1
            for r = 1:R
                if lambda(i, r) > 0
                    x(r) = max(Lpo(i, r) - rho_po(i, r), 0) / Li;
                    gx(r) = rho_po(i, r) * rho_i / ((1 - rho_i) * Li);
                end
            end
        end
        for p = 1:PIdx
            n = Dec(p, :);
            ntot = sum(n);
            if ntot == 0
                F(p, i) = 1;
                continue;
            end
            if any(n > 0 & lambda(i, :) <= 0)
                F(p, i) = 0; % class not visiting this station
                continue;
            end
            logmult = lgamma(ntot) - sum(lgamma(n + 1)); % (|n|-1)!/prod n_r!
            tot = 0;
            for r = 1:R
                if n(r) > 0 && gx(r) > 0
                    % n_r * gx_r * x_r^(n_r-1) * prod_{s~=r} x_s^{n_s}
                    lterm = log(n(r)) + log(gx(r));
                    ok = true;
                    for s = 1:R
                        es = n(s);
                        if s == r
                            es = es - 1;
                        end
                        if es > 0
                            if x(s) > 0
                                lterm = lterm + es * log(x(s));
                            else
                                ok = false;
                                break;
                            end
                        end
                    end
                    if ok
                        tot = tot + exp(logmult + lterm);
                    end
                end
            end
            F(p, i) = tot;
        end
    end
    fmax = max(F(:, i));
    if fmax > 0
        F(:, i) = F(:, i) / fmax;
    end
    F(1, i) = max(F(1, i), realmin); % f_i(0) stays positive after scaling
end
end

function [L, U] = me_cqn_convolve(M, R, N, PIdx, Dec, sz, F, c)
% Computes the normalising constant by convolving the f_i over the
% population lattice, and the per-station marginals by prefix/suffix
% convolutions; returns closed mean queue lengths and busy probabilities
rad = ones(1, R);
for r = 2:R
    rad(r) = rad(r - 1) * sz(r - 1);
end
% Prefix and suffix convolutions
Gpre = cell(1, M + 1);
Gsuf = cell(1, M + 2);
G0 = zeros(PIdx, 1);
G0(1) = 1;
Gpre{1} = G0;
for k = 1:M
    Gpre{k + 1} = me_cqn_convpair(Gpre{k}, F(:, k), PIdx, Dec, rad);
end
Gsuf{M + 1} = G0;
for k = M:-1:1
    Gsuf{k} = me_cqn_convpair(Gsuf{k + 1}, F(:, k), PIdx, Dec, rad);
end
Z = Gpre{M + 1}(PIdx);
L = zeros(M, R);
U = zeros(M, 1);
for i = 1:M
    Grest = me_cqn_convpair(Gpre{i}, Gsuf{i + 1}, PIdx, Dec, rad);
    % Marginal P_i(n) = f_i(n) * Grest(N - n) / Z
    for p = 1:PIdx
        if F(p, i) > 0
            n = Dec(p, :);
            q = 1 + (N - n) * rad';
            pin = F(p, i) * Grest(q) / Z;
            if p > 1
                U(i) = U(i) + pin;
            end
            for r = 1:R
                if n(r) > 0
                    L(i, r) = L(i, r) + n(r) * pin;
                end
            end
        end
    end
end
end

function G2 = me_cqn_convpair(G, f, PIdx, Dec, rad)
% Convolution of a partial normalising constant with one station term
% over the population lattice
G2 = zeros(PIdx, 1);
for p = 1:PIdx % station population n
    if f(p) == 0
        continue;
    end
    n = Dec(p, :);
    for q = 1:PIdx % remainder population m
        if G(q) == 0
            continue;
        end
        m = Dec(q, :);
        t = n + m;
        idx = 1 + t * rad';
        if idx <= PIdx && all(t <= Dec(PIdx, :))
            G2(idx) = G2(idx) + f(p) * G(q);
        end
    end
end
end
