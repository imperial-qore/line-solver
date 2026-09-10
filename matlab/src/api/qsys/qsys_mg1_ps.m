function result = qsys_mg1_ps(lambda, svc, svcparam, varargin)
% RESULT = QSYS_MG1_PS(LAMBDA, ALPHA, T, ...)
% RESULT = QSYS_MG1_PS(LAMBDA, BLST, M1, ...)
%
% Sojourn time distribution of the M/G/1 processor-sharing queue.
%
% Jobs arrive in a Poisson stream of rate LAMBDA at a single egalitarian
% processor-sharing server whose service requirement has Laplace-Stieltjes
% transform bhat(tau) and mean M1. Writing V(x) for the sojourn time of a
% tagged job of service requirement x and RHO = LAMBDA*M1 < 1, Ott (1984)
% and Yashkov (1983) express the conditional transform as
%
%   E[exp(-s V(x))] = (1-RHO) / D(s,x),
%
% where D(s,x) is the inverse Laplace transform, evaluated at x, of
%
%   f(tau;s) = [ (1-RHO)*tau^2 - (1-RHO)*LAMBDA*(1-bhat(tau))*tau
%                + s*RHO*tau - s*LAMBDA*(1-bhat(tau)) ]
%              / [ tau^2 * (tau - s - LAMBDA*(1-bhat(tau))) ].
%
% The transform is exact but implicit, since f must be inverted in tau. For
% phase-type service f(tau;s) is a proper rational function of tau, the
% double pole at tau = 0 cancels, and D(s,x) is obtained in closed form as a
% finite sum of residues (or, for repeated poles, from the matrix exponential
% of the companion realization). This makes the M/PH/1-PS queue, hence every
% service law LINE can fit with a phase-type distribution, exactly solvable.
% For a service transform supplied as a function handle, f is inverted in tau
% numerically on a Bromwich contour placed to the right of the dominant
% singularity tau*(s), the unique root of tau = s + LAMBDA*(1-bhat(tau)) in
% the right half plane, which the fixed-point iteration of that equation
% reaches at geometric rate RHO.
%
% The conditional sojourn time is atomic on the lattice t = (k+1)*x, for
% k = 0,1,2,...: processor sharing gives every job in the system the same amount
% of work, so if the k jobs present on arrival all outlive the tagged job and no
% arrival intervenes, the sojourn is exactly (k+1)*x. For exponential service
% the masses are A_k = (1-RHO)*RHO^k*exp(-k*MU*x)*exp(-LAMBDA*(k+1)*x), the
% k = 0 term being the probability of finding the system empty and sharing it
% with nobody, which is the only one that stays exact for general service. The
% k = 0 atom is removed before inverting in s; the remaining atoms make cdfCond
% jump and leave no density, so pdfCond is NaN on the lattice.
%
% Input:
%   lambda   - Poisson arrival rate, positive scalar
%   svc      - phase-type initial probability vector alpha (1,n), or a
%              function handle bhat(tau) returning the service LST, which
%              must accept complex arguments
%   svcparam - phase-type subgenerator T (n,n) when svc is a vector, or the
%              mean service time m1 when svc is a function handle
%
% Optional name-value pairs:
%   'x'       - service requirements to condition on (default [])
%   's'       - transform arguments at which to tabulate the LST (default [])
%   't'       - times at which to evaluate the sojourn time distribution
%               (default [])
%   'nterms'  - function evaluations per numerical Laplace inversion, odd
%               (default 41)
%   'pdf'     - service density handle, needed to remove the conditioning
%               when svc is a function handle (default [], filled in
%               automatically on the phase-type path)
%
% Output (struct):
%   rho          - utilization LAMBDA*M1
%   m1, m2       - first two moments of the service requirement (m2 is NaN
%                  when only a transform handle is supplied)
%   lstCond      - handle (s,x) -> E[exp(-s V(x))]
%   lstExcess    - handle (s,x) -> E[exp(-s (V(x)-x))], bounded at large s
%   lstUncond    - handle s -> E[exp(-s V)], by quadrature over the density
%   dominantRoot - handle s -> tau*(s)
%   x, s         - the requested grids
%   lstCondVal   - (numel(x),numel(s)) values of lstCond
%   lstUncondVal - (1,numel(s)) values of lstUncond
%   atomCond     - (1,numel(x)) atom (1-RHO)*exp(-LAMBDA*x) at t = x
%   atomUncond   - (1-RHO)*bhat(LAMBDA), mass of the unshared jobs
%   meanCond     - (1,numel(x)) exact conditional mean x/(1-RHO)
%   m2Cond       - (1,numel(x)) conditional second moment, from the
%                  transform derivatives
%   varCond      - (1,numel(x)) conditional variance
%   meanUncond   - exact unconditional mean M1/(1-RHO)
%   m2Uncond     - unconditional second moment
%   varUncond    - unconditional variance
%   t            - the requested time grid
%   pdfCond      - (numel(x),numel(t)) density of V(x), NaN on the lattice
%   cdfCond      - (numel(x),numel(t)) P(V(x) <= t), atom included
%   pdfUncond    - (1,numel(t)) density of V
%   cdfUncond    - (1,numel(t)) P(V <= t)
%
% References:
%   T. J. Ott, "The sojourn-time distribution in the M/G/1 queue with
%   processor sharing", J. Appl. Prob. 21(2), 1984, pp. 360-378.
%   S. F. Yashkov, "A derivation of response time distribution for an M/G/1
%   processor-sharing queue", Probl. Contr. Inform. Theory 12, 1983,
%   pp. 133-148.
%   Q. Zhen, C. Knessl, "Asymptotic expansions for the sojourn time
%   distribution in the M/G/1-PS queue", Math. Meth. Oper. Res. 74, 2011,
%   equations (2.2)-(2.5).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

opts = struct('x', [], 's', [], 't', [], 'nterms', 41, 'pdf', []);
if mod(numel(varargin), 2) ~= 0
    line_error(mfilename, 'optional arguments must be name-value pairs');
end
for i = 1:2:numel(varargin)
    name = varargin{i};
    if ~ischar(name) || ~isfield(opts, name)
        line_error(mfilename, 'unknown option, expected one of x, s, t, nterms, pdf');
    end
    opts.(name) = varargin{i+1};
end
if ~isscalar(lambda) || ~isfinite(lambda) || lambda <= 0
    line_error(mfilename, 'lambda must be a finite positive scalar');
end
nterms = opts.nterms;
if mod(nterms, 2) == 0 || nterms < 11
    line_error(mfilename, 'nterms must be an odd integer of at least 11');
end

isPH = ~isa(svc, 'function_handle');
if isPH
    alpha = svc(:)';
    T = svcparam;
    n = numel(alpha);
    if ~isequal(size(T), [n n])
        line_error(mfilename, 'T must be %d x %d to match alpha', n, n);
    end
    if abs(sum(alpha) - 1) > 1e-8 || any(alpha < -1e-12)
        line_error(mfilename, 'alpha must be a probability vector');
    end
    exitrate = -T * ones(n, 1);
    if any(exitrate < -1e-10) || any(diag(T) >= 0)
        line_error(mfilename, 'T must be a proper phase-type subgenerator');
    end
    m1 = alpha * ((-T) \ ones(n, 1));
    m2 = 2 * alpha * ((-T) \ ((-T) \ ones(n, 1)));
    bhat = @(tau) alpha * ((tau * eye(n) - T) \ exitrate);
    bpdf = @(y) alpha * expm(T * y) * exitrate;
    % Faddeev-LeVerrier gives det(tau*I-T) and the adjugate in one sweep, so
    % bhat(tau) = nb(tau)/db(tau) as polynomials of degree n-1 and n
    db = [1, zeros(1, n)];
    nb = zeros(1, n);
    Mk = eye(n);
    for k = 1:n
        nb(k) = alpha * Mk * exitrate;
        TM = T * Mk;
        db(k+1) = -trace(TM) / k;
        Mk = TM + db(k+1) * eye(n);
    end
else
    bhat = svc;
    m1 = svcparam;
    if ~isscalar(m1) || ~isfinite(m1) || m1 <= 0
        line_error(mfilename, 'the mean service time must be a finite positive scalar');
    end
    m2 = NaN;
    bpdf = opts.pdf;
    db = [];
    nb = [];
end

rho = lambda * m1;
if rho >= 1
    line_error(mfilename, sprintf('system is unstable: utilization %.6f >= 1', rho));
end

rootfun = @(s) local_root(s, lambda, bhat, rho);
if isPH
    denom = @(s, x) local_denom_ph(s, x, lambda, rho, db, nb);
    ymax = max(40 * m1, 40 / min(-real(eig(T))));
else
    shift = 0.25 * (1 - rho) / m1;
    denom = @(s, x) local_denom_gen(s, x, lambda, rho, bhat, rootfun, shift, nterms);
    ymax = 60 * m1;
end
lstCond = @(s, x) local_lst_cond(s, x, rho, denom, 0);
% transform of the excess V(x)-x, which is what stays bounded as s grows
lstExcess = @(s, x) local_lst_cond(s, x, rho, denom, 1);
hasPdf = ~isempty(bpdf);
% fixed panelled Gauss-Legendre rule for removing the conditioning: the nodes do
% not move with s, so the service density is sampled only once
[yq, wq] = local_quad_nodes(ymax);
if hasPdf
    bq = arrayfun(bpdf, yq);
else
    bq = [];
end
lstUncond = @(s) local_lst_uncond(s, rho, denom, yq, wq, bq);

result = struct();
result.rho = rho;
result.m1 = m1;
result.m2 = m2;
result.lstCond = lstCond;
result.lstExcess = lstExcess;
result.lstUncond = lstUncond;
result.dominantRoot = rootfun;
result.meanUncond = m1 / (1 - rho);
result.atomUncond = (1 - rho) * bhat(lambda);

x = opts.x(:)';
s = opts.s(:)';
t = opts.t(:)';
result.x = x;
result.s = s;
result.t = t;

result.lstCondVal = zeros(numel(x), numel(s));
for i = 1:numel(x)
    for j = 1:numel(s)
        result.lstCondVal(i, j) = real(lstCond(s(j), x(i)));
    end
end
result.lstUncondVal = NaN(1, numel(s));
if hasPdf
    for j = 1:numel(s)
        result.lstUncondVal(j) = real(lstUncond(s(j)));
    end
end

result.atomCond = (1 - rho) * exp(-lambda * x);
result.meanCond = x / (1 - rho);
result.m2Cond = zeros(1, numel(x));
result.varCond = zeros(1, numel(x));
for i = 1:numel(x)
    result.m2Cond(i) = local_second_moment(@(u) lstCond(u, x(i)), ...
        result.meanCond(i), result.meanUncond);
    result.varCond(i) = result.m2Cond(i) - result.meanCond(i)^2;
end
result.m2Uncond = NaN;
result.varUncond = NaN;
if hasPdf
    % integrate the conditional second moment, which is far better conditioned
    % than differentiating the quadrature that removes the conditioning
    m2q = zeros(1, numel(yq));
    for k = 1:numel(yq)
        m2q(k) = local_second_moment(@(u) lstCond(u, yq(k)), ...
            yq(k) / (1 - rho), result.meanUncond);
    end
    result.m2Uncond = sum(wq .* bq .* m2q);
    result.varUncond = result.m2Uncond - result.meanUncond^2;
end

if ~isempty(t) && ~isPH
    line_error(mfilename, ['the sojourn time distribution needs phase-type service, since ' ...
        'inverting a numerically inverted transform is unstable in double precision; ' ...
        'with a transform handle only the LST and its moments are available']);
end
result.pdfCond = zeros(numel(x), numel(t));
result.cdfCond = zeros(numel(x), numel(t));
for i = 1:numel(x)
    atom = result.atomCond(i);
    xi = x(i);
    % V(x) >= x with an atom at x, so invert the excess V(x)-x net of its atom
    gpdf = @(u) lstExcess(u, xi) - atom;
    gcdf = @(u) gpdf(u) / u;
    for j = 1:numel(t)
        if t(j) < xi
            continue
        elseif t(j) == xi
            result.cdfCond(i, j) = atom;
            continue
        end
        result.cdfCond(i, j) = real(local_ilt(gcdf, t(j) - xi, nterms)) + atom;
        % V(x) is atomic on the lattice (k+1)*x, where no density exists
        ratio = t(j) / xi;
        if abs(ratio - round(ratio)) < 1e-9
            result.pdfCond(i, j) = NaN;
        else
            result.pdfCond(i, j) = real(local_ilt(gpdf, t(j) - xi, nterms));
        end
    end
end
result.pdfUncond = NaN(1, numel(t));
result.cdfUncond = NaN(1, numel(t));
if hasPdf
    for j = 1:numel(t)
        if t(j) <= 0
            result.pdfUncond(j) = 0;
            result.cdfUncond(j) = 0;
            continue
        end
        result.pdfUncond(j) = real(local_ilt(lstUncond, t(j), nterms));
        result.cdfUncond(j) = real(local_ilt(@(u) lstUncond(u) / u, t(j), nterms));
    end
end
end

function v = local_lst_cond(s, x, rho, denom, excess)
% E[exp(-s V(x))] = (1-rho)/D(s,x), the excess flag returning instead the
% transform of V(x)-x, and D carrying its dominant exponential separately so
% that neither factor overflows at large s
[val, scale] = denom(s, x);
v = (1 - rho) * exp((excess * s - scale) * x) ./ val;
v(x == 0) = 1;
end

function v = local_lst_uncond(s, rho, denom, yq, wq, bq)
if isempty(bq)
    line_error('qsys_mg1_ps', ['the service density is required to remove the conditioning, ' ...
        'pass it as the pdf option']);
end
v = sum(wq .* bq .* local_lst_cond(s, yq, rho, denom, 0));
end

function tau = local_root(s, lambda, bhat, rho)
% unique root of tau = s + lambda*(1-bhat(tau)) in the right half plane
tau = s;
maxit = max(200, ceil(3 * log(1e-15) / log(max(rho, 1e-3))));
for it = 1:maxit
    taunew = s + lambda * (1 - bhat(tau));
    if abs(taunew - tau) <= 1e-14 * max(1, abs(taunew))
        tau = taunew;
        return
    end
    tau = taunew;
end
line_error('qsys_mg1_ps', 'the dominant root iteration did not converge');
end

function [val, scale] = local_denom_ph(s, x, lambda, rho, db, nb)
% D(s,x) = exp(scale*x)*val for phase-type service, exactly, from the residues
% of f(tau;s), whose double pole at the origin cancels, and for a whole vector
% of service requirements at once
n = numel(nb);
dm = db - [0, nb];
P = conv([1, -(s + lambda)], db) + [0, 0, lambda * nb];
A = conv([(1 - rho), s * rho, 0], db) - [0, conv([(1 - rho), s], lambda * dm)];
if norm(A(end-1:end)) > 1e-6 * max(1, norm(A))
    line_error('qsys_mg1_ps', 'the double pole at the origin did not cancel');
end
Ahat = A(1:n+1);
r = roots(P);
scale = max(real(r));
sep = abs(r - r.');
sep(1:numel(r)+1:end) = Inf;
if all(min(sep, [], 2) > 1e-7 * max(1, max(abs(r))))
    coef = polyval(Ahat, r) ./ polyval(polyder(P), r);
    val = reshape(exp(x(:) * (r - scale).') * coef, size(x));
else
    % repeated poles: use the companion realization of Ahat/P instead
    Pn = P / P(1);
    Ac = [-Pn(2:end); [eye(n) zeros(n, 1)]];
    e1 = [1; zeros(n, 1)];
    val = zeros(size(x));
    for k = 1:numel(x)
        val(k) = (Ahat / P(1)) * expm((Ac - scale * eye(n + 1)) * x(k)) * e1;
    end
end
end

function [val, scale] = local_denom_gen(s, x, lambda, rho, bhat, rootfun, shift, nterms)
% D(s,x) for a general service transform, by inverting f(tau;s) in tau on a
% contour placed just to the right of the dominant singularity
scale = real(rootfun(s)) + shift;
f = @(u) local_f(u + scale, s, lambda, rho, bhat);
val = zeros(size(x));
for k = 1:numel(x)
    val(k) = local_ilt(f, x(k), nterms);
end
end

function v = local_f(tau, s, lambda, rho, bhat)
bh = bhat(tau);
num = (1 - rho) * tau^2 - (1 - rho) * lambda * (1 - bh) * tau + s * rho * tau - s * lambda * (1 - bh);
v = num / (tau^2 * (tau - s - lambda * (1 - bh)));
end

function g = local_ilt(fun, t, nterms)
% Abate-Whitt Euler inversion, symmetrized so that complex-valued time
% functions are handled as well as real-valued ones
ne = floor((nterms - 1) / 2);
eta = [0.5, ones(1, ne), zeros(1, ne - 1), 2^(-ne)];
for k = 1:ne-1
    eta(2 * ne - k + 1) = eta(2 * ne - k + 2) + ...
        exp(gammaln(ne + 1) - ne * log(2) - gammaln(k + 1) - gammaln(ne - k + 1));
end
k = 0:2 * ne;
beta = ne * log(10) / 3 + 1i * pi * k;
eta = 10^(ne / 3) * (1 - mod(k, 2) * 2) .* eta;
g = 0;
for j = 1:numel(k)
    bj = beta(j) / t;
    g = g + 0.5 * eta(j) * (fun(bj) + fun(conj(bj)));
end
g = g / t;
end

function [y, w] = local_quad_nodes(ymax, npanel, ng)
% panelled Gauss-Legendre rule on [0,ymax], with the panels growing
% geometrically so that both ends of an exponentially decaying density are
% resolved. The rule is fixed, so it is identical in every codebase.
if nargin < 2
    npanel = 8;
end
if nargin < 3
    ng = 32;
end
% Golub-Welsch: nodes and weights from the Jacobi matrix of the Legendre family
kk = 1:ng-1;
bk = kk ./ sqrt(4 * kk.^2 - 1);
[V, D] = eig(diag(bk, 1) + diag(bk, -1));
[xg, idx] = sort(diag(D)');
wg = 2 * V(1, idx).^2;
edges = [0, ymax * 2.^(-npanel:0)];
y = zeros(1, (numel(edges) - 1) * ng);
w = zeros(1, numel(y));
for k = 1:numel(edges)-1
    a = edges(k);
    b = edges(k+1);
    y((k-1)*ng + (1:ng)) = 0.5 * (a + b) + 0.5 * (b - a) * xg;
    w((k-1)*ng + (1:ng)) = 0.5 * (b - a) * wg;
end
end

function m2v = local_second_moment(lst, meanref, meanscale)
% second moment from the transform curvature at the origin. The stencil is
% one-sided so that the transform is never sampled at negative arguments, where
% it need not converge, the step is scaled by the conditional mean but capped by
% the unconditional one so that it stays finite as x -> 0, and Richardson
% extrapolation over h and h/2 removes the leading truncation.
if meanref <= 0
    m2v = 0;
    return
end
h = min(1e-2 / meanref, 1 / meanscale);
m2v = (16 * local_d2_forward(lst, h / 2) - local_d2_forward(lst, h)) / 15;
end

function d = local_d2_forward(lst, h)
f = zeros(1, 6);
for j = 0:5
    f(j + 1) = real(lst(j * h));
end
d = (45 * f(1) - 154 * f(2) + 214 * f(3) - 156 * f(4) + 61 * f(5) - 10 * f(6)) / (12 * h^2);
end
