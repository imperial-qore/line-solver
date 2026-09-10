function out = infer_variational(spec, options)
% OUT = INFER_VARIATIONAL(SPEC, OPTIONS)
%
% Variational inference for Markovian queueing networks, following
% I. Perez, G. Casale, "Variational Inference for Markovian Queueing
% Networks", Advances in Applied Probability 53(3), 2021.
%
% The network trajectory is reparameterised by the transition counts
% Y^eta, eta=(i,j,c), so that the station marginals decouple:
%
%   x_{i,c}(t) = x_{i,c}(0) + sum_{eta in In(i,c)} Y^eta(t)
%                           - sum_{eta in Out(i,c)} Y^eta(t)
%
% The variational family is a product of inhomogeneous pure-birth
% processes, one per transition, with rate nu^eta(t,y), times a product
% of Gamma densities over the unknown service rates. The state space is
% expanded by adding DELTA to every feasible rate, so that queue lengths
% may go negative and the approximating measure stays mutually absolutely
% continuous with the target (Sec. 3 of the paper); the original model is
% recovered as DELTA -> 0.
%
% Each iteration performs, per transition, a backward pass for the
% Lagrange multipliers r^eta (Eq. 14) with multiplicative jumps at the
% observation epochs, the rate update
%   nu^eta(t,y) = exp(E log Xi^eta(t,y)) r^eta(t,y+1)/r^eta(t,y),
% and a forward pass of the master equation for the marginal. The
% conjugate Gamma posteriors are then refreshed from the expected number
% of firings and the expected exposure time of each station-class pair.
%
% Expectations over the other transitions are taken on a deterministic
% Halton lattice mapped through the inverse marginal c.d.f., so the
% estimator carries no random-number stream and is reproducible across
% the MATLAB, Java, Python and C++ implementations.
%
% SPEC fields (plain arrays; (m,r) pairs are flattened column-major, so
% that pair (m,r) sits at index (r-1)*M+m):
%   arcs      (T x 3)   [i j c] transitions; i=0 external source, j=0 sink
%   x0        (M x R)   initial queue lengths
%   sched     (M x 1)   0=INF, 1=shared server (PS/FCFS), 2=external
%   nservers  (M x 1)   number of servers
%   routeprob (T x 1)   routing probability p^c_{i,j} of each transition
%   arcparam  (T x 1)   index in 1..P of the rate governing the arc, 0 if known
%   arcrate   (T x 1)   known rate for arcs with arcparam==0, NaN otherwise
%   alpha0    (P x 1)   Gamma prior shapes
%   beta0     (P x 1)   Gamma prior rates
%   obsTimes  (K x 1)   observation epochs
%   obsData   (K x M*R) observed queue lengths, NaN where not observed
%   obsRange  (M x R)   support size of the uniform contamination
%   epsilon   scalar    probability that a reading is faulty
%   capacity  (M x R)   optional upper bound on the queue length, Inf by
%                       default. In a closed network this is the chain
%                       population, and clamping the load there keeps the
%                       expanded state space from crediting a station with
%                       more jobs than the network holds
%
% OPTIONS fields: tmax, dt, ngrid, ymax, nsamples, iter_max, tol, delta,
% floor, rate_max, verbose.
%
% OUT fields: alpha, beta, rates (posterior mean rates), meanServiceTime,
% bound (ELBO per iteration), Y (T x G x ymax+1 marginals), nu, tgrid,
% qlen (G x M*R expected queue lengths), iter, converged, tailmass.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

if nargin < 2
    options = struct();
end
[spec, opt] = ivSetup(spec, options);

M = size(spec.x0, 1);
R = size(spec.x0, 2);
narcs = size(spec.arcs, 1);
P = numel(spec.alpha0);
src = spec.arcs(:,1);
dst = spec.arcs(:,2);
cls = spec.arcs(:,3);

% signed incidence of each transition on class counts and station totals
sgnClass = zeros(narcs, M*R);
sgnStat = zeros(narcs, M);
for e = 1:narcs
    if dst(e) > 0
        k = (cls(e)-1)*M + dst(e);
        sgnClass(e,k) = sgnClass(e,k) + 1;
        sgnStat(e,dst(e)) = sgnStat(e,dst(e)) + 1;
    end
    if src(e) > 0
        k = (cls(e)-1)*M + src(e);
        sgnClass(e,k) = sgnClass(e,k) - 1;
        sgnStat(e,src(e)) = sgnStat(e,src(e)) - 1;
    end
end

G = opt.ngrid;
dt = opt.dt;
ymax = opt.ymax;
S = opt.nsamples;
ny = ymax + 1;
yvec = (0:ymax)';
tgrid = (0:(G-1)) * dt;

x0v = spec.x0(:)';
x0s = sum(spec.x0, 2)';

K = numel(spec.obsTimes);
obsIdx = zeros(K,1);
for k = 1:K
    obsIdx(k) = min(G, max(1, round(spec.obsTimes(k)/dt) + 1));
end

capStat = zeros(1,M);
for m = 1:M
    capStat(m) = sum(spec.capacity(m + (0:(R-1))*M));
end

% origin discipline and server count of each transition
arcSched = zeros(narcs,1);
arcServers = ones(narcs,1);
for e = 1:narcs
    if src(e) > 0
        arcSched(e) = spec.sched(src(e));
        arcServers(e) = spec.nservers(src(e));
    else
        arcSched(e) = 2;
    end
end

alpha = spec.alpha0(:);
beta = spec.beta0(:);

Y = zeros(narcs, G, ny);
nu = zeros(narcs, G, ny);
slack = zeros(narcs, G, ny);
gexp = zeros(narcs, G, ny);
hexp = ones(narcs, G, ny);

% initial homogeneous rates, from the mean observed occupancy when known
for e = 1:narcs
    lam = ivRateMean(spec, alpha, beta, e);
    if src(e) > 0
        u0 = ivUps(opt.xbar((cls(e)-1)*M+src(e)), opt.xbars(src(e)), ...
            arcServers(e), arcSched(e), spec.capacity((cls(e)-1)*M+src(e)), ...
            capStat(src(e)));
    else
        u0 = 1;
    end
    nue = zeros(G, ny);
    nue(:,1:ymax) = max(opt.delta, lam * u0);
    nu(e,:,:) = nue;
    Y(e,:,:) = ivForward(nue, dt, opt);
end

bound = zeros(opt.iter_max, 1);
alphaTrace = zeros(P, opt.iter_max);
betaTrace = zeros(P, opt.iter_max);
converged = false;
iter = 0;

for it = 1:opt.iter_max
    iter = it;

    for e = 1:narcs
        % deterministic lattice samples of every transition count from the
        % current marginals; the own contribution is removed inside
        % ivRateMoments, so the conditional expectations are under Q_{-eta}
        Ys = ivSampleAll(Y, S, opt);
        [ge, he, obsw] = ivRateMoments(spec, opt, e, Ys, Y, sgnClass, sgnStat, ...
            x0v, x0s, capStat, arcSched, arcServers, alpha, beta, obsIdx, yvec, true);

        Ye = reshape(Y(e,:,:), G, ny);
        r = ivBackward(ge, he, reshape(slack(e,:,:), G, ny), Ye, obsIdx, obsw, dt, opt);

        % Eq. (15). A vanishing multiplier marks a count the future
        % observations rule out; the rate there is zero, which is exactly
        % what keeps the forward pass from placing mass on it.
        nue = zeros(G, ny);
        den = r(:,1:ymax);
        num = r(:,2:ny);
        ratio = zeros(G, ymax);
        pos = den > 0;
        ratio(pos) = num(pos) ./ den(pos);
        nue(:,1:ymax) = he(:,1:ymax) .* ratio;
        nue(~isfinite(nue) | nue < 0) = 0;
        sl = zeros(G, ny);
        if isfinite(opt.rate_max)
            over = nue > opt.rate_max;
            if any(over(:))
                sl(over) = max(opt.floor, Ye(over)) .* log(nue(over)/opt.rate_max);
                nue(over) = opt.rate_max;
            end
        end
        nu(e,:,:) = nue;
        slack(e,:,:) = sl;
        Y(e,:,:) = ivForward(nue, dt, opt);
    end

    % conjugate Gamma updates: the shape gains the expected number of
    % firings, the rate the expected exposure time of the station-class
    % pair that the parameter governs
    Ys = ivSampleAll(Y, S, opt);
    firings = zeros(P,1);
    exposure = zeros(P,1);
    seen = zeros(P, M*R);
    for e = 1:narcs
        p = spec.arcparam(e);
        if p == 0
            continue
        end
        % expected number of firings of the transition over the horizon,
        % taken from the marginal itself, which is exact, rather than by
        % quadrature of the intensity, which a near-deterministic marginal
        % makes inaccurate
        Ye = reshape(Y(e,:,:), G, ny);
        firings(p) = firings(p) + (Ye(G,:) - Ye(1,:)) * yvec;
        kclass = (cls(e)-1)*M + src(e);
        if seen(p, kclass) == 0
            seen(p, kclass) = 1;
            ue = zeros(G,1);
            for g = 1:G
                Yg = reshape(Ys(:,g,:), narcs, S);
                Xall = repmat(x0v', 1, S) + sgnClass' * Yg;
                Xsta = repmat(x0s', 1, S) + sgnStat' * Yg;
                ue(g) = mean(ivUps(Xall(kclass,:), Xsta(src(e),:), ...
                    spec.nservers(src(e)), spec.sched(src(e)), ...
                    spec.capacity(kclass), capStat(src(e))));
            end
            exposure(p) = exposure(p) + ivTrapz(ue, dt);
        end
    end
    alpha = spec.alpha0(:) + firings;
    beta = spec.beta0(:) + exposure;
    alphaTrace(:,it) = alpha;
    betaTrace(:,it) = beta;

    % the bound is evaluated at the state the iteration ended in, so the rate
    % moments are recomputed against the updated marginals rather than reused
    % from the sweep that produced them
    for e = 1:narcs
        [ge, he] = ivRateMoments(spec, opt, e, Ys, Y, sgnClass, sgnStat, ...
            x0v, x0s, capStat, arcSched, arcServers, alpha, beta, obsIdx, yvec, false);
        gexp(e,:,:) = ge;
        hexp(e,:,:) = he;
    end

    bound(it) = ivBound(spec, opt, alpha, beta, Y, nu, gexp, hexp, Ys, ...
        sgnClass, x0v, obsIdx, dt);

    if opt.verbose > 0
        line_printf('\ninfer_variational: iteration %d, lower bound %.6f, max rate %.3f', ...
            it, bound(it), max(nu(:)));
    end
    % The rate update solves a stationarity condition rather than
    % maximising the bound in a block, so the bound need not ascend;
    % convergence is judged on the bound AND on the rate posteriors.
    if it > 1
        crit = abs(bound(it) - bound(it-1)) / max(1, abs(bound(it-1)));
        if P > 0
            prev = alphaTrace(:,it-1) ./ betaTrace(:,it-1);
            crit = max(crit, max(abs(alpha./beta - prev) ./ max(1e-12, prev)));
        end
        if crit <= opt.tol
            converged = true;
            break
        end
    end
end

tailmass = 0;
for e = 1:narcs
    tailmass = max(tailmass, max(Y(e,:,ny)));
end
if tailmass > 1e-6
    line_warning(mfilename, 'Transition-count truncation ymax=%d carries mass %.3e, increase options.ymax.', ymax, tailmass);
end

qlen = repmat(x0v, G, 1);
for e = 1:narcs
    my = reshape(Y(e,:,:), G, ny) * yvec;
    qlen = qlen + my * sgnClass(e,:);
end

out = struct();
out.alpha = alpha;
out.beta = beta;
out.rates = alpha ./ beta;
out.meanServiceTime = beta ./ alpha;
out.bound = bound(1:iter);
out.alphaTrace = alphaTrace(:,1:iter);
out.betaTrace = betaTrace(:,1:iter);
out.Y = Y;
out.nu = nu;
out.tgrid = tgrid;
out.qlen = qlen;
out.iter = iter;
out.converged = converged;
out.tailmass = tailmass;
end

%% ------------------------------------------------------------------ %%

function [spec, opt] = ivSetup(spec, options)
% Validate the specification and fill in the default options.

req = {'arcs','x0','sched','nservers','routeprob','arcparam','arcrate', ...
    'alpha0','beta0','obsTimes','obsData','obsRange','epsilon'};
for i = 1:numel(req)
    if ~isfield(spec, req{i})
        line_error(mfilename, sprintf('spec.%s is required.', req{i}));
    end
end

M = size(spec.x0, 1);
R = size(spec.x0, 2);
narcs = size(spec.arcs, 1);
if size(spec.arcs,2) ~= 3
    line_error(mfilename, 'spec.arcs must have three columns [i j c].');
end
if any(spec.arcs(:,1) == 0 & spec.arcs(:,2) == 0)
    line_error(mfilename, 'A transition cannot be external at both ends.');
end
if numel(spec.sched) ~= M || numel(spec.nservers) ~= M
    line_error(mfilename, 'spec.sched and spec.nservers must have one entry per station.');
end
if numel(spec.routeprob) ~= narcs || numel(spec.arcparam) ~= narcs || numel(spec.arcrate) ~= narcs
    line_error(mfilename, 'spec.routeprob, spec.arcparam and spec.arcrate must have one entry per transition.');
end
if numel(spec.alpha0) ~= numel(spec.beta0)
    line_error(mfilename, 'spec.alpha0 and spec.beta0 must have the same length.');
end
for e = 1:narcs
    if spec.arcparam(e) == 0 && ~(spec.arcrate(e) > 0)
        line_error(mfilename, sprintf('Transition %d has no parameter and no positive known rate.', e));
    end
end
if size(spec.obsData,1) ~= numel(spec.obsTimes)
    line_error(mfilename, 'spec.obsData must have one row per observation epoch.');
end
if size(spec.obsData,2) ~= M*R
    line_error(mfilename, 'spec.obsData must have M*R columns.');
end

spec.obsTimes = spec.obsTimes(:);
spec.routeprob = spec.routeprob(:);
spec.arcparam = spec.arcparam(:);
spec.arcrate = spec.arcrate(:);
spec.sched = spec.sched(:);
spec.nservers = spec.nservers(:);
spec.alpha0 = spec.alpha0(:);
spec.beta0 = spec.beta0(:);
spec.obsRange = reshape(spec.obsRange, 1, M*R);
if ~isfield(spec, 'capacity') || isempty(spec.capacity)
    spec.capacity = Inf(1, M*R);
else
    spec.capacity = reshape(spec.capacity, 1, M*R);
end

opt = struct();
opt.verbose = 0;
opt.iter_max = 20;
opt.tol = 1e-3;
opt.nsamples = 200;
opt.delta = 1e-3;
opt.floor = 1e-4;
opt.rate_max = [];
opt.rate_cap_factor = 10;
opt.unifmax = 30;
opt.unif_tol = 1e-12;
opt.unif_max_terms = 2000;
opt.tmax = [];
opt.dt = [];
opt.ngrid = [];
opt.ymax = [];
fn = fieldnames(options);
for i = 1:numel(fn)
    opt.(fn{i}) = options.(fn{i});
end

if isempty(opt.tmax)
    if isempty(spec.obsTimes)
        line_error(mfilename, 'options.tmax is required when there are no observations.');
    end
    opt.tmax = max(spec.obsTimes);
end
if opt.tmax <= 0
    line_error(mfilename, 'options.tmax must be positive.');
end
if isempty(opt.ngrid) && isempty(opt.dt)
    opt.ngrid = 201;
end
if isempty(opt.ngrid)
    opt.ngrid = round(opt.tmax/opt.dt) + 1;
end
opt.ngrid = max(2, round(opt.ngrid));
opt.dt = opt.tmax / (opt.ngrid - 1);

% mean occupancy used to size the truncation and the initial rates
xbar = spec.x0(:)';
for k = 1:(M*R)
    col = spec.obsData(:,k);
    col = col(~isnan(col));
    if ~isempty(col)
        xbar(k) = mean(col);
    end
end
xbars = zeros(1,M);
for m = 1:M
    xbars(m) = sum(xbar(m + (0:(R-1))*M));
end
opt.xbar = xbar;
opt.xbars = xbars;

if isempty(opt.ymax)
    fmax = 0;
    for e = 1:narcs
        if spec.arcparam(e) > 0
            lam = spec.routeprob(e) * spec.alpha0(spec.arcparam(e)) / spec.beta0(spec.arcparam(e));
        else
            lam = spec.routeprob(e) * spec.arcrate(e);
        end
        i = spec.arcs(e,1);
        if i > 0
            u = ivUps(xbar((spec.arcs(e,3)-1)*M+i), xbars(i), spec.nservers(i), ...
                spec.sched(i), spec.capacity((spec.arcs(e,3)-1)*M+i));
        else
            u = 1;
        end
        fmax = max(fmax, lam * u * opt.tmax);
    end
    opt.ymax = max(20, ceil(2*fmax + 5*sqrt(max(1,fmax))));
end
opt.ymax = max(2, round(opt.ymax));

% Cap on the variational rates. A count that the observations rule out
% leaves a vanishing multiplier in the denominator of Eq. (15), so the
% ratio is unbounded at the boundary of the excluded region; the paper
% constrains the rate there and carries the constraint into the backward
% equation through the slack multiplier. The default admits a rate that
% would traverse the whole count range ten times over the horizon.
if isempty(opt.rate_max)
    opt.rate_max = opt.rate_cap_factor * opt.ymax / opt.tmax;
end
end

%% ------------------------------------------------------------------ %%

function [ge, he, obsw] = ivRateMoments(spec, opt, e, Ys, Y, sgnClass, sgnStat, ...
    x0v, x0s, capStat, arcSched, arcServers, alpha, beta, obsIdx, yvec, wantObs)
% Conditional rate moments of one transition, and its observation jumps.
% Returns E[Xi | Y^eta=y] and exp(E[log Xi | Y^eta=y]) on the time grid, both
% taken under Q with the transition's own contribution removed, plus the
% multiplicative jump each observation carries in the backward pass.
[narcs, G, ny] = size(Y);
S = size(Ys, 3);
M = size(spec.x0, 1);
src = spec.arcs(:,1);
cls = spec.arcs(:,3);
lam = ivRateMean(spec, alpha, beta, e);
loglam = ivRateLogMean(spec, alpha, beta, e);
if src(e) > 0
    kclass = (cls(e)-1)*M + src(e);
else
    kclass = 0;
end
ge = zeros(G, ny);
he = zeros(G, ny);
obsw = ones(numel(spec.obsTimes), ny);
sgnC = sgnClass;
sgnC(e,:) = 0;
sgnS = sgnStat;
sgnS(e,:) = 0;
for g = 1:G
    Yg = reshape(Ys(:,g,:), narcs, S);
    Aall = x0v' + sgnC' * Yg;
    if kclass > 0
        Asta = x0s' + sgnS' * Yg;
        xic = Aall(kclass,:) + sgnClass(e,kclass) * yvec;
        xis = Asta(src(e),:) + sgnStat(e,src(e)) * yvec;
        ups = ivUps(xic, xis, arcServers(e), arcSched(e), ...
            spec.capacity(kclass), capStat(src(e)));
    else
        ups = ones(ny, S);
    end
    % E[Xi] and exp(E[log Xi]) of the SAME rate Xi = delta + lam*Ups. Writing
    % the second as exp(E[log lam]) exp(E[log(Ups + delta/E[lam])]) keeps the
    % two consistent wherever Ups is deterministic, which is what stops the
    % backward equation from developing a gradient away from the
    % empty-station boundary.
    ge(g,:) = opt.delta + lam * mean(ups, 2)';
    he(g,:) = exp(loglam + mean(log(ups + opt.delta/lam), 2))';
    if wantObs
        kk = find(obsIdx == g);
        for q = 1:numel(kk)
            obsw(kk(q),:) = ivObsWeight(spec.obsData(kk(q),:), spec.obsRange, ...
                spec.epsilon, Aall, sgnClass(e,:), yvec, opt);
        end
    end
end
end

%% ------------------------------------------------------------------ %%

function u = ivUps(xic, xis, nservers, sched, cap, capstat)
% Load factor Upsilon of a transition leaving a station-class pair.
% sched: 0 = infinite server, 1 = shared server (PS/FCFS), 2 = external.
if nargin < 5, cap = Inf; end
if nargin < 6, capstat = Inf; end
if sched == 2
    u = ones(size(xic));
    return
end
xic = min(cap, max(0, xic));
if sched == 0
    u = xic;
    return
end
xis = min(capstat, max(0, xis));
u = zeros(size(xic));
pos = xis > 0;
u(pos) = xic(pos) ./ xis(pos) .* min(nservers, xis(pos));
end

%% ------------------------------------------------------------------ %%

function lam = ivRateMean(spec, alpha, beta, e)
% E[lambda_eta] under the current Gamma posterior.
p = spec.arcparam(e);
if p == 0
    lam = spec.routeprob(e) * spec.arcrate(e);
else
    lam = spec.routeprob(e) * alpha(p) / beta(p);
end
end

function lg = ivRateLogMean(spec, alpha, beta, e)
% E[log lambda_eta] under the current Gamma posterior.
p = spec.arcparam(e);
if p == 0
    lg = log(spec.routeprob(e) * spec.arcrate(e));
else
    lg = log(spec.routeprob(e)) + psi(alpha(p)) - log(beta(p));
end
end

%% ------------------------------------------------------------------ %%

function ys = ivSampleAll(Y, S, opt)
% Deterministic lattice samples of every transition count.
[narcs, G, ny] = size(Y);
ys = zeros(narcs, G, S);
for e = 1:narcs
    ys(e,:,:) = ivSample(reshape(Y(e,:,:), G, ny), S, e, opt);
end
end

function ys = ivSample(q, S, e, opt) %#ok<INUSD>
% Inverse-c.d.f. samples of a marginal on a Halton lattice. Each
% transition uses its own prime base, so the samples of distinct
% transitions are jointly equidistributed rather than comonotone.
[G, ny] = size(q);
base = ivPrime(e);
u = zeros(1,S);
for s = 1:S
    u(s) = ivRadicalInverse(s, base);
end
[us, perm] = sort(u);
ys = zeros(G, S);
for g = 1:G
    c = cumsum(q(g,:));
    if c(ny) > 0
        c = c / c(ny);
    end
    c(ny) = 1;
    idx = zeros(1,S);
    j = 1;
    for s = 1:S
        while j < ny && c(j) < us(s)
            j = j + 1;
        end
        idx(s) = j - 1;
    end
    ys(g,perm) = idx;
end
end

function p = ivPrime(k)
% k-th prime, k >= 1.
p = 2;
n = 0;
c = 1;
while n < k
    c = c + 1;
    isp = true;
    d = 2;
    while d*d <= c
        if mod(c,d) == 0
            isp = false;
            break
        end
        d = d + 1;
    end
    if isp
        n = n + 1;
        p = c;
    end
end
end

function r = ivRadicalInverse(i, base)
% Van der Corput radical inverse of i in the given base.
r = 0;
f = 1/base;
while i > 0
    r = r + f * mod(i, base);
    i = floor(i/base);
    f = f / base;
end
end

%% ------------------------------------------------------------------ %%

function w = ivObsWeight(obsRow, obsRange, epsilon, Aall, sgnE, yvec, opt)
% Multiplicative jump carried by an observation in the backward pass,
% w(y) = exp(E_{Q_{-eta}}[log p(o | x)]) as a function of the own count.
ny = numel(yvec);
S = size(Aall, 2); %#ok<NASGU>
acc = zeros(ny, S);
obs = find(~isnan(obsRow));
for q = 1:numel(obs)
    k = obs(q);
    x = Aall(k,:) + sgnE(k) * yvec;
    hit = (x == obsRow(k));
    feas = (x >= 0) & (x <= obsRange(k));
    p = hit * (1-epsilon) + (~hit & feas) * (epsilon / max(1, obsRange(k)));
    acc = acc + log(opt.floor + p);
end
w = exp(mean(acc, 2))';
end

%% ------------------------------------------------------------------ %%

function r = ivBackward(ge, he, sl, Ye, obsIdx, obsw, dt, opt)
% Backward pass for the Lagrange multipliers. Equation (14),
%   d r(t,y)/dt = r(t,y) E[Xi(t,y)] - r(t,y+1) exp(E log Xi(t,y)),
% is LINEAR in r, so on a grid cell with frozen coefficients it is the
% action of a matrix exponential. The generator
%   B(y,y) = -E[Xi],  B(y,y+1) = +exp(E log Xi)
% has non-positive row sums by Jensen, so exp(B dt) is sub-stochastic and
% uniformization evaluates it without the stiffness that an explicit rule
% suffers when exp(E log Xi) is orders of magnitude below E[Xi]. Only the
% ratios r(y+1)/r(y) are used downstream, so r is rescaled at every step.
[G, ny] = size(ge);
r = zeros(G, ny);
v = ones(1, ny);
kk = find(obsIdx == G);
for q = 1:numel(kk)
    v = v .* max(0, obsw(kk(q),:));
end
v = ivRescale(v);
r(G,:) = v;
for g = (G-1):-1:1
    gv = max(0, ge(g,:));
    hv = max(0, he(g,:) .* ivDamp(sl(g,:), Ye(g,:), opt));
    hv = min(hv, gv);
    Lam = max(gv);
    if Lam > 0
        ncell = max(1, ceil(Lam*dt/opt.unifmax));
        h = dt / ncell;
        pd = gv / Lam;
        pu = hv / Lam;
        for c = 1:ncell
            v = ivBackUniformize(v, pd, pu, Lam*h, opt);
        end
        v = ivRescale(v);
    end
    kk = find(obsIdx == g);
    for q = 1:numel(kk)
        v = v .* max(0, obsw(kk(q),:));
        v = ivRescale(v);
    end
    r(g,:) = v;
end
end

function v = ivBackUniformize(v0, pd, pu, lt, opt)
% One uniformization step of the backward sub-generator, with pd = E[Xi]
% and pu = exp(E log Xi) both divided by the uniformization constant.
ny = numel(v0);
w = exp(-lt);
v = w * v0;
u = v0;
cum = w;
n = 1;
while (1 - cum) > opt.unif_tol && n < opt.unif_max_terms
    un = u .* (1 - pd);
    un(1:ny-1) = un(1:ny-1) + u(2:ny) .* pu(1:ny-1);
    u = un;
    w = w * lt / n;
    v = v + w * u;
    cum = cum + w;
    n = n + 1;
end
end

function v = ivRescale(v0)
% Normalise by the largest entry; the multipliers enter only through
% their ratios, so any positive scaling is immaterial.
m = max(v0);
if m > 0 && isfinite(m)
    v = v0 / m;
else
    v = v0;
end
end

function damp = ivDamp(sl, Ye, opt)
% Slack multiplier of the rate cap; unity when the cap is inactive.
z = sl ./ max(opt.floor, Ye);
damp = (1 + z) ./ exp(z);
damp(sl == 0) = 1;
damp(~isfinite(damp) | damp < 0) = 0;
end

%% ------------------------------------------------------------------ %%

function q = ivForward(nue, dt, opt)
% Forward master equation of an inhomogeneous pure-birth process, solved
% cell by cell with uniformization so that the marginal stays a proper
% distribution at every grid point.
[G, ny] = size(nue);
q = zeros(G, ny);
v = zeros(1, ny);
v(1) = 1;
q(1,:) = v;
for g = 1:(G-1)
    rates = max(0, nue(g,:));
    Lam = max(rates);
    if Lam <= 0
        q(g+1,:) = v;
        continue
    end
    ncell = max(1, ceil(Lam*dt/opt.unifmax));
    h = dt / ncell;
    p = rates / Lam;
    for c = 1:ncell
        v = ivUniformize(v, p, Lam*h, opt);
    end
    v = max(0, v);
    s = sum(v);
    if s > 0
        v = v / s;
    end
    q(g+1,:) = v;
end
end

function v = ivUniformize(v0, p, lt, opt)
% One uniformization step of the pure-birth chain with jump
% probabilities p and Poisson parameter lt.
ny = numel(v0);
w = exp(-lt);
v = w * v0;
u = v0;
cum = w;
n = 1;
while (1 - cum) > opt.unif_tol && n < opt.unif_max_terms
    un = u .* (1 - p);
    un(2:ny) = un(2:ny) + u(1:ny-1) .* p(1:ny-1);
    u = un;
    w = w * lt / n;
    v = v + w * u;
    cum = cum + w;
    n = n + 1;
end
end

%% ------------------------------------------------------------------ %%

function v = ivTrapz(f, dt)
% Trapezoidal integral of a grid function.
f = f(:);
n = numel(f);
if n < 2
    v = 0;
    return
end
v = dt * (sum(f) - 0.5*f(1) - 0.5*f(n));
end

function kl = ivKLGamma(a, b, a0, b0)
% KL(Gamma(a,b) || Gamma(a0,b0)) with rate parameterisation.
kl = (a - a0)*psi(a) - gammaln(a) + gammaln(a0) + a0*(log(b) - log(b0)) + a*(b0 - b)/b;
end

%% ------------------------------------------------------------------ %%

function b = ivBound(spec, opt, alpha, beta, Y, nu, gexp, hexp, Ys, sgnClass, x0v, obsIdx, dt)
% Evidence lower bound: path term, observation term and the divergence of
% the rate posteriors from their priors.
[narcs, G, ny] = size(Y);
S = size(Ys, 3);
b = 0;
for e = 1:narcs
    Ye = reshape(Y(e,:,:), G, ny);
    nue = reshape(nu(e,:,:), G, ny);
    ge = reshape(gexp(e,:,:), G, ny);
    he = reshape(hexp(e,:,:), G, ny);
    term = nue - ge;
    pos = nue > 0;
    tmp = zeros(G, ny);
    tmp(pos) = nue(pos) .* log(nue(pos) ./ max(opt.floor, he(pos)));
    term = term - tmp;
    b = b + ivTrapz(sum(Ye .* term, 2), dt);
end
for k = 1:numel(spec.obsTimes)
    g = obsIdx(k);
    Yg = reshape(Ys(:,g,:), narcs, S);
    Xall = repmat(x0v', 1, S) + sgnClass' * Yg;
    acc = zeros(1, S);
    row = spec.obsData(k,:);
    obs = find(~isnan(row));
    for q = 1:numel(obs)
        j = obs(q);
        x = Xall(j,:);
        hit = (x == row(j));
        feas = (x >= 0) & (x <= spec.obsRange(j));
        p = hit * (1-spec.epsilon) + (~hit & feas) * (spec.epsilon / max(1, spec.obsRange(j)));
        acc = acc + log(opt.floor + p);
    end
    b = b + mean(acc);
end
for p = 1:numel(alpha)
    b = b - ivKLGamma(alpha(p), beta(p), spec.alpha0(p), spec.beta0(p));
end
end
