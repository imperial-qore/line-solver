function [G, lG, Q] = pfqn_pas_is(N, mu, H, options)
% [G, LG, Q] = PFQN_PAS_IS(N, MU, H, OPTIONS)
%
% Importance-sampling (IS) estimate of the normalizing constant of a SINGLE
% communicating class of a cyclic two-station pass-and-swap (P&S) queueing 
% network with swap graph H. This is the Monte-Carlo 
% counterpart of the exact microstate convolution PFQN_PAS_NC: it estimates the
% same per-communicating-class constant G_C but scales to populations where the
% exact enumeration of the feasible orderings becomes expensive.
%
% Model. Two OI/P&S stations (1 = upstream, 2 = downstream of the cycle) hold
% all N jobs (no delay). With a non-empty swap graph the ordered-state chain is
% reducible; the recurrent communicating class is the set of splits of the
% orderings that are non-decreasing w.r.t. the placement partial order induced
% by H (Comte & Dorsman, 2021, arXiv:2009.12299). Writing D for that set of
% orderings and Phi_m for the balanced-fairness balance function of station m,
%   G_C = sum_{c in D} sum_{k=0}^{ell} Phi_1(c_{1..k}) Phi_2(c_{ell..k+1}),
% where c_{1..k} is the length-k prefix placed at station 1 and c_{ell..k+1} the
% reversed suffix placed at station 2. Along a fixed ordering q the balanced-
% fairness value is the ordered product of reciprocal rank rates,
%   Phi_m(q) = prod_{p=1}^{|q|} 1 / mu_m(n(q_{1..p})),   n(.) = prefix counts,
% evaluated at the per-class COUNT vector of each prefix (OI property P1 makes
% mu permutation-invariant, i.e. a function of the counts -- not of the support
% alone, which differs as soon as any class holds two or more jobs).
%
% Auto-normalized IS (notebook generator IS_3). Orderings c are drawn from D by
% placing, at each step, a uniformly random placement-order-minimal present
% class; the draw probability p(c) is the product of the reciprocal branching
% factors. Then, for any coefficient xi,
%   G_C[xi] = E_{C~p}[ (sum_k xi(C,k) Phi_1(C_{1..k}) Phi_2(C_{ell..k+1})) / p(C) ],
% and E[xi] = G_C[xi]/G_C[1] reuses the SAME samples for numerator and
% denominator (auto-normalized IS; A. Owen, MCM notes; convergence per
% Agapiou et al. 2017). Taking xi = number of class-r jobs in the prefix yields
% the mean queue length of class r at station 1.
%
% Parameters:
%   N       - (1 x R) closed population vector (the macrostate), finite.
%   mu      - cell {1 x 2} of function handles. mu{m}(n) returns the total OI
%             rank rate of station m given the per-class occupancy (count)
%             vector n (1 x R); OI, so it depends only on supp(n), i.e. the sum
%             of the capacities of the servers compatible with the present
%             classes. This is exactly the svcRateFun stored on an OI/PAS node.
%   H       - (R x R) swap-graph adjacency. An ordering is feasible iff it is
%             non-decreasing w.r.t. H: class a may not precede class b whenever
%             H(b,a) ~= 0. The empty/all-zero graph reduces D to all orderings
%             (pure OI); the estimate then targets the OI constant of PFQN_NCOI.
%   options - solver options (optional). Fields used:
%               .samples  number of IS samples (default 1e4);
%               .seed     RNG seed for reproducibility (optional);
%               .verbose  print progress (default false);
%               .qlen     estimate the queue lengths too (default true). False
%                         estimates ONLY G: the prefix-count coefficients are
%                         neither allocated nor accumulated and Q comes back
%                         zero. The ordering is drawn from the same stream
%                         either way, so G is unchanged to the last bit -- this
%                         is for the callers that want G(N - e_r) and read
%                         nothing else from it.
%
% Returns:
%   G  - IS estimate of the communicating-class normalizing constant G_C.
%   lG - log(G).
%   Q  - (2 x R) IS estimate of the mean per-class queue length; Q(1,:) at
%        station 1, Q(2,:) = N - Q(1,:) at station 2.
%
% Example (two OI/P&S queues, R=5, star swap graph):
%   mu1 = [1 2 0.5]; nb1 = {1,2,3,[1 3],[2 3]};
%   mu2 = [1 2];     nb2 = {1,2,[1 2],1,2};
%   rate = @(nb,mu,n) sum(mu(unique([nb{n>0}])));
%   mu = {@(n) rate(nb1,mu1,n), @(n) rate(nb2,mu2,n)};
%   H  = [0 0 0 1 0; 0 0 0 0 1; 0 0 0 1 1; 0 0 0 0 0; 0 0 0 0 0];
%   [G,lG,Q] = pfqn_pas_is([1 1 1 3 3], mu, H, struct('samples',1e5));
%
% See also PFQN_PAS_NC, PFQN_NCOI, PAS_PLACEMENT, PAS_SWAP2ORDER.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4, options = struct(); end
if nargin < 3, H = []; end
if ~iscell(mu), mu = {mu}; end
if numel(mu) ~= 2
    line_error(mfilename, 'pfqn_pas_is models a two-station pass-and-swap tandem: mu must have exactly two rate functions.');
end

R = numel(N);
N = round(N(:)');
if any(~isfinite(N))
    line_error(mfilename, 'pfqn_pas_is requires finite (closed) populations.');
end
if isempty(H), H = zeros(R, R); end
H = (H ~= 0);
if ~isequal(size(H), [R R])
    line_error(mfilename, 'H must be a (R x R) swap-graph adjacency matrix.');
end

nsamples = 1e4;
if isfield(options, 'samples') && ~isempty(options.samples)
    nsamples = round(options.samples);
end
if isfield(options, 'seed') && ~isempty(options.seed)
    rng(options.seed);
end
verbose = isfield(options, 'verbose') && ~isempty(options.verbose) && options.verbose;
wantQ = true;
if isfield(options, 'qlen') && ~isempty(options.qlen)
    wantQ = logical(options.qlen);
end

ell = sum(N);

% Placement-order logic (feasible orderings) is isolated in PAS_PLACEMENT:
% placeable(x) returns the classes drawable next given remaining counts x.
[~, placeable] = pas_placement(H);

if ell == 0
    G = 1; lG = 0; Q = zeros(2, R);
    return
end

if wantQ
    nCoef = R + 1;             % xi = [1, n_{1,1}, ..., n_{1,R}]
else
    nCoef = 1;                 % xi = [1] alone; G needs no prefix counts
end
accum = zeros(1, nCoef);
Nmat = N;                      % row template for feasibility masking

for s = 1:nsamples
    % ---- draw an ordering c from D (auto-normalized IS, generator IS_3) ----
    x = Nmat;                  % remaining per-class jobs to place
    c = zeros(1, ell);
    logp = 0;                  % log p(c) = -sum log(branching factor)
    for ppos = 1:ell
        avail = placeable(x);                  % placement-order-minimal present classes
        na = numel(avail);
        if na == 0
            line_error(mfilename, 'swap graph induces no feasible ordering (cyclic placement order).');
        end
        pick = avail(randi(na));
        c(ppos) = pick;
        logp = logp - log(na);
        x(pick) = x(pick) - 1;
    end
    p_c = exp(logp);

    % ---- split-convolution sample value for every coefficient at once ------
    % Precompute prefix balance Phi_1 up to each cut and suffix balance Phi_2.
    % Phi1cut(k+1) = Phi_1(c_{1..k}); Phi2cut(k+1) = Phi_2(reversed c_{k+1..ell})
    Phi1 = ones(1, ell + 1);
    if wantQ
        cnt1 = zeros(ell + 1, R);
    end
    phi = 1; occ = zeros(1, R);
    for k = 1:ell
        cls = c(k);
        occ(cls) = occ(cls) + 1;
        phi = phi / mu{1}(occ);
        Phi1(k + 1) = phi;
        if wantQ
            cnt1(k + 1, :) = occ;
        end
    end
    % Station 2 receives the reversed suffix c(ell), c(ell-1), ..., c(k+1). Scan
    % positions ell..1 accumulating support; Phi2cut(k) is the balance of the
    % reversed suffix c_{ell..k} (length ell-k+1), so at cut k station 2 holds
    % c_{ell..k+1} whose balance is Phi2cut(k+1) (and 1 for the empty suffix).
    Phi2cut = ones(1, ell + 1);
    occ2 = zeros(1, R); phi = 1;
    for k = ell:-1:1
        cls = c(k);
        occ2(cls) = occ2(cls) + 1;
        phi = phi / mu{2}(occ2);
        Phi2cut(k) = phi;
    end

    % sample_val(coef) = sum_{k=0}^{ell} xi_coef(k) * Phi1(k) * Phi2cut-at-k
    sv = zeros(1, nCoef);
    for k = 0:ell
        w = Phi1(k + 1) * Phi2cut_at(Phi2cut, k, ell);
        sv(1) = sv(1) + w;                     % xi = 1
        if wantQ && k > 0
            sv(2:end) = sv(2:end) + w * cnt1(k + 1, :);   % xi = n_{1,r}
        end
    end
    accum = accum + sv / p_c;

    if verbose && mod(s, max(1, floor(nsamples / 10))) == 0
        line_printf('\npfqn_pas_is: %d/%d samples', s, nsamples);
    end
end

est = accum / nsamples;
G = est(1);
lG = log(G);
Q = zeros(2, R);
if wantQ
    if G > 0
        Q(1, :) = est(2:end) / G;
    end
    Q(2, :) = N - Q(1, :);
end
end

function v = Phi2cut_at(Phi2cut, k, ell)
% Balance of the reversed suffix c_{ell..k+1} (length ell-k) at cut k.
if k >= ell
    v = 1;                     % empty suffix
else
    v = Phi2cut(k + 1);        % accumulated from position ell down to k+1
end
end
