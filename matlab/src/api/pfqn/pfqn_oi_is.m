function [G, lG, Q] = pfqn_oi_is(N, mu, options)
% [G, LG, Q] = PFQN_OI_IS(N, MU, OPTIONS)
%
% Importance-sampling (IS) estimate of the normalizing constant of a closed
% two-station order-independent (OI) tandem. This is PFQN_PAS_IS with an EMPTY
% swap graph: it does not restrict the sampled orderings to a pass-and-swap
% communicating class but samples over ALL microstates (every ordering of the
% present jobs is feasible), thereby recovering the plain OI normalizing
% constant of PFQN_NCOI (which it estimates by Monte Carlo rather than exact
% balanced-fairness convolution).
%
% Model. Two OI stations (1 = upstream, 2 = downstream) hold all N jobs (no
% delay). The ordered-state chain is irreducible, so the communicating class is
% the full set of orderings D = all permutations of the job multiset. Writing
% Phi_m for the balanced-fairness balance function of station m,
%   G = sum_{c in D} sum_{k=0}^{ell} Phi_1(c_{1..k}) Phi_2(c_{ell..k+1}),
% with the ordered product Phi_m(q) = prod_{p=1}^{|q|} 1/mu_m(n(q_{1..p})),
% n(.) the per-class COUNT vector of the prefix (Comte-Dorsman product form:
% pi(c) = (1/G) prod_j 1/mu(c_1..c_j) prod_s sigma_s^{N_s-n_s}/(N_s-n_s)!).
% The argument is the prefix MULTISET, not its support: OI property P1 asks
% only that mu be permutation-invariant, i.e. a function of the counts. The
% two coincide for a compatibility rate that reads only which classes are
% present, and differ for any count-dependent OI rate -- an INF station
% (mu(n) = sum_r n_r sigma_r) above all, which is exactly what
% SOLVER_NC_PAS_IS_ANALYZER feeds in as station 2 of a Delay + OI cycle.
%
% Auto-normalized IS. Orderings c are drawn by placing, at each step, a
% uniformly random present class (no placement constraint); the draw probability
% p(c) is the product of the reciprocal branching factors. Then
%   G[xi] = E_{C~p}[ (sum_k xi(C,k) Phi_1(C_{1..k}) Phi_2(C_{ell..k+1})) / p(C) ],
% and E[xi] = G[xi]/G[1] reuses the same samples (auto-normalized IS). Taking
% xi = number of class-r jobs in the prefix yields the class-r mean queue length
% at station 1.
%
% Parameters:
%   N       - (1 x R) closed population vector (the macrostate), finite.
%   mu      - cell {1 x 2} of function handles. mu{m}(n) returns the total OI
%             rank rate of station m for the per-class occupancy (count) vector
%             n (1 x R). Permutation-invariant (OI property P1), but NOT
%             necessarily a function of supp(n) alone. This is the
%             svcRateFun stored on an OI station.
%   options - solver options (optional). Fields used:
%               .samples  number of IS samples (default 1e4);
%               .seed     RNG seed for reproducibility (optional);
%               .verbose  print progress (default false);
%               .qlen     estimate the queue lengths too (default true). False
%                         estimates ONLY G: the prefix-count coefficients are
%                         neither allocated nor accumulated and Q comes back
%                         zero. The ordering is drawn from the same stream
%                         either way, so G is unchanged to the last bit.
%
% Returns:
%   G  - IS estimate of the OI normalizing constant (== PFQN_NCOI with Z=0).
%   lG - log(G).
%   Q  - (2 x R) IS estimate of the mean per-class queue length; Q(1,:) at
%        station 1, Q(2,:) = N - Q(1,:) at station 2.
%
% Example (two OI queues, R=5, no swap graph):
%   mu1 = [1 2 0.5]; nb1 = {1,2,3,[1 3],[2 3]};
%   mu2 = [1 2];     nb2 = {1,2,[1 2],1,2};
%   rate = @(nb,muv,n) sum(muv(unique([nb{n>0}])));
%   mu = {@(n) rate(nb1,mu1,n), @(n) rate(nb2,mu2,n)};
%   [G,lG,Q] = pfqn_oi_is([1 1 1 3 3], mu, struct('samples',1e5));
%
% See also PFQN_PAS_IS, PFQN_NCOI, PFQN_OI_FNC.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3, options = struct(); end
if ~iscell(mu), mu = {mu}; end
if numel(mu) ~= 2
    line_error(mfilename, 'pfqn_oi_is models a two-station OI tandem: mu must have exactly two rate functions.');
end

R = numel(N);
N = round(N(:)');
if any(~isfinite(N))
    line_error(mfilename, 'pfqn_oi_is requires finite (closed) populations.');
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

for s = 1:nsamples
    % ---- draw an ordering c by placing a uniformly random present class -----
    x = N;                     % remaining per-class jobs to place
    c = zeros(1, ell);
    logp = 0;
    for ppos = 1:ell
        avail = find(x > 0);   % OI: every present class is placeable
        na = numel(avail);
        pick = avail(randi(na));
        c(ppos) = pick;
        logp = logp - log(na);
        x(pick) = x(pick) - 1;
    end
    p_c = exp(logp);

    % ---- prefix balance Phi_1 and per-class counts at each cut -------------
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
    % ---- reversed-suffix balance Phi_2 (station 2 = reversed suffix) -------
    Phi2cut = ones(1, ell + 1);
    occ2 = zeros(1, R); phi = 1;
    for k = ell:-1:1
        cls = c(k);
        occ2(cls) = occ2(cls) + 1;
        phi = phi / mu{2}(occ2);
        Phi2cut(k) = phi;
    end

    % ---- split-convolution sample value for every coefficient -------------
    sv = zeros(1, nCoef);
    for k = 0:ell
        if k >= ell
            w2 = 1;
        else
            w2 = Phi2cut(k + 1);
        end
        w = Phi1(k + 1) * w2;
        sv(1) = sv(1) + w;
        if wantQ && k > 0
            sv(2:end) = sv(2:end) + w * cnt1(k + 1, :);
        end
    end
    accum = accum + sv / p_c;

    if verbose && mod(s, max(1, floor(nsamples / 10))) == 0
        line_printf('\npfqn_oi_is: %d/%d samples', s, nsamples);
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
