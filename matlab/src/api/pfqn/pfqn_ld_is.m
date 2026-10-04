function [G, lG] = pfqn_ld_is(L, N, Z, mu, options)
% [G, LG] = PFQN_LD_IS(L, N, Z, MU, OPTIONS)
%
% Importance-sampling (IS) estimate of the normalizing constant of a closed
% LOAD-DEPENDENT product-form queueing network. This is the load-dependent
% counterpart of PFQN_PAS_IS / PFQN_OI_IS: the same sample-an-ordering estimator,
% with the order-independent rank rate replaced by the load-dependent capacity.
%
% Identity. Every product-form station's balance function is the sum, over the
% orderings q of a given per-class count vector n, of an ordered product of a
% per-position factor:
%   F_i(n) = |n|!/prod_r(n_r!) * prod_r L(i,r)^{n_r} / prod_{k=1}^{|n|} mu_i(k)
%          = sum_{q: |q|=n} prod_{p=1}^{|n|} L(i,q_p) / mu_i(p),
% since the multiset has |n|!/prod_r(n_r!) orderings and each contributes the
% same ordered product. The delay (infinite-server) node is the special case
% mu_Z(k) = k, giving F_Z(n) = prod_r Z_r^{n_r}/n_r!; a single-server queue is
% mu_i(k) = 1; a c-server queue is mu_i(k) = min(k,c).
%
% Consequently, writing ell = sum(N) and letting a "cut vector" split an ordering
% c of all ell jobs into S contiguous segments (one per station),
%   G(N) = sum_{c} sum_{cuts} prod_{m=1}^{S} w_m(seg_m),
%   w_m(q) = prod_{p=1}^{|q|} L(m,q_p) / mu_m(p),
% because summing over the orderings of each segment independently reproduces
% prod_m F_m(n_m), and each count split (n_1,...,n_S) is realized exactly once.
%
% Estimator. An ordering c is drawn by placing, at each step, a uniformly random
% present class; p(c) is the product of the reciprocal branching factors. For the
% sampled c the inner sum over ALL cut vectors is computed exactly by the
% dynamic program
%   A_0(0) = 1,   A_m(k) = sum_{j=0}^{k} A_{m-1}(j) * w_m(c_{j+1..k}),
% so S(c) = A_S(ell) in O(S*ell^2) time (no cut enumeration). Then
%   G = E_{C~p}[ S(C) / p(C) ]
% is unbiased and is estimated by the sample mean.
%
% Parameters:
%   L       - (M x R) per-class service demands at the M queueing stations.
%   N       - (1 x R) closed population vector, finite.
%   Z       - (1 x R) aggregated think time (delay) demand; [] or zeros if none.
%   mu      - load-dependent capacities. Either an (M x ell) matrix with
%             mu(i,k) the capacity of station i holding k jobs, or a cell
%             {1 x M} of function handles mu{i}(k), or [] for the
%             load-independent case mu(i,k) = 1 (see PFQN_IS).
%   options - solver options (optional). Fields used:
%               .samples  number of IS samples (default 1e4);
%               .seed     RNG seed for reproducibility (optional).
%
% Returns:
%   G  - IS estimate of the normalizing constant G(N).
%   lG - log(G).
%
% Example (2 queues + delay, load-dependent):
%   L = [0.5 0.3; 0.2 0.4];  N = [3 2];  Z = [1 1];
%   mu = [1 2 2 2 2; 1 1 1 1 1];        % station 1 is a 2-server queue
%   [G,lG] = pfqn_ld_is(L, N, Z, mu, struct('samples',1e5));
%
% See also PFQN_IS, PFQN_OI_IS, PFQN_PAS_IS, PFQN_NCLD, PFQN_NC.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5, options = struct(); end
if nargin < 4, mu = []; end
if nargin < 3, Z = []; end

[M, R] = size(L);
N = round(N(:)');
if numel(N) ~= R
    line_error(mfilename, 'L must have as many columns as N has classes.');
end
if any(~isfinite(N))
    line_error(mfilename, 'pfqn_ld_is requires finite (closed) populations.');
end
if isempty(Z), Z = zeros(1, R); end
Z = Z(:)';
ell = sum(N);

nsamples = 1e4;
if isfield(options, 'samples') && ~isempty(options.samples)
    nsamples = round(options.samples);
end
if isfield(options, 'seed') && ~isempty(options.seed)
    rng(options.seed);
end

if ell == 0
    G = 1; lG = 0;
    return
end

% ---- assemble the station list: M queues, plus the delay as mu_Z(k)=k -------
% D(m,r) is the per-class demand of station m; B(m,k) its capacity at k jobs.
hasZ = any(Z > 0);
S = M + double(hasZ);
D = zeros(S, R);
B = ones(S, ell);
for i = 1:M
    D(i, :) = L(i, :);
    if isempty(mu)
        B(i, :) = 1;                        % load-independent single server
    elseif iscell(mu)
        for k = 1:ell, B(i, k) = mu{i}(k); end
    else
        B(i, 1:min(ell, size(mu, 2))) = mu(i, 1:min(ell, size(mu, 2)));
        if size(mu, 2) < ell                % extend with the last capacity
            B(i, size(mu,2)+1:ell) = mu(i, end);
        end
    end
end
if hasZ
    D(S, :) = Z;
    B(S, :) = 1:ell;                        % delay: mu_Z(k) = k
end
if any(B(:) <= 0)
    line_error(mfilename, 'load-dependent capacities must be strictly positive.');
end

acc = 0;
for s = 1:nsamples
    % ---- draw an ordering c (uniformly random present class at each step) ---
    x = N;
    c = zeros(1, ell);
    logp = 0;
    for ppos = 1:ell
        avail = find(x > 0);
        na = numel(avail);
        pick = avail(randi(na));
        c(ppos) = pick;
        logp = logp - log(na);
        x(pick) = x(pick) - 1;
    end

    % ---- exact inner sum over all cut vectors, by dynamic programming -------
    % A(k) = weight of assigning the first k jobs of c to the stations seen so
    % far; W(j+1,k) = w_m(c_{j+1..k}) accumulated incrementally over k.
    A = zeros(1, ell + 1);
    A(1) = 1;                                % A(k+1) indexes k jobs placed
    for m = 1:S
        Anew = zeros(1, ell + 1);
        for j = 0:ell
            if A(j + 1) == 0, continue, end
            w = 1;
            Anew(j + 1) = Anew(j + 1) + A(j + 1);        % empty segment
            for k = j+1:ell
                w = w * D(m, c(k)) / B(m, k - j);        % position within segment
                if w == 0, break, end
                Anew(k + 1) = Anew(k + 1) + A(j + 1) * w;
            end
        end
        A = Anew;
    end
    acc = acc + A(ell + 1) * exp(-logp);
end

G = acc / nsamples;
lG = log(G);
end
