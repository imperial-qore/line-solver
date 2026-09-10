function [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_nc_pas_is_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER,METHOD] = SOLVER_NC_PAS_IS_ANALYZER(SN, OPTIONS)
%
% Importance-sampling (IS) normalizing-constant analysis of a closed two-station
% pass-and-swap (P&S) tandem with a non-empty swap graph (Casale, Comte &
% Dorsman, 2026). With a genuine swap graph the ordered-state chain is reducible
% (Comte & Dorsman, 2021, arXiv:2009.12299); the recurrent communicating class
% carries a per-class product form pi(c) = Phi_1(c_1) Phi_2(c_2)/G_C, and its
% constant G_C is estimated by the auto-normalized IS routine PFQN_PAS_IS. This
% is the Monte-Carlo counterpart of SOLVER_NC_OI_ANALYZER for the case that the
% exact OI convolution does not apply (non-empty swap graph).
%
% Station 1 is the upstream P&S queue (prefix of the ordering), station 2 the
% downstream queue (reversed suffix); the swap graph is a class-level property
% read from the upstream station. Both stations must be OI/P&S with unit per-
% class visits (the balanced-fairness framework). Mean per-class queue lengths
% come directly from PFQN_PAS_IS (auto-normalized IS, same samples for numerator
% and denominator). Per-class throughput uses the balanced-fairness ratio
% X_r = G(N - e_r)/G(N) with common random numbers across the N and N-e_r runs
% for variance reduction; utilization, response time, and system time follow the
% same conventions as SOLVER_NC_OI_ANALYZER.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
Tstart = tic;
iter = 1;
method = 'is';         % importance sampling, specialized to OI/P&S stations

M = sn.nstations;
K = sn.nclasses;

% ---- the route's own premises (class switching, open chains, a rate lattice,
% a fork, two stations, a missing rank rate, unit visits) are decided by
% NC_OI_REFUSAL, which SolverNC's support gate asks too, so a pair the report
% offers is a pair this analyzer accepts ------------------------------------
pasReason = nc_oi_refusal(sn, 'pas');
if ~isempty(pasReason)
    line_error(mfilename, pasReason);
end
N = round(sn.njobs(:)');

% ---- classify stations and read swap graph ---------------------------------
svc = cell(M, 1);
swapG = cell(M, 1);
for ist = 1:M
    if sn.sched(ist) == SchedStrategy.INF
        % An infinite server IS an order-independent queue: every job is in
        % service, so mu(c) = sum_i sigma_{c_i} is permutation invariant and
        % satisfies P1-P3. Admitting it here is what lets the canonical
        % IS + OI cyclic model be sampled, which is otherwise rejected for
        % having a station that is not literally declared OI/PAS.
        infrates = sn.rates(ist, :);
        svc{ist} = @(c) sum(infrates(c));
        swapG{ist} = [];
        continue
    end
    if sn.sched(ist) ~= SchedStrategy.PAS && sn.sched(ist) ~= SchedStrategy.OI
        line_error(mfilename, 'solver_nc_pas_is requires both stations to be OI/PAS/INF (station %d is not).', ist);
    end
    ind = sn.stationToNode(ist);
    if ind < 1 || ind > numel(sn.nodeparam) || ~isstruct(sn.nodeparam{ind})
        line_error(mfilename, 'station %d has no OI/PAS node parameters.', ist);
    end
    svc{ist} = sn.nodeparam{ind}.svcRateFun;
    if isempty(svc{ist})
        line_error(mfilename, 'OI/PAS station %d has no service rate function; set it via setService(@(c) ...).', ist);
    end
    if isfield(sn.nodeparam{ind}, 'swapGraph')
        swapG{ist} = sn.nodeparam{ind}.swapGraph;
    else
        swapG{ist} = [];
    end
end

% derive the global placement-order DAG H (recurrent communicating class) from
% the P&S dynamics via pas_swap2order; see _kb/06-solver-catalog.md (NC section)
G1 = swapG{1}; if isempty(G1), G1 = zeros(K, K); end
G2 = swapG{2}; if isempty(G2), G2 = zeros(K, K); end
H = pas_swap2order({G1, G2}, {svc{1}, svc{2}}, ones(1, K));

% ---- per-class visits (chain == class); require unit visits ----------------
V = zeros(M, K);
for r = 1:K
    c = find(sn.chains(:, r));
    vis = sn.visits{c};
    for ist = 1:M
        isf = sn.stationToStateful(ist);
        V(ist, r) = vis(isf, r);
    end
    vref = V(sn.refstat(r), r);
    if vref > 0
        V(:, r) = V(:, r) / vref;
    end
end
% ---- OI rank-rate handles on a per-class count vector ----------------------
% svcRateFun(c) takes an ordered microstate list; for an OI station it is
% permutation-invariant, so evaluate it on a canonical microstate for count n.
rate1 = @(n) svc{1}(pas_is_microstate(n));
rate2 = @(n) svc{2}(pas_is_microstate(n));
mu = {rate1, rate2};

% ---- IS options; fix a base seed so the N and N-e_r runs share randoms -----
isopt = options;
if ~isfield(isopt, 'samples') || isempty(isopt.samples)
    if isfield(options, 'iter_max') && ~isempty(options.iter_max) && options.iter_max > 1
        isopt.samples = options.iter_max;
    else
        isopt.samples = 1e4;
    end
end
if ~isfield(isopt, 'seed') || isempty(isopt.seed)
    if isfield(options, 'seed') && ~isempty(options.seed)
        isopt.seed = options.seed;
    else
        isopt.seed = 23456;
    end
end

% ---- normalizing constant and mean queue lengths at population N -----------
% An EMPTY placement DAG means no swap ever reorders the queue, so the sampled
% orderings are unrestricted and the constant is the plain OI one. Call
% PFQN_OI_IS for that case rather than PFQN_PAS_IS with H=0: they compute the
% same quantity, but the OI kernel is the one that says so, and it does not
% carry the swap bookkeeping the degenerate graph would make dead weight.
isOIonly = isempty(H) || ~any(H(:));
if isOIonly
    method = 'is';
    isKernel = @(NN, mm, oo) pfqn_oi_is(NN, mm, oo);
else
    isKernel = @(NN, mm, oo) pfqn_pas_is(NN, mm, H, oo);
end
[G, lG, Qpas] = isKernel(N, mu, isopt);

Q = zeros(M, K);
Q(1, :) = Qpas(1, :);
Q(2, :) = Qpas(2, :);

% ---- per-class throughput X_r = G(N - e_r)/G(N) (common random numbers) ----
% Only the constant is read here, so these runs skip the prefix-count
% coefficients: same stream, same G, none of the queue-length bookkeeping.
isoptG = isopt;
isoptG.qlen = false;
X = zeros(1, K);
for r = 1:K
    if N(r) > 0
        er = zeros(1, K); er(r) = 1;
        Gr = isKernel(N - er, mu, isoptG);
        if G > 0
            X(r) = Gr / G;
        end
    end
end

% ---- throughput, utilization, response time --------------------------------
T = zeros(M, K);
U = zeros(M, K);
R = zeros(M, K);
for ist = 1:M
    for r = 1:K
        T(ist, r) = X(r) * V(ist, r);
    end
end
for ist = 1:M
    S = sn.nservers(ist);
    if ~isfinite(S) || S <= 0, S = 1; end
    for r = 1:K
        if N(r) > 0
            er = zeros(1, K); er(r) = 1;
            muR = mu{ist}(er);            % rank rate with only class r present
            if muR > 0
                U(ist, r) = T(ist, r) / muR / S;
            end
        end
    end
end
for ist = 1:M
    for r = 1:K
        if T(ist, r) > 0
            R(ist, r) = Q(ist, r) / T(ist, r);
        end
    end
end

C = zeros(1, K);
for r = 1:K
    if X(r) > 0
        C(r) = N(r) / X(r);
    end
end

runtime = toc(Tstart);
end

% ==========================================================================
function c = pas_is_microstate(n)
% Canonical ordered microstate holding n_r copies of class r (n a count
% vector). For an order-independent station the rank rate is invariant to the
% ordering, so this representative suffices to evaluate svcRateFun(c).
c = repelem(1:numel(n), round(n));
end
