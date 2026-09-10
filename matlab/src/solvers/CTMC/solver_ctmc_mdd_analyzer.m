function [QN,UN,RN,TN,CN,XN,info] = solver_ctmc_mdd_analyzer(sn, options, model)
% [QN,UN,RN,TN,CN,XN,INFO] = SOLVER_CTMC_MDD_ANALYZER(SN, OPTIONS)
% [QN,UN,RN,TN,CN,XN,INFO] = SOLVER_CTMC_MDD_ANALYZER(SN, OPTIONS, MODEL)
% Stationary analysis of a closed single-class network, or of a stochastic
% Petri net, whose CTMC state space is held in a decision diagram and solved by
% level aggregation, after
% A.S. Miner, G. Ciardo, S. Donatelli, "Using the exact state space of a Markov
% model to compute approximate stationary measures", SIGMETRICS 2000.
%
% This is the 'mdd' method of SolverCTMC. It never forms the |S|-state
% generator: the reachable set is stored in an MDD and K coupled level-CTMCs
% are iterated to a fixed point, so the memory cost is O(sum_k |M_k|) rather
% than O(|S|). The saving grows with the number of stations, and is negative
% at K=3, where the diagram compresses nothing.
%
% -- Exactness
% The single approximation is Pr{i_k | alpha} = Pr{i_k | p}. It is EXACT on
% product-form networks (paper Sec. 5), which covers exponential service under
% any work-conserving discipline and general service at PS or IS stations
% (BCMP types 2 and 3). It is an approximation otherwise, notably phase-type
% service at FCFS or LCFS, where errors of a fraction of a percent on the mean
% queue lengths have been observed.
%
% -- Petri nets
% Passing MODEL routes a net holding Places and Transitions through SPN_MDD
% instead of MDD_DESCRIPTOR: the levels are then (place, class) pairs plus one
% phase level per phase-type mode, and the measures come back per place. The
% approximation is the same Eq. 5 as for a queueing network, and it is exact on
% a product-form net, which SolverNC's 'rec' method (SOLVER_NC_SPN_ANALYZER)
% solves exactly and far more cheaply -- the aggregation earns its place on the
% nets that have NO product form.
%
% A net carries no per-station service rate, so MDD_MCD returns only the level
% marginals. The token throughput is then assembled here from the mode rates
% and those marginals, under the SAME independence across levels that the
% aggregation already assumes: it is the method's own approximation applied
% once more, not a second one layered on top.
%
% -- Output
% QN,UN,RN,TN : per station, mean queue length, utilization, response time and
%               throughput
% CN,XN       : per class, system response time and throughput
% INFO        : struct with the diagram, the descriptor, the level sizes and
%               the iteration count
%
% See also: mdd_mcd, mdd_descriptor, mdd_ps, mdd_reachset, spn_mdd, SolverCTMC.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(options), options = struct(); end
if nargin < 3, model = []; end
if ~isfield(options, 'verbose'), options.verbose = 0; end

% The model-shape rules live in SOLVER_CTMC_MDD_SUPPORTS, which
% SolverCTMC.supportsModelMethod also asks: the analyzer must refuse exactly
% what the report refuses, and one predicate with two callers is what keeps
% the two from drifting apart. The routing chain P and the local-state
% encoding KIND are read off the same call rather than recomputed here.
[shapeOk, shapeReason, P, kind] = solver_ctmc_mdd_supports(sn, options);
if ~shapeOk
    line_error(mfilename, shapeReason);
end

if any(sn.nodetype == NodeType.Place)
    if isempty(model)
        line_error(mfilename, ['a stochastic Petri net is read from the model object, not from ' ...
            'the network structure; call SOLVER_CTMC_MDD_ANALYZER(SN, OPTIONS, MODEL)']);
    end
    [QN,UN,RN,TN,CN,XN,info] = i_spn(model, sn, options);
    return
end

M = sn.nstations;
R = sn.nclasses;
N = sn.njobs(1);

% ---- service laws and the station-to-station routing chain
mu = zeros(1, M); servers = zeros(1, M); proc = cell(1, M);
for i = 1:M
    mu(i) = sn.rates(i, 1);                 % 1/E[S] by LINE convention
    servers(i) = sn.nservers(i);
    nph = 1;
    if isfield(sn, 'phases') && ~isempty(sn.phases), nph = sn.phases(i, 1); end
    if nph > 1
        proc{i} = {full(sn.proc{i}{1}{1}), full(sn.proc{i}{1}{2})};
    end
end

% ---- the encoding the disciplines present admit, decided by the predicate
sched = cell(1, M);
for i = 1:M, sched{i} = sn.sched(i); end

switch kind
    case 'ps'
        desc = mdd_ps(mu, P, servers, N, struct('proc', {proc}));
    otherwise
        desc = mdd_descriptor(mu, P, servers, N, ...
            struct('proc', {proc}, 'sched', {sched}));
end

mdd = mdd_reachset(desc.domain, desc.init, desc.nextfun);
% The level iteration is an INNER numerical solve and needs a far tighter
% tolerance than the reported means: mdd_mcd verifies the population invariant
% at 1e-6, so a loose tolerance converges short of the fixed point and trips
% that guard. options.iter_tol is the solver-level fixed-point tolerance
% (default 1e-4, sized for AMVA outer loops) and must NOT be reused here; the
% level knobs are taken from options.config instead.
mcdopt = struct();
if isfield(options, 'config') && isstruct(options.config)
    if isfield(options.config, 'mdd_maxiter') && ~isempty(options.config.mdd_maxiter)
        mcdopt.maxiter = options.config.mdd_maxiter;
    end
    if isfield(options.config, 'mdd_tol') && ~isempty(options.config.mdd_tol)
        mcdopt.tol = options.config.mdd_tol;
    end
end
out = mdd_mcd(mdd.toStruct(), desc, mcdopt);

% ---- pack the analyzer contract
QN = zeros(M, R); UN = zeros(M, R); RN = zeros(M, R); TN = zeros(M, R);
for i = 1:M
    QN(i, 1) = out.QLen(i);
    UN(i, 1) = out.U(i);
    TN(i, 1) = out.X(i);
    if TN(i, 1) > 0
        RN(i, 1) = QN(i, 1) / TN(i, 1);     % Little's law at the station
    end
end

% system throughput at the reference station, per unit visit
ref = sn.refstat(1);
vis = 1;
if isfield(sn, 'visits') && ~isempty(sn.visits) && numel(sn.visits) >= 1 && ~isempty(sn.visits{1})
    v = sn.visits{1};
    % visits is indexed by STATEFUL node, refstat by station: they coincide only
    % when every stateful node is a station.
    isf = sn.stationToStateful(ref);
    if size(v, 1) >= isf && v(isf, 1) > 0, vis = v(isf, 1); end
end
XN = zeros(1, R); CN = zeros(1, R);
XN(1) = TN(ref, 1) / vis;
if XN(1) > 0, CN(1) = N / XN(1); end

info.mdd = mdd;
info.desc = desc;
info.levelSizes = out.levelSizes;
info.iters = out.iters;
info.numStates = mdd.cardinality();
info.encoding = kind;
% Structural certificate: when no diagram node is shared, conditioning on the
% node equals conditioning on the whole path and the single approximation is an
% identity, so the answer is EXACT without a reference solve. False means "not
% certified" rather than "approximate": a product-form model is exact however
% much its diagram shares.
info.noAggregation = out.noAggregation;
info.pathsPerLevel = out.pathsPerLevel;

if options.verbose
    line_printf(['\nCTMC-mdd: %d levels, |S| = %d held as %d level states (%.1fx), ' ...
        '%d fixed-point sweeps\n'], desc.K, info.numStates, sum(out.levelSizes), ...
        info.numStates / sum(out.levelSizes), out.iters);
    if out.noAggregation
        line_printf('  no node is shared, so this result is EXACT (certified structurally)\n');
    else
        line_printf('  max paths per node = %g, so exactness rests on product form\n', ...
            max(out.pathsPerLevel));
    end
end
end

% ------------------------------------------------------------------------
function [QN,UN,RN,TN,CN,XN,info] = i_spn(model, sn, options)
% The Petri-net route: SPN_MDD supplies the reachable set and the Kronecker
% descriptor, MDD_MCD aggregates, and the measures are read back per place.
[mdds, desc, spninfo] = spn_mdd(model, struct('verbose', options.verbose > 1));

mcdopt = struct();
if isfield(options, 'config') && isstruct(options.config)
    if isfield(options.config, 'mdd_maxiter') && ~isempty(options.config.mdd_maxiter)
        mcdopt.maxiter = options.config.mdd_maxiter;
    end
    if isfield(options.config, 'mdd_tol') && ~isempty(options.config.mdd_tol)
        mcdopt.tol = options.config.mdd_tol;
    end
end
out = mdd_mcd(mdds, desc, mcdopt);

M = sn.nstations; R = sn.nclasses; L = spninfo.nplacelevels;
QN = zeros(M, R); UN = zeros(M, R); RN = zeros(M, R); TN = zeros(M, R);

places = spninfo.places;
for pp = 1:numel(places)
    ist = sn.nodeToStation(places(pp));
    if ist < 1, continue; end
    for k = 1:R
        QN(ist, k) = out.QLen((pp - 1) * R + k);
        UN(ist, k) = QN(ist, k);        % a Place is an INF station: U = Q
    end
end

% Mode throughputs from the level marginals. P(m_l = v) is read off the level
% chain; the enabling degree of a mode is then treated as independent across
% its input levels, which is Eq. 5 of the paper applied once more rather than a
% fresh approximation.
pl = i_levelmarginals(out, mdds, L);
md = spninfo.modes;
for e = 1:numel(md)
    if md(e).nph > 1, continue; end     % no single rate; read the phase level
    x = md(e).D1(1) * i_meanservers(pl, md(e), L);
    % FIRING EVENTS, not tokens: SN_PN_AVG_RATES converts the Place rows to a
    % token rate afterwards, exactly as it does for the explicit CTMC path, and
    % weighting here as well would count a weighted arc twice.
    for l = 1:L
        if md(e).enab(l) > 0
            pp = floor((l - 1) / R) + 1; k = mod(l - 1, R) + 1;
            ist = sn.nodeToStation(places(pp));
            if ist >= 1
                TN(ist, k) = TN(ist, k) + x;
            end
        end
    end
end
for i = 1:M
    for k = 1:R
        if TN(i, k) > 0, RN(i, k) = QN(i, k) / TN(i, k); end
    end
end
% RN and CN below are provisional: SN_PN_AVG_RATES recomputes them on the token
% rate once the caller has converted TN.

XN = zeros(1, R); CN = zeros(1, R);
for k = 1:R
    ref = sn.refstat(k);
    if ref >= 1 && ref <= M, XN(k) = TN(ref, k); end
    Nk = sum(QN(:, k));
    if XN(k) > 0 && Nk > 0, CN(k) = Nk / XN(k); end
end

info = struct('mdd', spninfo.mdd, 'desc', desc, 'spn', spninfo, ...
    'levelSizes', out.levelSizes, 'iters', out.iters, 'marginal', {pl});
end

% ------------------------------------------------------------------------
function pl = i_levelmarginals(out, mdds, L)
% P(level l = v) for each place level, from the converged level chains.
% MDD_MCD works in the paper's orientation, paper level k <-> level K+1-k.
K = mdds.K;
pl = cell(1, L);
for l = 1:L
    k = K + 1 - l;
    p = zeros(1, mdds.domain(l));
    rows = out.Mrows{k}; pk = out.pik{k};
    for r = 1:size(rows, 1)
        p(rows(r, 2) + 1) = p(rows(r, 2) + 1) + pk(r);
    end
    s = sum(p);
    if s > 0, p = p / s; end
    pl{l} = p;
end
end

% ------------------------------------------------------------------------
function n = i_meanservers(pl, mde, L)
% E[min(enabling degree, servers)] under independence across the input levels.
lv = [];
for l = 1:L
    if mde.enab(l) > 0, lv(end + 1) = l; end %#ok<AGROW>
end
if isempty(lv), n = 1; return; end
kmax = Inf;
for l = lv, kmax = min(kmax, floor((numel(pl{l}) - 1) / mde.enab(l))); end
if isfinite(mde.srv), kmax = min(kmax, mde.srv); end
n = 0;
for k = 1:kmax
    ge = 1;                              % P(degree >= k) = prod_l P(m_l >= k*I_l)
    for l = lv
        thr = k * mde.enab(l) + 1;
        if thr > numel(pl{l}), ge = 0; break; end
        ge = ge * sum(pl{l}(thr:end));
    end
    n = n + ge;                          % E[min(deg,srv)] = sum_k P(min >= k)
end
end
