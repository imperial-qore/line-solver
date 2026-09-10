function [res, hitprob, missprob, delayedprob, it, sn] = da_cacheqn_retrieval(sn, netfun, options)
% [RES,HITPROB,MISSPROB,DELAYEDPROB,IT,SN] = DA_CACHEQN_RETRIEVAL(SN, NETFUN, OPTIONS)
%
% Decomposition-aggregation driver for a CLOSED integrated cache-queueing model
% whose Cache node has a delayed-hit retrieval system (Cache.setRetrievalSystem).
% Adapts da_cacheqn: the cache is relabelled as a ClassSwitch (read -> hit / miss)
% and the finite-population delayed-hit coalescing is made to emerge from the
% closed AMVA by giving the retrieval (fetch) station a load-dependent COALESCING
% service rate. With k would-be-miss jobs at the fetch station spanning
% d(k) = n_eff*(1-(1-1/n_eff)^k) distinct uncached items (n_eff = nitems -
% totalcapacity), the single fetch server releases a whole coalesced batch per
% fetch, so the job-completion rate is mu(k) = (1/F)*k/d(k), i.e.
% lldscaling(fetch,k) = k/d(k). The distinct-fetch (backend) throughput then
% saturates at 1/F and the delayed-hit fraction = would-be-miss - fetcher/X
% emerges from the finite population, rather than from an open-arrival closed form.
%
% Returns the netsolve result RES, and per-cache hitprob/missprob/delayedprob
% (1 x nclasses on the read class), the iteration count IT and the mutated SN.
%
% LIMITATIONS (FURTHER WORK NEEDED - EXPERIMENTAL):
%  - hitprob/missprob are the true cache probabilities (hit = P(item cached),
%    miss = 1 - hit) and are accurate (uniform: exact m/n; validated vs LDES).
%    delayedprob is returned as 0: the finite-population delayed-hit FRACTION is
%    NOT recovered analytically (it folds into miss). Deriving the exact
%    closed-population coalescing split remains an open problem.
%  - The coalescing THROUGHPUT benefit is captured only in DIRECTION and is
%    UNDERSTATED: vs LDES it recovers the right sign of the throughput gain but a
%    smaller magnitude (e.g. LN ~+21% where LDES shows ~+31%). The load-dependent
%    lldscaling(k)=k/d(k) closure is calibrated for small closed populations (the
%    LN sublayer regime) and drifts for large N.
%  - Absolute LCQ throughput may carry the pre-existing LN(MVA) approximation
%    error for backend-bottleneck models (present with or without retrieval).
%  - Single fetch station (single-backend) only.
%  - NO CLOSED-RETRIEVAL EXAMPLE ships in the suite: every retrieval_* example is
%    OPEN. On a hand-built closed model this driver readily hits a reducible /
%    singular routing and a zero read-rate denominator, so a plain closed model
%    can return an empty or unreliable table. Treat the closed path as
%    experimental and validate any closed model against LDES before trusting it.
%    The counterpart is ported to Python (api/da/cacheqn_retrieval.py), the
%    JAR (Da_cacheqn_retrieval) and C++ (cpp/include/line/api/da/
%    da_cacheqn_retrieval.h), which share the same caveats. Both C++ solver-level
%    analyzers are now in as well: solver_nc_cacheqn_retrieval.h and
%    solver_mva_cacheqn_retrieval.h, the latter wired into mvaDispatch branch 2
%    on the no-Source test. Verified 2026-07-29 that the closed MVA path agrees
%    with the closed NC path to 1e-9 on the Delay->Cache->Fetch fixture, which is
%    the check that closed the SolverMVA refusal.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

I = sn.nnodes;
K = sn.nclasses;

statefulNodes = find(sn.isstateful)';
statefulNodesClasses = [];
for ind = statefulNodes
    statefulNodesClasses(end+1:end+K) = ((ind-1)*K+1):(ind*K);
end

caches = find(sn.nodetype == NodeType.Cache);
if numel(caches) ~= 1
    line_error(mfilename, 'da_cacheqn_retrieval requires exactly one Cache node.');
end
ci = caches(1);
ch = sn.nodeparam{ci};

% --- retrieval configuration ---
rk = keys(ch.retrievalSystemQueueIndices);
readClass = double(rk(1)) + 1;                 % 1-indexed read class
queueNodes = double(ch.retrievalSystemQueueIndices{rk(1)});
if numel(queueNodes) ~= 1
    line_error(mfilename, 'da_cacheqn_retrieval currently supports a single-station (single-backend) retrieval system.');
end
fetchNode = queueNodes(1);
fetchStation = sn.nodeToStation(fetchNode);
nitems = ch.nitems;
totcap = sum(ch.itemcap);
n_eff = max(1, nitems - totcap);               % distinct uncacheable items

% closed population
Npop = round(sum(sn.njobs(~isinf(sn.njobs))));

% mean fetch service F of the read class at the fetch station (per-item identical here)
rcls0 = ch.retrievalClasses(1, readClass);
F = 1 / sn.rates(fetchStation, rcls0);

% --- load-dependent COALESCING rate on the fetch station: lldscaling(k)=k/d(k) ---
alpha = ones(1, max(1, Npop));
for k = 1:Npop
    d = n_eff * (1 - (1 - 1/n_eff)^k);
    alpha(k) = k / d;
end
if isempty(sn.lldscaling)
    sn.lldscaling = ones(sn.nstations, Npop);
elseif size(sn.lldscaling,2) < Npop
    sn.lldscaling(:, end+1:Npop) = 1;
end
sn.lldscaling(fetchStation, 1:Npop) = alpha;

% relabel the cache as a class switch
sn.nodetype(ci) = NodeType.ClassSwitch;

hitClass = ch.hitclass;
missClass = ch.missclass;
retrievalClasses = ch.retrievalClasses;        % (item x class) -> retrieval class
pread = ch.pread{readClass};
pread = pread(:).' / sum(pread);

% --- fixed point on the read-class cache arrival rate ---
hitprob = zeros(1, K);
missprob = zeros(1, K);
delayedprob = zeros(1, K);
res = struct();

lambda0 = zeros(1, K);
lambda0(readClass) = 1;                          % seed
fpopts = options;
fpopts.config.da_norm = @(d) norm(d, 1);
[~, it] = da_fpi(@da_sweep, lambda0, fpopts);

    function [xnew, xref] = da_sweep(x, itnum) %#ok<INUSD>
        lambda = x;

        % isolated cache occupancy -> per-item uncached (would-be-miss) prob pi0_i
        [gamma, lambda_cache, ~] = da_cache_isolate(ch, lambda);
        [~, ~, ~, pi0] = cache_miss_fpi(gamma, ch.itemcap, lambda_cache);   % 1 x nitems
        pi0 = reshape(pi0, 1, []);
        nonhit = sum(pread .* pi0);              % aggregate would-be-miss prob
        hp = 1 - nonhit;

        % see _kb/09-ldes-and-cache.md (da_cacheqn_retrieval) for rationale
        r = readClass;
        sn.rtnodes((ci-1)*K+r, :) = 0;
        sn.rtnodes((ci-1)*K+r, (ci-1)*K+hitClass(r)) = hp;
        for i = 1:nitems
            rcls = retrievalClasses(i, r);
            if rcls > 0
                sn.rtnodes((ci-1)*K+r, (fetchNode-1)*K+rcls) = pread(i) * pi0(i);
                sn.rtnodes((fetchNode-1)*K+rcls, :) = 0;
                sn.rtnodes((fetchNode-1)*K+rcls, (ci-1)*K+missClass(r)) = 1;
                sn.rtnodes((ci-1)*K+rcls, :) = 0;   % drop the unused Cache->Retr edge
            end
        end
        sn.rt = dtmc_stochcomp(sn.rtnodes, statefulNodesClasses);

        [visits, nodevisits, sn] = sn_refresh_visits(sn, sn.chains, sn.rt, sn.rtnodes);
        sn.visits = visits;
        sn.nodevisits = nodevisits;

        res = netfun(sn);

        % throughput -> new read arrival rate at the cache
        nv = cellsum(nodevisits);
        c = find(sn.chains(:, r));
        inchain = find(sn.chains(c, :));
        refnode = sn.stationToNode(sn.refstat(r));
        if sn.refclass(c) > 0
            denom = nv(refnode, sn.refclass(c));
        else
            denom = nv(refnode, r);
        end
        Xr = sum(res.XN(inchain));
        lambda(r) = Xr * nv(ci, r) / denom;

        % see _kb/09-ldes-and-cache.md (da_cacheqn_retrieval) for rationale
        hitprob(r) = hp;
        missprob(r) = nonhit;
        delayedprob(r) = 0;

        xnew = lambda;
        xref = x;
    end
end
