function logNstates = ctmc_state_space_logsize(sn, options)
% LOGNSTATES = CTMC_STATE_SPACE_LOGSIZE(SN, OPTIONS)
%
% Worst-case log-size of the CTMC state space induced by SN. The estimate is
% the product of four factors, summed in log space:
%
%   1. job placements: stars-and-bars C(n_k+M-1,M-1) per class, over the
%      stations that do NOT keep an ordered buffer, with open classes
%      truncated at the cutoff;
%   2. buffer orderings: a station outside the share family keeps the CLASS
%      SEQUENCE of the jobs it holds, so with K>1 classes each occupancy
%      vector is as many states as its sequences. Bounded by the number of
%      sequences of length at most sum(n_k) over K symbols, once per such
%      station;
%   3. service phases: sn.phasessz raised to the number of jobs that can be
%      in service concurrently at each station;
%   4. routing state: one pointer over the outgoing links per (node,class)
%      using RROBIN or WRROBIN.
%
% This is the quantity fed to CTMC_MEMORY_GATE. It is exposed separately so
% that SolverAUTO can screen CTMC out of a candidate ranking without building
% the chain; see SolverCTMC.isStateSpaceTractable.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(options)
    options = SolverCTMC.defaultOptions();
end

M = sn.nstations;
K = sn.nclasses;
NK = sn.njobs;

cutoff = options.cutoff;
if isempty(cutoff) || all(isinf(cutoff(:)))
    % Same default the CTMC analyzer installs for open/mixed models.
    cutoff = ceil(6000^(1/(M*K)));
else
    cutoff = max(cutoff(isfinite(cutoff)));
end

shareSched = [SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.DPS, ...
    SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, ...
    SchedStrategy.GPSPRIO, SchedStrategy.LPS];

% ORDERED BUFFERS, computed EXACTLY and with POPULATION CONSERVED across them.
% A station outside the share family keeps the class SEQUENCE of the jobs it
% holds. Three things this gets right that a bound did not:
%  - placement and ordering are counted TOGETHER, so they cannot double count;
%  - a job placed at one ordered station is NOT available to another (charging
%    each station the full population priced a 2-station PAS model at
%    1957*1957 = 3829849 against a true 5040);
%  - CUTOFF TRUNCATES AN OPEN CLASS'S POPULATION IN THE NETWORK, exactly as the
%    plain stars-and-bars term has always treated it, so open classes are
%    conserved too. Treating them as per-station independent priced
%    prio_hol_open (3 ordered stations, cutoff 1) at 32768 and refused a model
%    that solves in a fraction of a second.
% A finite station capacity bounds the buffer TOTAL, sum_k m_k <= cap, not each
% class separately; per-class bounds come from sn.classcap.
if isfield(sn,'issignal') && ~isempty(sn.issignal)
    buffered = find(~logical(sn.issignal(:)'));
else
    buffered = 1:K;
end
ordIdx = [];
if numel(buffered) > 1
    for i = 1:M
        if sn.sched(i) == SchedStrategy.EXT || any(sn.sched(i) == shareSched)
            continue
        end
        ordIdx(end+1) = i; %#ok<AGROW>
    end
end
nOrd = numel(ordIdx);

logNstates = 0;
nkEff = zeros(1,K);
for k = 1:K
    if isinf(NK(k)), nkEff(k) = cutoff; else, nkEff(k) = NK(k); end
end

% A ZERO per-class capacity means the class is DISABLED at that station, so it
% never occupies a slot there and the placement term must spread it over the
% stations that admit it, not over all M. ld_whittle_bandwidth disables each of
% its three PS routes for the other two classes; counting all M=4 priced it at
% nchoosek(9,6)^3 = 592704 states, 7852 GB under the quadratic byte model, and
% the gate refused a model whose true space is 7^3 = 343. The ordered branch
% below already reads classcap through capsPer; only the placement terms were
% blind to it.
%
% ZERO VISITS IS THE SECOND WAY IN, see CTMC_NEVER_OCCUPIES below.
hasCcap = isfield(sn,'classcap') && ~isempty(sn.classcap);
never = ctmc_never_occupies(sn, M, K);
admits = @(i, k) ~(hasCcap && i <= size(sn.classcap,1) ...
    && k <= size(sn.classcap,2) && sn.classcap(i,k) == 0) && ~never(i,k);
admitting = @(k, stations) nnz(arrayfun(@(i) admits(i,k), stations));

% A class no station admits holds no job anywhere, so its population is zero
% rather than the cutoff. Leaving it at the cutoff made the class occupy ordered
% buffers it can never enter: test_classswitch_ctmc is a Source -> FCFS Queue ->
% Sink model whose second open class exists only as the switch target on the
% queue's departure arc, and pricing it as cutoff many jobs of BOTH classes in
% the FCFS buffer read 184755 states against the 11 the generator enumerates,
% predicted 17.62 GB, and had the gate refuse a model that solves in 0.05 s.
for k = 1:K
    if nkEff(k) > 0 && admitting(k, 1:M) == 0
        nkEff(k) = 0;
    end
end

if nOrd == 0
    for k = 1:K
        Mk = max(admitting(k, 1:M), 1);
        logNstates = logNstates + gammaln(1+nkEff(k)+Mk-1) - gammaln(1+Mk-1) - gammaln(1+nkEff(k));
    end
else
    remIdx = setdiff(1:M, ordIdx);
    mRem = numel(remIdx);
    for k = setdiff(1:K, buffered)
        Mk = max(admitting(k, 1:M), 1);
        logNstates = logNstates + gammaln(1+nkEff(k)+Mk-1) - gammaln(1+Mk-1) - gammaln(1+nkEff(k));
    end
    njb = floor(nkEff(buffered));
    capsPer = zeros(nOrd, numel(njb));
    capTot = inf(1, nOrd);
    for a = 1:nOrd
        i = ordIdx(a);
        for bi = 1:numel(buffered)
            c = njb(bi);
            if isfield(sn,'classcap') && ~isempty(sn.classcap) && i <= size(sn.classcap,1) ...
                    && buffered(bi) <= size(sn.classcap,2) && isfinite(sn.classcap(i,buffered(bi)))
                c = min(c, floor(sn.classcap(i,buffered(bi))));
            end
            if ~admits(i, buffered(bi))
                c = 0;
            end
            capsPer(a,bi) = c;
        end
        if isfield(sn,'cap') && ~isempty(sn.cap) && i <= numel(sn.cap) && isfinite(sn.cap(i))
            capTot(a) = floor(sn.cap(i));
        end
    end
    if prod(njb+1) <= 1e6
        logNstates = logNstates + ctmc_logOrderedJoint(capsPer, capTot, njb, mRem);
    else
        Kb = numel(njb); T = sum(njb); logKb = log(Kb);
        logNstates = logNstates + nOrd * ((T+1)*logKb - log(Kb-1) + log1p(-exp(-(T+1)*logKb)));
        for k = buffered
            Mk = admitting(k, remIdx);
            if Mk >= 1
                logNstates = logNstates + gammaln(1+nkEff(k)+Mk-1) - gammaln(1+Mk-1) - gammaln(1+nkEff(k));
            end
        end
    end
end

if isfield(sn,'phasessz') && ~isempty(sn.phasessz)
    for i = 1:min(M, size(sn.phasessz,1))
        for k = 1:min(K, size(sn.phasessz,2))
            p = sn.phasessz(i,k);
            if ~isfinite(p) || p <= 1
                continue
            end
            if sn.sched(i) == SchedStrategy.EXT
                m = 1;
            elseif any(sn.sched(i) == shareSched)
                m = nkEff(k);
            else
                m = min(nkEff(k), sn.nservers(i));
            end
            if ~isfinite(m)
                m = nkEff(k);
            end
            logNstates = logNstates + gammaln(1+m+p-1) - gammaln(1+p-1) - gammaln(1+m);
        end
    end
end

if isfield(sn,'routing') && ~isempty(sn.routing) && isfield(sn,'connmatrix') && ~isempty(sn.connmatrix)
    for ind = 1:min(size(sn.routing,1), size(sn.connmatrix,1))
        nout = nnz(sn.connmatrix(ind,:));
        if nout <= 1
            continue
        end
        nrr = sum(sn.routing(ind,:) == RoutingStrategy.RROBIN | ...
                  sn.routing(ind,:) == RoutingStrategy.WRROBIN);
        if nrr > 0
            logNstates = logNstates + nrr * log(nout);
        end
    end
end
end

function never = ctmc_never_occupies(sn, M, K)
% NEVER = CTMC_NEVER_OCCUPIES(SN, M, K)
%
% (M,K) logical: station I provably never holds a class-K job.
%
% THE STATE SPACE GENERATOR ALREADY KNOWS THIS, and this file claims to bound
% what the generator enumerates, so it has to read the same field.
% State.spaceGeneratorNodes sets capacityc(ind,r) = 0 whenever the chain visit
% ratio sn.visits{c}(isf,r) is zero, because a class the routing never delivers
% to a station cannot be there. Only an SPN PLACE escapes, keeping its declared
% initial marking (State.initialOccupancy), so a Place is left alone here too.
%
% Reading only sn.classcap missed all of it: refreshCapacity zeroes a
% (station,class) pair when its SERVICE RATE is NaN, which says the class is
% disabled at the station, not that nothing routes to it as that class. A class
% that exists only as a switch target on some departure arc has a perfectly good
% service rate everywhere and a visit ratio of zero everywhere.
%
% VISIT ROWS ARE STATEFUL-NODE INDEXED, not station indexed (see
% sn_refresh_visits), so the station index is mapped through stationToStateful
% first. Anything this cannot resolve -- no visits, no mapping, no chain for the
% class -- leaves the pair admitted, which keeps the estimate an upper bound.

% A pair the chain never enters carries floating-point noise (~1e-17) rather
% than an exact 0; ctmc_ssg uses the same tolerance for the same reason.
VISIT_ZERO_TOL = 1e-12;

never = false(M, K);
if ~isfield(sn,'visits') || isempty(sn.visits) || ~isfield(sn,'chains') || isempty(sn.chains)
    return
end
if ~isfield(sn,'stationToStateful') || numel(sn.stationToStateful) < M
    return
end
hasNodeType = isfield(sn,'nodetype') && ~isempty(sn.nodetype) ...
    && isfield(sn,'stationToNode') && ~isempty(sn.stationToNode);
for i = 1:M
    if hasNodeType && i <= numel(sn.stationToNode)
        nd = sn.stationToNode(i);
        if nd >= 1 && nd <= numel(sn.nodetype) && sn.nodetype(nd) == NodeType.Place
            continue
        end
    end
    isf = sn.stationToStateful(i);
    for k = 1:K
        c = find(sn.chains(:,k), 1);
        if isempty(c) || c > numel(sn.visits) || isempty(sn.visits{c})
            continue
        end
        V = sn.visits{c};
        if isf < 1 || isf > size(V,1) || k > size(V,2)
            continue
        end
        v = V(isf,k);
        if isfinite(v) && abs(v) < VISIT_ZERO_TOL
            never(i,k) = true;
        end
    end
end
end

% ------------------------------------------------------------------------
function lv = ctmc_logOrderedJoint(capsPer, capTot, njb, mRem)
% Log count of (placement, ordering) configurations over ALL order-preserving
% stations at once, population conserved. State is the per-class jobs still to
% place, held in a linear array over the box 0..njb (no containers.Map).
Kb = numel(njb);
dims = njb + 1;
nstate = prod(dims);
L = -inf(1, nstate);
L(idxOf(njb, dims)) = 0;              % start: everything still to place
for a = 1:size(capsPer,1)
    Lnew = -inf(1, nstate);
    caps = capsPer(a,:); ct = capTot(a);
    for si = 1:nstate
        if ~isfinite(L(si)), continue; end
        rem = subOf(si, dims);
        av = min(caps, rem);
        mgrid = cell(1,Kb);
        for k = 1:Kb, mgrid{k} = 0:av(k); end
        [gg{1:Kb}] = ndgrid(mgrid{:});
        for gi = 1:numel(gg{1})
            m = zeros(1,Kb);
            for k = 1:Kb, m(k) = gg{k}(gi); end
            t = sum(m);
            if isfinite(ct) && t > ct, continue; end
            v = L(si) + gammaln(t+1) - sum(gammaln(m+1));
            di = idxOf(rem - m, dims);
            if isfinite(Lnew(di))
                mx = max(Lnew(di), v);
                Lnew(di) = mx + log(exp(Lnew(di)-mx) + exp(v-mx));
            else
                Lnew(di) = v;
            end
        end
    end
    L = Lnew;
end
terms = [];
for si = 1:nstate
    if ~isfinite(L(si)), continue; end
    rem = subOf(si, dims);
    v = L(si);
    if mRem >= 1
        v = v + sum(gammaln(rem+mRem) - gammaln(rem+1) - gammaln(mRem));
    elseif any(rem > 0)
        continue
    end
    terms(end+1) = v; %#ok<AGROW>
end
if isempty(terms), lv = -inf; return; end
top = max(terms);
lv = top + log(sum(exp(terms - top)));
end

function ix = idxOf(v, dims)
ix = 1; mult = 1;
for k = 1:numel(dims)
    ix = ix + v(k)*mult;
    mult = mult * dims(k);
end
end

function v = subOf(ix, dims)
v = zeros(1, numel(dims)); r = ix - 1;
for k = 1:numel(dims)
    v(k) = mod(r, dims(k));
    r = floor(r / dims(k));
end
end
