function sn = solver_ctmc_pools(sn)
% SN = SOLVER_CTMC_POOLS(SN)
%
% Rewrite every heterogeneous-server station (Queue.addServerType) into the
% form the CTMC state space enumerates exactly, following the LDES semantics.
%
% A class-r job in service at a pooled station occupies one server of ONE
% compatible pool t, and is served by that pool's law (setHeteroService), or by
% the station's own law for r when the pool declares none. The per-class
% service process is therefore replaced by a block-diagonal phase-type law with
% one block per compatible pool, taken in ascending pool order: a phase of the
% class-r server block identifies both the pool and the service phase. The
% pool bookkeeping, which the event handler State.afterEventStationPool and the
% enumerator State.fromMarginalPool read, is stored in
% sn.nodeparam{ind}.ctmcpool:
%   ntypes, count(t), compat(t,r)     pools, their sizes, class compatibility
%   policy                            HeteroSchedPolicy of the station
%   pools{r}                          compatible pools of r, ascending
%   off{r}(k), len{r}(k)              offset and length of block k of class r
%   alpha{r}{k}, exit{r}{k}, D0{r}{k} entry vector, exit rates, D0 of block k
%   fsfrate(t,r)                      1/mean of the law of r at pool t (FSF)
%   alfsorder                         pools by ascending number of classes (ALFS)
%   rotate, perms, varpos             ALIS/FAIRNESS pool order, see below
%
% ALIS and FAIRNESS keep a global pool order: a job picks the first pool of
% that order among those with a free compatible server and that pool moves to
% the back, but only when it had more than one candidate. When some class has
% two or more compatible pools the order is part of the state: it is the index
% into ctmcpool.perms held in the shared trailing local variable (nvars column
% 2*R+1), at position ctmcpool.varpos of the local-variable block.
%
% The rewrite is idempotent (a station already carrying ctmcpool is skipped),
% so a struct returned by solver_ctmc can be handed back to it. It also adds
% the PHASE synchronizations a pool law with more than one phase needs, and
% rebuilds sn.state for the pooled stations in the new layout.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isfield(sn,'nodeparam') || isempty(sn.nodeparam)
    return
end
R = sn.nclasses;
snOld = sn;
changed = false(1, sn.nnodes);
for ind = 1:sn.nnodes
    if ~sn.isstation(ind) || numel(sn.nodeparam) < ind || isempty(sn.nodeparam{ind}) ...
            || ~isstruct(sn.nodeparam{ind}) || ~isfield(sn.nodeparam{ind},'nservertypes') ...
            || sn.nodeparam{ind}.nservertypes <= 0 || isfield(sn.nodeparam{ind},'ctmcpool')
        continue
    end
    ist = sn.nodeToStation(ind);
    % PAS/OI stations model heterogeneous compatible servers through the OI
    % rank rate (svcRateFun), not through server pools.
    if sn.sched(ist) == SchedStrategy.PAS || sn.sched(ist) == SchedStrategy.OI
        continue
    end
    sn = sub_pool_station(sn, ind, ist, R);
    changed(ind) = true;
end
if ~any(changed)
    return
end
sn.phasessz = max(sn.phases, ones(size(sn.phases)));
sn.phasessz(sn.nodeToStation(sn.nodetype == NodeType.Join),:) = sn.phases(sn.nodeToStation(sn.nodetype == NodeType.Join),:);
if isfield(sn,'markidx') && ~isempty(sn.markidx)
    sn.phasessz(sn.markidx > 1) = 1;
end
sn.phaseshift = [zeros(size(sn.phases,1),1), cumsum(sn.phasessz,2)];
local = sn.nnodes + 1;
for ind = find(changed)
    ist = sn.nodeToStation(ind);
    for r = 1:R
        if sn.phases(ist,r) > 1 && snOld.phases(ist,r) <= 1 && ~sub_has_phase_sync(sn.sync, ind, r)
            % refreshSync only emits a PHASE action where the station's own law has phases
            sn.sync{end+1,1} = struct('active',cell(1),'passive',cell(1));
            sn.sync{end,1}.active{1} = Event(EventType.PHASE, ind, r);
            sn.sync{end,1}.passive{1} = Event(EventType.LOCAL, local, r, 1.0);
        end
    end
    isf = sn.nodeToStateful(ind);
    if isfield(sn,'state') && numel(sn.state) >= isf && ~isempty(sn.state{isf})
        sn.state{isf} = sub_rebuild_state(snOld, sn, ind, ist, snOld.state{isf}(1,:));
    end
end
end

function sn = sub_pool_station(sn, ind, ist, R)
np = sn.nodeparam{ind};
name = sn.nodenames{ind};
T = np.nservertypes;
compat = logical(np.servercompat);
count = np.serverspertype(:)';
policy = HeteroSchedPolicy.ORDER;
if isfield(np,'heteroschedpolicy') && ~isempty(np.heteroschedpolicy)
    policy = np.heteroschedpolicy;
end
heteroproc = cell(T, R);
if isfield(np,'heteroproc') && ~isempty(np.heteroproc)
    heteroproc = np.heteroproc;
end
switch sn.sched(ist)
    case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRIO, SchedStrategy.LCFS, SchedStrategy.LCFSPRIO, SchedStrategy.SIRO}
    otherwise
        line_error(mfilename, sprintf(['SolverCTMC supports heterogeneous server pools under FCFS, HOL, FCFSPRIO, LCFS, LCFSPRIO ' ...
            'and SIRO; station ''%s'' uses %s.'], name, SchedStrategy.toText(sn.sched(ist))));
end
sub_refuse_features(sn, ind, ist, R, name);
if size(sn.nvars,2) >= 2*R+1 && sn.nvars(ind,2*R+1) > 0
    line_error(mfilename, sprintf(['Station ''%s'' combines heterogeneous server pools with another feature that owns ' ...
        'its trailing local variable (polling, BAS blocking, breakdown); SolverCTMC cannot represent both.'], name));
end
pools = cell(1,R);
off = cell(1,R);
len = cell(1,R);
alpha = cell(1,R);
exitr = cell(1,R);
D0b = cell(1,R);
fsfrate = zeros(T,R);
for r = 1:R
    base = sn.proc{ist}{r};
    baseOk = ~isempty(base) && iscell(base) && ~any(isnan(base{1}(:)));
    hasPoolLaw = false;
    for t = 1:T
        hasPoolLaw = hasPoolLaw || (compat(t,r) && ~isempty(heteroproc{t,r}));
    end
    if ~baseOk && ~hasPoolLaw
        continue % the class is not served here
    end
    pools{r} = find(compat(:,r))';
    if isempty(pools{r})
        line_error(mfilename, sprintf(['Station ''%s'' declares no server pool compatible with class ''%s'', so a job of ' ...
            'that class would wait forever.'], name, sn.classnames{r}));
    end
    if baseOk
        sub_require_ph(base, sn.procid(ist,r), name, sn.classnames{r}, 'its default service');
    end
    blocks = {};
    shift = 0;
    for k = 1:numel(pools{r})
        t = pools{r}(k);
        law = heteroproc{t,r};
        if isempty(law)
            if ~baseOk
                line_error(mfilename, sprintf(['Server pool ''%s'' of station ''%s'' accepts class ''%s'' but declares no ' ...
                    'law for it, and the station''s own service for the class is disabled.'], np.servertypenames{t}, name, sn.classnames{r}));
            end
            law = base;
        else
            sub_require_ph(law, [], name, sn.classnames{r}, sprintf('pool ''%s''', np.servertypenames{t}));
        end
        D0 = law{1};
        ex = -sum(D0,2);
        a = map_pie(law);
        a = a(:)' / sum(a);
        off{r}(k) = shift;
        len{r}(k) = size(D0,1);
        alpha{r}{k} = a;
        exitr{r}{k} = ex;
        D0b{r}{k} = D0;
        fsfrate(t,r) = 1 / map_mean(law);
        blocks{end+1} = D0; %#ok<AGROW>
        shift = shift + size(D0,1);
    end
    D0x = blkdiag(blocks{:});
    exx = -sum(D0x,2);
    piex = zeros(1, size(D0x,1));
    piex(1:len{r}(1)) = alpha{r}{1};
    sn.proc{ist}{r} = {D0x, exx * piex};
    sn.pie{ist}{r} = piex;
    sn.mu{ist}{r} = -diag(D0x);
    sn.phi{ist}{r} = exx ./ (-diag(D0x));
    sn.phases(ist,r) = size(D0x,1);
    sn.procid(ist,r) = ProcessType.PH;
    if isfield(sn,'isph') && ~isempty(sn.isph)
        sn.isph(ist,r) = true;
    end
end
ncls = sum(compat,2)';
[~, alfsorder] = sort(ncls, 'ascend'); % MATLAB sort is stable
rotate = (policy == HeteroSchedPolicy.ALIS || policy == HeteroSchedPolicy.FAIRNESS) ...
    && any(cellfun(@numel, pools) >= 2);
perms_t = zeros(0, T);
varpos = 0;
if rotate
    perms_t = sortrows(perms(1:T)); % identity first
    sn.nvars(ind, 2*R+1) = 1;
    varpos = sum(sn.nvars(ind, 1:(2*R))) + 1;
end
pool = struct();
pool.ntypes = T;
pool.count = count;
pool.compat = compat;
pool.policy = policy;
pool.pools = pools;
pool.off = off;
pool.len = len;
pool.alpha = alpha;
pool.exit = exitr;
pool.D0 = D0b;
pool.fsfrate = fsfrate;
pool.alfsorder = alfsorder;
pool.rotate = rotate;
pool.perms = perms_t;
pool.varpos = varpos;
sn.nodeparam{ind}.ctmcpool = pool;
% the pools are the server bank, as in LDES
sn.nservers(ist) = sum(count);
end

function sub_require_ph(law, procid, name, cname, what)
if ~isempty(procid) && (procid == ProcessType.MAP || procid == ProcessType.MMPP2 || procid == ProcessType.MMAP ...
        || procid == ProcessType.ME || procid == ProcessType.RAP)
    line_error(mfilename, sprintf(['SolverCTMC serves a heterogeneous server pool with a renewal phase-type law; class ''%s'' at ' ...
        'station ''%s'' has a %s law in %s.'], cname, name, ProcessType.toText(procid), what));
end
D0 = law{1};
D1 = law{2};
ex = -sum(D0,2);
a = map_pie(law);
a = a(:)' / max(sum(a), GlobalConstants.FineTol);
offd = D0 - diag(diag(D0));
if any(offd(:) < -GlobalConstants.FineTol) || any(ex < -GlobalConstants.FineTol) ...
        || norm(D1 - ex*a, 1) > GlobalConstants.CoarseTol * max(1, norm(D1,1))
    line_error(mfilename, sprintf(['SolverCTMC serves a heterogeneous server pool with a renewal phase-type law; class ''%s'' at ' ...
        'station ''%s'' has a correlated or matrix-exponential law in %s.'], cname, name, what));
end
end

function sub_refuse_features(sn, ind, ist, R, name)
why = '';
if ~isempty(sn.lldscaling) && size(sn.lldscaling,1) >= ist && any(sn.lldscaling(ist,:) ~= 1)
    why = 'load-dependent service';
elseif ~isempty(sn.cdscaling) && numel(sn.cdscaling) >= ist && ~isempty(sn.cdscaling{ist})
    why = 'class-dependent service';
elseif isfield(sn,'jdscaling') && ~isempty(sn.jdscaling) && numel(sn.jdscaling) >= ist && ~isempty(sn.jdscaling{ist})
    why = 'joint-dependent service';
elseif isfield(sn,'gdscaling') && ~isempty(sn.gdscaling)
    why = 'global dependence';
elseif isfield(sn,'hasbreakdown') && ~isempty(sn.hasbreakdown) && numel(sn.hasbreakdown) >= ind && sn.hasbreakdown(ind) == 1
    why = 'server breakdowns';
elseif isfield(sn,'retrialProc') && ~isempty(sn.retrialProc) && any(~cellfun(@isempty, sn.retrialProc(ist,:)))
    why = 'retrial';
elseif isfield(sn,'balkingStrategy') && ~isempty(sn.balkingStrategy) && any(sn.balkingStrategy(ist,:) > 0)
    why = 'balking';
elseif isfield(sn,'impatienceClass') && ~isempty(sn.impatienceClass) && any(sn.impatienceClass(ist,:) > 0)
    why = 'reneging';
elseif ~isempty(sn.isbasblocking) && numel(sn.isbasblocking) >= ind && sn.isbasblocking(ind) == 1
    why = 'BAS blocking';
elseif isfield(sn,'immfeed') && ~isempty(sn.immfeed) && any(sn.immfeed(ist,:))
    why = 'immediate feedback';
elseif isfield(sn,'replyblock') && ~isempty(sn.replyblock) && size(sn.replyblock,1) >= ind && any(sn.replyblock(ind,:) > 0)
    why = 'synchronous calls';
elseif isfield(sn.nodeparam{ind},'serverparallelism') && any(sn.nodeparam{ind}.serverparallelism > 1)
    why = 'server parallelism';
elseif isfield(sn,'issignal') && ~isempty(sn.issignal) && any(sn.issignal) && ~isempty(sn.rtnodes)
    for s = find(sn.issignal(:)')
        if any(sn.rtnodes(:, (ind-1)*R + s) > 0)
            why = 'signals';
            break
        end
    end
end
if ~isempty(why)
    line_error(mfilename, sprintf('SolverCTMC does not combine heterogeneous server pools with %s (station ''%s'').', why, name));
end
end

function tf = sub_has_phase_sync(sync, ind, r)
tf = false;
for a = 1:numel(sync)
    ev = sync{a}.active{1};
    if ev.event == EventType.PHASE && ev.node == ind && ev.class == r
        tf = true;
        return
    end
end
end

function st = sub_rebuild_state(snOld, sn, ind, ist, oldrow)
% Replay the jobs of the declared initial state as arrivals into an empty pooled
% station: the jobs in service first, by class, then the waiting jobs from the
% oldest to the newest. Each arrival takes the most likely entry phase and, under
% RAIS, the first candidate pool.
R = sn.nclasses;
Kold = snOld.phasessz(ist,:);
V = sum(snOld.nvars(ind,:));
Vnew = sum(sn.nvars(ind,:));
srvOld = oldrow((end-sum(Kold)-V+1):(end-V));
bufOld = oldrow(1:(end-sum(Kold)-V));
varOld = oldrow((end-V+1):end);
seq = [];
for r = 1:R
    seq = [seq, r*ones(1, sum(srvOld((snOld.phaseshift(ist,r)+1):(snOld.phaseshift(ist,r)+Kold(r)))))]; %#ok<AGROW>
end
if sn.sched(ist) == SchedStrategy.SIRO
    for r = 1:R
        seq = [seq, r*ones(1, bufOld(r))]; %#ok<AGROW>
    end
else
    seq = [seq, fliplr(bufOld(bufOld > 0))]; % rightmost is the oldest
end
pool = sn.nodeparam{ind}.ctmcpool;
if sn.sched(ist) == SchedStrategy.SIRO
    Wb = R;
else
    Wb = max(1, numel(seq));
end
varNew = zeros(1, Vnew);
varNew(1:V) = varOld;
if pool.rotate
    varNew(pool.varpos) = 1;
    if pool.varpos <= V
        line_error(mfilename, 'The pool order variable overlaps an existing local variable.');
    end
end
st = [zeros(1, Wb), zeros(1, sum(sn.phasessz(ist,:))), varNew];
for j = 1:numel(seq)
    [outs, ~, outp] = State.afterEventStationPool(sn, ind, ist, st, EventType.ARV, seq(j), false, [], R, Vnew, NaN, ...
        Inf(sn.nstations,1), Inf(sn.nstations,R));
    if isempty(outs)
        line_error(mfilename, sprintf('The initial state of station ''%s'' cannot be placed on its server pools.', sn.nodenames{ind}));
    end
    [~, best] = max(outp);
    st = outs(best,:);
end
end
