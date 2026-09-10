function [pi,SSq,arvRates,depRates,tranSysState,tranSync,sn,dlyRates,startRates,preemptRates,tranStartTag,tranPreemptTag]=solver_ssa(sn, init_state, options, eventCache)
% [PI,SSQ,ARVRATES,DEPRATES,TRANSYSSTATE,QN,DLYRATES,STARTRATES,PREEMPTRATES,TRANSTARTTAG,TRANPREEMPTTAG]=SOLVER_SSA(QN,OPTIONS)
%
% STARTRATES and PREEMPTRATES are the derived START/PREEMPT rates per unique
% state, laid out like ARVRATES/DEPRATES: how fast the transitions enabled in
% that state start a class-r service at a stateful node, or push a class-r job
% in service there back into the buffer.
%
% TRANSTARTTAG{e} and TRANPREEMPTTAG{e} are the tags of the transition that
% actually FIRED at step e, as rows [statefulIndex class]. The trace needs
% them explicitly: @SolverSSA/sample.m rebuilds events from sn.sync{tranSync(e)},
% which by construction carries no derived tag.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% by default the jobs are all initialized in the first valid state

% Impatience and server-count support checks. The rules live in
% SOLVER_CTMC_STATE_SUPPORTS, shared with SolverCTMC because both drive the
% same State.afterEventStation; SolverSSA.supportsModelMethod asks it too, so
% the report and this run cannot differ.
[stateOk, stateWhy] = solver_ctmc_state_supports(sn, 'SolverSSA');
if ~stateOk
    line_error(mfilename, stateWhy);
end

if ~isfield(options,'seed')
    options.seed = 23000;
end
% Handle parallel computing toolbox gracefully - get worker index
if isMATLABReleaseOlderThan("R2022b")
    % Use labindex for older MATLAB versions
    try
        if ~isempty(getCurrentTask())
            lab_idx = labindex();  %#ok<DLABINDEX>
        else
            lab_idx = 1;
        end
    catch
        line_warning(mfilename,'Parallel Computing Toolbox not available or not running in parallel mode. Using labindex = 1.');
        lab_idx = 1;
    end
else
    % Use spmdIndex for R2022b and newer (labindex is deprecated)
    try
        lab_idx = spmdIndex;
        if isempty(lab_idx) || lab_idx == 0
            lab_idx = 1;
        end
    catch
        line_warning(mfilename,'Parallel Computing Toolbox not available or not running in parallel mode. Using labindex = 1.');
        lab_idx = 1;
    end
end
Solver.resetRandomGeneratorSeed(options.seed + lab_idx - 1);

%% generate local state spaces
%nstations = sn.nstations;
nstateful = sn.nstateful;
%init_nserver = sn.nservers; % restore Inf at delay nodes
R = sn.nclasses;
N = sn.njobs';
nnodes = sn.nnodes;
sync = sn.sync;
gsync = sn.gsync;

line_debug('SSA solver starting: nstateful=%d, nclasses=%d, njobs=%s, samples=%d', nstateful, R, mat2str(N), options.samples);
csmask = sn.csmask;

cutoff = options.cutoff;
if isscalar(cutoff)
    cutoff = cutoff * ones(sn.nstations, sn.nclasses);
end

%%
Np = N';
capacityc = zeros(sn.nnodes, sn.nclasses);
original_classcap = sn.classcap; % preserve original classcap for class switching scenarios
for ind=1:sn.nnodes
    if sn.isstation(ind) % place jobs across stations
        ist = sn.nodeToStation(ind);
        %isf = sn.nodeToStateful(ind);
        for r=1:sn.nclasses %cut-off open classes to finite capacity
            c = find(sn.chains(:,r));
            % Check if visits is 0, but also preserve capacity for classes that can
            % receive jobs via class switching (indicated by non-zero original classcap)
            if isfield(sn,'fjclassmap') && ~isempty(sn.fjclassmap) && length(sn.fjclassmap) >= r && sn.fjclassmap(r) > 0
                % see _kb/06-solver-catalog.md for rationale (native fork-join)
                capacityc(ind,r) = original_classcap(ist,r);
            elseif ~isempty(sn.visits{c}) && sn.visits{c}(ist,r) == 0 && original_classcap(ist,r) == 0
                capacityc(ind,r) = 0;
            elseif ~isempty(sn.proc) && ~isempty(sn.proc{ist}{r}) && any(any(isnan(sn.proc{ist}{r}{1}))) && sn.nodetype(ind) ~= NodeType.Place % disabled (but not Place nodes)
                capacityc(ind,r) = 0;
            else
                if isinf(N(r))
                    capacityc(ind,r) =  min(cutoff(ist,r), sn.classcap(ist,r));
                else
                    % closed classes: enumerate up to the chain population, but never
                    % beyond the class capacity at this station (finite-buffer stations)
                    capacityc(ind,r) =  min(sum(sn.njobs(sn.chains(c,:))), sn.classcap(ist,r));
                end
            end
        end
        % never raise the station capacity above its configured total capacity
        capacity_sum = min(sum(capacityc(ind,:)), sn.cap(ist));
        if sn.sched(ist) == SchedStrategy.PAS
            % see _kb/06-solver-catalog.md for rationale (SSA PAS capacity)
            capacity_sum = sn.cap(ist);
        end
        if isinf(sn.nservers(ist))
            sn.nservers(ist) = capacity_sum;
        end
        sn.cap(ist,:) = capacity_sum;
        sn.classcap(ist,:) = capacityc(ind,:);
    end
end
% see _kb/06-solver-catalog.md for rationale (SSA G-network signal capacity)
if isfield(sn,'issignal') && ~isempty(sn.issignal) && any(sn.issignal)
    for ii = 1:sn.nstations
        if sn.sched(ii) ~= SchedStrategy.EXT
            sn.classcap(ii, sn.issignal(:)') = 0;
        end
    end
end

% see _kb/06-solver-catalog.md for rationale (SSA heterogeneous servers)
for ind = 1:sn.nnodes
    if sn.isstation(ind) && isfield(sn,'nodeparam') && numel(sn.nodeparam) >= ind ...
            && ~isempty(sn.nodeparam{ind}) && isstruct(sn.nodeparam{ind}) ...
            && isfield(sn.nodeparam{ind},'nservertypes') && sn.nodeparam{ind}.nservertypes > 0
        ist = sn.nodeToStation(ind);
        np = sn.nodeparam{ind};
        served = [];
        for r = 1:sn.nclasses
            if ~isempty(sn.proc{ist}{r}) && ~any(any(isnan(sn.proc{ist}{r}{1}))) && sn.rates(ist,r) > 0
                served(end+1) = r; %#ok<AGROW>
            end
        end
        if numel(served) > 1
            line_error(mfilename,'SolverSSA supports heterogeneous servers only for single-class stations. Use SolverJMT or SolverLDES for multi-class heterogeneous servers.');
        end
        if numel(served) == 1
            r = served;
            srvrates = [];
            for t = 1:np.nservertypes
                if np.servercompat(t,r) && np.heterorates(t,r) > 0
                    srvrates = [srvrates, repmat(np.heterorates(t,r), 1, np.serverspertype(t))]; %#ok<AGROW>
                end
            end
            c = numel(srvrates);
            mu_base = sn.rates(ist,r);
            if c > 0 && mu_base > 0
                if isempty(sn.lldscaling)
                    sn.lldscaling = ones(sn.nstations, max([c, sum(sn.njobs(isfinite(sn.njobs))), 1]));
                elseif size(sn.lldscaling,2) < c
                    sn.lldscaling(:, (size(sn.lldscaling,2)+1):c) = repmat(sn.lldscaling(:,end), 1, c-size(sn.lldscaling,2));
                end
                for n = 1:size(sn.lldscaling,2)
                    mun = sum(srvrates(1:min(n,c)));
                    sn.lldscaling(ist,n) = mun / (mu_base * min(n,c));
                end
            end
        end
    end
end

% see _kb/06-solver-catalog.md for rationale (SSA FCR)
fcrOn = isfield(sn,'nregions') && sn.nregions > 0;
if fcrOn
    Kfcr = sn.nclasses;
    fcrMembers   = cell(sn.nregions,1);
    fcrMemberMask= cell(sn.nregions,1);
    fcrClassCap  = cell(sn.nregions,1);
    fcrGlobalCap = inf(sn.nregions,1);
    fcrMemCap    = inf(sn.nregions,1);
    fcrSz        = cell(sn.nregions,1);
    fcrA         = cell(sn.nregions,1);
    fcrb         = cell(sn.nregions,1);
    for f = 1:sn.nregions
        Rmat = sn.region{f};                          % M x (K+1)
        % membership: any job-count cap OR the region memory budget set on the
        % station row (a memory-only region has all job-count entries at -1)
        memvecFCR = -ones(sn.nstations,1);
        if isfield(sn,'regionmaxmem') && numel(sn.regionmaxmem) >= f && ~isempty(sn.regionmaxmem{f})
            memvecFCR = sn.regionmaxmem{f}(:);
        end
        mask = sn_region_members(sn, f, Rmat, memvecFCR);
        fcrMemberMask{f} = mask;
        fcrMembers{f} = find(mask);
        ccap = inf(1,Kfcr);
        for r = 1:Kfcr
            cv = Rmat(fcrMembers{f}, r); cv = cv(cv ~= -1);
            if ~isempty(cv); ccap(r) = min(cv); end
        end
        fcrClassCap{f} = ccap;
        gv = Rmat(fcrMembers{f}, Kfcr+1); gv = gv(gv ~= -1);
        if ~isempty(gv); fcrGlobalCap(f) = min(gv); end
        if isfield(sn,'regionmaxmem') && numel(sn.regionmaxmem) >= f && ~isempty(sn.regionmaxmem{f})
            mv = sn.regionmaxmem{f}(fcrMembers{f}); mv = mv(mv ~= -1);
            if ~isempty(mv); fcrMemCap(f) = min(mv); end
        end
        fcrSz{f} = sn.regionsz(f,:);
        if isfield(sn,'regionlincon') && size(sn.regionlincon,1) >= f && ~isempty(sn.regionlincon{f,1})
            fcrA{f} = sn.regionlincon{f,1};
            fcrb{f} = sn.regionlincon{f,2};
        end
    end
    % see _kb/06-solver-catalog.md for rationale (SSA FCR)
    fcrRule = false(sn.nregions, Kfcr);
    if isfield(sn,'regionrule') && ~isempty(sn.regionrule)
        for f = 1:sn.nregions
            for r = 1:Kfcr
                fcrRule(f,r) = sn.regionrule(f,r) ~= DropStrategy.DROP;
            end
        end
    end
    fcrBuf = cell(sn.nregions,1);
    for f = 1:sn.nregions
        fcrBuf{f} = zeros(1,0);
    end
end

%%
if any(isinf(Np))
    Np(isinf(Np)) = 0;
end

init_state_hashed = ones(1,nstateful); % pick the first state in init_state{i}

%%
arvRatesSamples = zeros(options.samples,nstateful,R);
depRatesSamples = zeros(options.samples,nstateful,R);
% see _kb/09-ldes-and-cache.md (SSA: delayed-hit rate from the merge transition).
% Cache row layout is [srv(R) | var(V)], so the srv block ends V columns from
% the right; keep V per stateful index rather than assuming V == 0.
dlyRatesSamples = zeros(options.samples,nstateful,R);
% Derived START/PREEMPT rates, sampled exactly like the three above: the rate
% at which the transitions enabled in the current state start a class-r
% service, or push a class-r job in service back into the buffer.
startRatesSamples = zeros(options.samples,nstateful,R);
preemptRatesSamples = zeros(options.samples,nstateful,R);
cacheVarW = -ones(nstateful,1);
for ind=1:sn.nnodes
    if sn.nodetype(ind) == NodeType.Cache && sn.isstateful(ind)
        cacheVarW(sn.nodeToStateful(ind)) = sum(sn.nvars(ind,:));
    end
end
A = length(sync);
G = length(gsync);
samples_collected = 1;
nir = {};
% fill stateCell with initial states
cur_state = cell(nstateful,1); % cell array with current stateful node states
for ind=1:sn.nnodes
    if sn.isstateful(ind)
        isf = sn.nodeToStateful(ind);
        cur_state{isf} = init_state{isf}(init_state_hashed(isf),:);
        if sn.isstation(ind)
            ist = sn.nodeToStation(ind);
            [~,nir{ist}] = State.toMarginal(sn, ind, init_state{isf}(init_state_hashed(isf),:));
            nir{ist} = nir{ist}(:);
        end
    end
end
cur_state_1 = cur_state;
% generate state vector
state = cell2mat(cur_state');
% create function to determine lengths of stateful node states
statelen = cellfun(@length, cur_state);
% data structures to save transient information - pre-allocate for all samples
nSamples = options.samples;
tranSync = zeros(nSamples,1);
% tags of the transition that fires at each step, as rows [statefulIndex class]
tranStartTag = cell(nSamples,1);
tranPreemptTag = cell(nSamples,1);
tranState = zeros(1+length(state), nSamples);
tranState(1:(1+length(state)),1) = [0, state]';
SSq = zeros(length(cell2mat(nir')), nSamples);
SSq(:,1) = cell2mat(nir');
local = sn.nnodes+1;
last_node_a = 0; % active in the last occurred synchronization
last_node_p = 0; % passive in the last occurred synchronization
for act=1:A
    node_a{act} = sync{act}.active{1}.node;
    node_p{act} = sync{act}.passive{1}.node;
    class_a{act} = sync{act}.active{1}.class;
    class_p{act} = sync{act}.passive{1}.class;
    event_a{act} = sync{act}.active{1}.event;
    event_p{act} = sync{act}.passive{1}.event;
    outprob_a{act} = [];
    outprob_p{act} = [];
    start_a{act} = [];
    preempt_a{act} = [];
    start_p{act} = [];
    preempt_p{act} = [];
    % see _kb/06-solver-catalog.md for rationale (SSA immfeed self-loop)
    immfeed_selfloop{act} = false;
    if event_a{act}==EventType.DEP && node_p{act}==node_a{act} ...
            && node_a{act}>=1 && node_a{act}<=sn.nnodes && sn.isstation(node_a{act}) ...
            && isfield(sn,'immfeed') && ~isempty(sn.immfeed)
        istA_if = sn.nodeToStation(node_a{act});
        if istA_if>=1 && istA_if<=size(sn.immfeed,1) ...
                && class_p{act}>=1 && class_p{act}<=size(sn.immfeed,2) ...
                && sn.immfeed(istA_if, class_p{act})
            immfeed_selfloop{act} = true;
        end
    end
end
enabled_next_states = cell(1,A);

%% Start main simulation loop
isSimulation = true; % allow state vector to grow, e.g. for FCFS buffers
% see _kb/06-solver-catalog.md for rationale (SSA preamble ordering)
aectx = State.afterEventInit(sn);
samples_collected = 1;
cur_time = 0;
use_inline = true; % true = stable version, false = dev version

% Global (Whittle) rate scaling declared through setGlobalDependence. It reads
% the FULL population matrix, so it is a constant within one state and factors
% out of the per-transition rates, exactly as in SOLVER_CTMC.
hasGD = isfield(sn,'gdscaling') && ~isempty(sn.gdscaling);
gdNow = [];

try
    while samples_collected < options.samples && cur_time <= options.timespan(2) && ~lineTimeoutExceeded(options)
        %% This section corresponds to solver_ssa_findenabled in Java
        %% Inlined for performance reasons
        if use_inline
            enabled_sync = []; % row is action label, col1=rate, col2=new state
            enabled_rates = [];
            enabled_fcr = zeros(0,4); % [region class dest isSwitch] FCR marker per transition
            % derived tags of each enabled transition, as rows [statefulIndex class]
            enabled_tagS = {};
            enabled_tagP = {};
            ctr = 1;
            A = length(sync);
            G = length(gsync);
            if hasGD
                gdNow = solver_ssa_gdfactor(sn, cur_state);
            end
            % FCR: current aggregate per-class population of each region, used by
            % the arrival gate below to block entries that would exceed a cap.
            if fcrOn
                xcurFCR = cell(sn.nregions,1);
                for f = 1:sn.nregions
                    xf = zeros(1,sn.nclasses);
                    for i = fcrMembers{f}
                        ind_i = sn.stationToNode(i);
                        isf_i = sn.stationToStateful(i);
                        [~, nir_i] = State.toMarginalAggr(sn, ind_i, cur_state{isf_i});
                        xf = xf + nir_i(:)';
                    end
                    xcurFCR{f} = xf;
                end
            end
            for act=1:A
                isf_a = sn.nodeToStateful(node_a{act});

                enabled_next_states{act} = cur_state;
                update_cond_a = true;
                if update_cond_a
                    [enabled_next_states{act}{isf_a}, rate_a{act}, outprob_a{act}, eventCache, start_a{act}, preempt_a{act}] =  State.afterEvent(sn, node_a{act}, cur_state{isf_a}, event_a{act}, class_a{act}, isSimulation, eventCache, aectx, immfeed_selfloop{act});
                end

                if isempty(enabled_next_states{act}{isf_a}) || isempty(rate_a{act})
                    continue
                end

                % PHASE matters as much as DEP: refreshSync emits phase moves as
                % active station events, so phase-type service would otherwise
                % advance unscaled.
                if hasGD && sn.isstation(node_a{act}) && (event_a{act} == EventType.DEP || event_a{act} == EventType.PHASE)
                    rate_a{act} = rate_a{act} * gdNow(sn.nodeToStation(node_a{act}), class_a{act});
                end

                for ia=1:size(enabled_next_states{act}{isf_a},1) % for all possible new states, check if they are enabled
                    % if the transition cannot occur
                    if isnan(rate_a{act}(ia)) || rate_a{act}(ia) == 0 % handles degenerate rate values
                        % set the transition with a zero rate so that it is
                        % never selected
                        rate_a{act}(ia) = 1e-38; % ~ zero in 32-bit precision
                    end

                    if enabled_next_states{act}{isf_a}(ia,:) == -1 % hash not found
                        continue
                    end
                    % A delayed hit is the ONLY cache transition that empties the
                    % node: the request merges onto the in-flight fetch and is held
                    % in block B, so it departs later in the hit class and is
                    % otherwise indistinguishable there from a true hit.
                    isMergeA = false;
                    if cacheVarW(isf_a) >= 0 && event_a{act} == EventType.READ
                        rowA = enabled_next_states{act}{isf_a}(ia,:);
                        preA = cur_state{isf_a};
                        eA = numel(rowA) - cacheVarW(isf_a);
                        eP = numel(preA) - cacheVarW(isf_a);
                        if eA >= R && eP >= R
                            isMergeA = (sum(rowA((eA-R+1):eA)) - sum(preA((eP-R+1):eP))) == -1;
                        end
                    end
                    update_cond_p = true; %samples_collected == 1 || ((node_p{act} == last_node_a || node_p{act} == last_node_p)) || isempty(outprob_a{act}) || isempty(outprob_p{act});

                    if rate_a{act}(ia)>0
                        if node_p{act} ~= local
                            if node_p{act} == node_a{act} %self-loop, active and passive are the same
                                isf_p = isf_a;
                                if update_cond_p
                                    [enabled_next_states{act}{isf_p}, ~, outprob_p{act}, eventCache, start_p{act}, preempt_p{act}] =  State.afterEvent(sn, node_p{act}, enabled_next_states{act}{isf_p}, event_p{act}, class_p{act}, isSimulation, eventCache, aectx);
                                end
                            else % departure
                                isf_p = sn.nodeToStateful(node_p{act});
                                if update_cond_p
                                    [enabled_next_states{act}{isf_p}, ~, outprob_p{act}, eventCache, start_p{act}, preempt_p{act}] =  State.afterEvent(sn, node_p{act}, enabled_next_states{act}{isf_p}, event_p{act}, class_p{act}, isSimulation, eventCache, aectx);
                                end
                            end
                            if ~isempty(enabled_next_states{act}{isf_p})
                                if sn.isstatedep(node_a{act},3)
                                    prob_sync_p{act} = sync{act}.passive{1}.prob(cur_state, enabled_next_states{act}); %state-dependent
                                else
                                    prob_sync_p{act} = sync{act}.passive{1}.prob;
                                end
                            else
                                prob_sync_p{act} = 0;
                            end
                        end
                        if ~isempty(enabled_next_states{act}{isf_a})
                            if node_p{act} == local
                                prob_sync_p{act} = 1;
                            end
                            if ~isnan(rate_a{act})
                                if all(~cellfun(@isempty,enabled_next_states{act}))
                                    % see _kb/06-solver-catalog.md for rationale (SSA FCR)
                                    blockFCR = false;
                                    fcrMark = [0 0 0 0]; % [region class dest isSwitch]
                                    if fcrOn && node_p{act} ~= local && node_p{act} <= sn.nnodes
                                        jp = sn.nodeToStation(node_p{act});
                                        if jp > 0
                                            ja = sn.nodeToStation(node_a{act});
                                            cc = class_p{act};
                                            for f = 1:sn.nregions
                                                mmask = fcrMemberMask{f};
                                                if mmask(jp) && (ja <= 0 || ja > numel(mmask) || ~mmask(ja))
                                                    xn = xcurFCR{f}; xn(cc) = xn(cc) + 1;
                                                    if xn(cc) > fcrClassCap{f}(cc) || sum(xn) > fcrGlobalCap(f) ...
                                                            || (xn * fcrSz{f}(:) > fcrMemCap(f)) ...
                                                            || (~isempty(fcrA{f}) && any(fcrA{f} * xn(:) > fcrb{f}(:)))
                                                        if fcrRule(f,cc)
                                                            fcrMark = [f cc node_p{act} 0]; % park in FIFO
                                                        else
                                                            fcrMark = [f cc node_p{act} 2]; % DROP: destroyed
                                                        end
                                                        break;
                                                    end
                                                elseif mmask(jp) && ja > 0 && ja <= numel(mmask) && mmask(ja) ...
                                                        && cc ~= class_a{act}
                                                    if fcrRule(f,cc)
                                                        fcrMark = [f cc node_p{act} 1]; % exit + gated re-entry
                                                    else
                                                        fcrMark = [f cc node_p{act} 3]; % exit + gated re-entry, DROP on refusal
                                                    end
                                                    break;
                                                end
                                            end
                                        end
                                    end
                                    if fcrMark(1) > 0
                                        % see _kb/06-solver-catalog.md for rationale (SSA FCR)
                                        enabled_next_states{act}{isf_p} = cur_state{isf_p};
                                    end
                                    if event_a{act} == EventType.DEP && ~blockFCR
                                        node_a_sf{act} = isf_a;
                                        node_p_sf{act} = isf_p;
                                        depRatesSamples(samples_collected,node_a_sf{act},class_a{act}) = depRatesSamples(samples_collected,node_a_sf{act},class_a{act}) + outprob_a{act} * outprob_p{act} * rate_a{act}(ia) * prob_sync_p{act};
                                        arvRatesSamples(samples_collected,node_p_sf{act},class_p{act}) = arvRatesSamples(samples_collected,node_p_sf{act},class_p{act}) + outprob_a{act} * outprob_p{act} * rate_a{act}(ia) * prob_sync_p{act};
                                    end
                                    % simulate also self-loops as we need to log them
                                    %if any(~cellfun(@isequal,new_state{act},cur_state))
                                    if node_p{act} < local && ~sn.csmask(class_a{act}, class_p{act}) && sn.nodetype(node_p{act})~=NodeType.Source && (rate_a{act}(ia) * prob_sync_p{act} >0)
                                        line_error(mfilename,sprintf('Error: state-dependent routing at node %d (%s) violates the class switching mask (node %d -> node %d, class %d -> class %d).', node_a{act}, sn.nodenames{node_a{act}}, node_a{act}, node_p{act}, class_a{act}, class_p{act}));
                                    end
                                    if isMergeA && ~blockFCR
                                        % Weighted like enabled_rates, NOT like
                                        % depRates: for a simulated READ the item
                                        % is already sampled from pread, so folding
                                        % outprob_a back in would count p(k) twice.
                                        dlyRatesSamples(samples_collected,isf_a,class_a{act}) = ...
                                            dlyRatesSamples(samples_collected,isf_a,class_a{act}) + rate_a{act}(ia) * prob_sync_p{act};
                                    end
                                    % START/PREEMPT tags of this arc, weighted like
                                    % enabled_rates: they annotate the transition
                                    % itself, so their rate is its rate. Written for
                                    % EVERY action, not only for departures -- a
                                    % retrial or a polling switchover starts service
                                    % without being a DEP, and the arrival half of a
                                    % departure is where most starts happen.
                                    tagS_here = zeros(0,2);
                                    tagP_here = zeros(0,2);
                                    if ~blockFCR
                                        w_tag = rate_a{act}(ia) * prob_sync_p{act};
                                        if ~isempty(start_a{act}) && ia <= size(start_a{act},1)
                                            startRatesSamples(samples_collected,isf_a,:) = reshape(startRatesSamples(samples_collected,isf_a,:),1,R) + w_tag * start_a{act}(ia,:);
                                            preemptRatesSamples(samples_collected,isf_a,:) = reshape(preemptRatesSamples(samples_collected,isf_a,:),1,R) + w_tag * preempt_a{act}(ia,:);
                                            tagS_here = [tagS_here; sub_tagrows(isf_a, start_a{act}(ia,:))]; %#ok<AGROW>
                                            tagP_here = [tagP_here; sub_tagrows(isf_a, preempt_a{act}(ia,:))]; %#ok<AGROW>
                                        end
                                        if node_p{act} ~= local && ~isempty(start_p{act})
                                            startRatesSamples(samples_collected,isf_p,:) = reshape(startRatesSamples(samples_collected,isf_p,:),1,R) + w_tag * start_p{act}(1,:);
                                            preemptRatesSamples(samples_collected,isf_p,:) = reshape(preemptRatesSamples(samples_collected,isf_p,:),1,R) + w_tag * preempt_p{act}(1,:);
                                            tagS_here = [tagS_here; sub_tagrows(isf_p, start_p{act}(1,:))]; %#ok<AGROW>
                                            tagP_here = [tagP_here; sub_tagrows(isf_p, preempt_p{act}(1,:))]; %#ok<AGROW>
                                        end
                                    end
                                    if ~blockFCR
                                        enabled_rates(ctr) = rate_a{act}(ia) * prob_sync_p{act};
                                        enabled_sync(ctr) = act;
                                        enabled_fcr(ctr,:) = fcrMark;
                                        enabled_tagS{ctr} = tagS_here;
                                        enabled_tagP{ctr} = tagP_here;
                                        ctr = ctr + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            gctr_start = ctr;

            for gact=1:G % event at node ind with global side-effects
                gind = gsync{gact}.active{1}.node; % node index for global event
                [enabled_next_states{A+gact}, outrate, outprob] = State.afterGlobalEvent(sn, gind, cur_state, gsync{gact}, isSimulation);
                for ia=find(outrate .* outprob)
                    enabled_rates(ctr) = outrate(ia) * outprob(ia);
                    enabled_sync(ctr) = A+gact;
                    % an SPN transition holds no server: no service starts there
                    enabled_tagS{ctr} = zeros(0,2);
                    enabled_tagP{ctr} = zeros(0,2);
                    ctr = ctr + 1;

                    % Record departure/arrival rates for FIRE events at Places
                    if gsync{gact}.active{1}.event == EventType.FIRE
                        mode = gsync{gact}.active{1}.mode;
                        % Get enabling/firing conditions to determine affected classes
                        enabling_m = sn.nodeparam{gind}.enabling{mode};
                        firing_m = sn.nodeparam{gind}.firing{mode};

                        for j=1:length(gsync{gact}.passive)
                            pev = gsync{gact}.passive{j};
                            % Decode linear index to (node, class) - pev.node is a linear index from find() on enabling/firing matrix
                            [pev_node, pev_class] = ind2sub([sn.nnodes, R], pev.node);
                            if pev.event == EventType.PRE
                                % Departure from input Place (consuming tokens)
                                if pev_node <= length(sn.nodeToStateful) && ~isnan(sn.nodeToStateful(pev_node)) && sn.nodeToStateful(pev_node) > 0
                                    ep_isf = sn.nodeToStateful(pev_node);
                                    % Record departures for the specific class from this PRE event
                                    depRatesSamples(samples_collected, ep_isf, pev_class) = ...
                                        depRatesSamples(samples_collected, ep_isf, pev_class) + outrate(ia) * outprob(ia);
                                end
                            elseif pev.event == EventType.POST
                                % Arrival at output Place (producing tokens)
                                if pev_node <= length(sn.nodeToStateful) && ~isnan(sn.nodeToStateful(pev_node)) && sn.nodeToStateful(pev_node) > 0
                                    fp_isf = sn.nodeToStateful(pev_node);
                                    % Record arrivals for the specific class from this POST event
                                    arvRatesSamples(samples_collected, fp_isf, pev_class) = ...
                                        arvRatesSamples(samples_collected, fp_isf, pev_class) + outrate(ia) * outprob(ia);
                                end
                            end
                        end
                    end
                end
            end

            % see _kb/06-solver-catalog.md for rationale (native fork-join)
            FJ = 0;
            if isfield(sn,'fjsync') && ~isempty(sn.fjsync)
                FJ = length(sn.fjsync);
            end
            for fjact=1:FJ
                [fjStates, fjrate, fjprob] = State.afterFJEvent(sn, sn.fjsync{fjact}, cur_state, isSimulation);
                if ~isempty(fjStates)
                    enabled_next_states{A+G+fjact} = fjStates{1};
                    enabled_rates(ctr) = fjrate(1) * fjprob(1);
                    enabled_sync(ctr) = A+G+fjact;
                    ctr = ctr + 1;
                    fjentry = sn.fjsync{fjact};
                    isf_fork = sn.nodeToStateful(fjentry.fork);
                    depRatesSamples(samples_collected, isf_fork, fjentry.class) = ...
                        depRatesSamples(samples_collected, isf_fork, fjentry.class) + fjrate(1) * fjprob(1);
                    for b=1:length(fjentry.branchheads)
                        isf_bh = sn.nodeToStateful(fjentry.branchheads(b));
                        arvRatesSamples(samples_collected, isf_bh, fjentry.auxclasses(b)) = ...
                            arvRatesSamples(samples_collected, isf_bh, fjentry.auxclasses(b)) + fjrate(1) * fjprob(1);
                    end
                end
            end
        else
            [enabled_next_states,enabled_rates,enabled_sync,gctr_start,depRatesSamples,arvRatesSamples,outprob_a,outprob_p,rate_a, eventCache,startRatesSamples,preemptRatesSamples,enabled_tagS,enabled_tagP] = solver_ssa_findenabled(sn,node_a,enabled_next_states,cur_state,outprob_a,event_a,class_a,isSimulation,node_p,local,outprob_p,event_p,class_p,sync,gsync,depRatesSamples,samples_collected,arvRatesSamples,last_node_a,last_node_p,eventCache,startRatesSamples,preemptRatesSamples);
        end
        %% Gillespie direct method
        tot_rate = sum(enabled_rates);
        cum_rate = cumsum(enabled_rates) / tot_rate;
        selected_transition = 1 + max([0,find( rand > cum_rate )]); % select action

        % Update record of last active/passive pair
        if isempty(enabled_sync)
            line_error(mfilename,'SSA simulation entered a deadlock before collecting all samples, no synchronization is enabled.');
        end
        if selected_transition < gctr_start
            % regular event pair
            last_node_a = node_a{enabled_sync(selected_transition)};
            last_node_p = node_p{enabled_sync(selected_transition)};
        else % global event
            last_node_a = NaN;
            last_node_p = NaN;
        end

        %% Update paddings
        % see _kb/06-solver-catalog.md for rationale (SSA state left-padding)
        for ind=1:sn.nnodes
            if sn.isstation(ind)
                isf = sn.nodeToStateful(ind);
                deltalen = length(cur_state{isf}) - statelen(isf);
                if deltalen>0
                    statelen(isf) = length(cur_state{isf});
                    % Rows before this node's block; row 1 of tranState is dt, so
                    % the block starts at shift+2. Keyed on the STATEFUL index,
                    % never on the node index: a station whose node is first need
                    % not be the first stateful node, and sum(x(1:0)) is already 0.
                    shift = sum(statelen(1:isf-1));
                    % INSERT deltalen rows at the LEFT of the block -- the buffer
                    % is right-aligned, so a widened row's fresh slots are its
                    % leftmost ones and every earlier sample must gain exactly
                    % those. Resuming the tail at (shift+1+deltalen) instead DROPPED
                    % deltalen-1 rows of real history and grew the matrix by one row
                    % rather than by deltalen. At deltalen == 1 the two agree, which
                    % is why every per-class-count buffer was unaffected; the
                    % preemptive families store (class,phase) PAIRS and grow by 2,
                    % so their history was silently re-encoded and unique() then
                    % conflated states that differ only in which class holds the
                    % server. On Source -> FCFSPRPRIO -> Sink at lambda = 0.08 per
                    % class that reported TN = [0.018 0.142] for an exact
                    % [0.08 0.08], the total right and the split wrong.
                    pad = zeros(deltalen, size(tranState,2));
                    tranState = [tranState(1:(shift+1), :); pad ; tranState((shift+2):end, :)];
                end
            end
        end

        %% Simulate the time increment
        state = cell2mat(cur_state');
        dt = -(log(rand)/tot_rate);
        cur_time = cur_time + dt;

        %% Save simulation output data
        tranState(1:(1+length(state)),samples_collected) = [dt, state]';
        tranSync(samples_collected,1) = enabled_sync(selected_transition);
        % the derived tags of the transition that fired; sn.sync carries none,
        % so the trace would otherwise have no way to report them
        if selected_transition <= numel(enabled_tagS)
            tranStartTag{samples_collected} = enabled_tagS{selected_transition};
            tranPreemptTag{samples_collected} = enabled_tagP{selected_transition};
        end
        for ind=1:sn.nnodes
            if sn.isstation(ind)
                isf = sn.nodeToStateful(ind);
                ist = sn.nodeToStation(ind);
                [~,nir{ist}] = State.toMarginal(sn, ind, cur_state{isf});
                nir{ist}=nir{ist}(:);
            end
        end
        SSq(:,samples_collected) = cell2mat(nir');

        %% Update current state and sample counter
        cur_state_1 = cur_state;
        cur_state = enabled_next_states{enabled_sync(selected_transition)};

        %% FCR WAITQ bookkeeping (see _kb/06-solver-catalog.md, SSA FCR)
        if fcrOn && use_inline
            pk = [0 0 0 0];
            if selected_transition <= size(enabled_fcr,1)
                pk = enabled_fcr(selected_transition,:);
            end
            if pk(1) > 0 && pk(4) == 0
                % blocked entry: park the (class, destination) token
                fcrBuf{pk(1)}(end+1) = (pk(3)-1)*R + pk(2);
            end
            % pk(4)==2: DROP, the refused job was destroyed (active part only)
            % release cascade over all regions
            [cur_state, fcrBuf, eventCache] = fcr_release(sn, cur_state, fcrBuf, ...
                fcrMembers, fcrClassCap, fcrGlobalCap, fcrMemCap, fcrSz, fcrA, fcrb, ...
                isSimulation, eventCache, aectx);
            if pk(1) > 0 && (pk(4) == 1 || pk(4) == 3)
                % class-switching hop: gate the re-entry after the cascade
                f_ = pk(1); cls_ = pk(2); dest_ = pk(3);
                xf_ = fcr_regionpop(sn, cur_state, fcrMembers{f_});
                xn_ = xf_; xn_(cls_) = xn_(cls_) + 1;
                admitted = false;
                if ~fcr_violates(xn_, fcrClassCap{f_}, fcrGlobalCap(f_), fcrMemCap(f_), fcrSz{f_}, fcrA{f_}, fcrb{f_})
                    isf_d = sn.nodeToStateful(dest_);
                    [ns_, ~, ~, eventCache] = State.afterEvent(sn, dest_, cur_state{isf_d}, EventType.ARV, cls_, isSimulation, eventCache, aectx);
                    if ~isempty(ns_)
                        cur_state{isf_d} = ns_(1,:);
                        admitted = true;
                    end
                end
                if ~admitted && pk(4) == 1
                    fcrBuf{f_}(end+1) = (dest_-1)*R + cls_;
                end
                % pk(4)==3 refused: DROP, the switching job is destroyed
            end
        end

        samples_collected = samples_collected + 1;

        %% Print progress
        print_progress(options,samples_collected,cur_time);
    end
    % The counter row is closed here rather than newline-terminated:
    % line_printf already ends an open row, so an explicit newline was a
    % SECOND one and showed as a blank row before the completion banner.
    LineStatus.close();
catch ME
    getReport(ME)
end

% Trim pre-allocated arrays to actual number of samples collected
samples_collected = samples_collected - 1;  % Adjust for the increment at end of loop
tranState = tranState(:, 1:samples_collected);
tranSync = tranSync(1:samples_collected, :);
tranStartTag = tranStartTag(1:samples_collected);
tranPreemptTag = tranPreemptTag(1:samples_collected);
SSq = SSq(:, 1:samples_collected);

% see _kb/06-solver-catalog.md for rationale (SSA warmup discard)
warmupfrac = 0.0;
if isfield(options, 'config') && isfield(options.config, 'warmupfrac')
    warmupfrac = max(0.0, min(0.99, options.config.warmupfrac));
end
if warmupfrac > 0 && samples_collected > 1
    nDrop = floor(warmupfrac * samples_collected);
    if nDrop > 0 && nDrop < samples_collected
        keep = (nDrop+1):samples_collected;
        tranState = tranState(:, keep);
        tranSync = tranSync(keep, :);
        tranStartTag = tranStartTag(keep);
        tranPreemptTag = tranPreemptTag(keep);
        SSq = SSq(:, keep);
        if exist('arvRatesSamples', 'var') && ~isempty(arvRatesSamples) ...
                && size(arvRatesSamples, 1) >= samples_collected
            arvRatesSamples = arvRatesSamples(keep, :, :);
        end
        if exist('depRatesSamples', 'var') && ~isempty(depRatesSamples) ...
                && size(depRatesSamples, 1) >= samples_collected
            depRatesSamples = depRatesSamples(keep, :, :);
        end
        if exist('dlyRatesSamples', 'var') && ~isempty(dlyRatesSamples) ...
                && size(dlyRatesSamples, 1) >= samples_collected
            dlyRatesSamples = dlyRatesSamples(keep, :, :);
        end
        if size(startRatesSamples, 1) >= samples_collected
            startRatesSamples = startRatesSamples(keep, :, :);
            preemptRatesSamples = preemptRatesSamples(keep, :, :);
        end
        samples_collected = numel(keep);
    end
end

tranState = tranState';


[u,ui,uj] = unique(tranState(:,2:end),'rows');
statesz = cellfun(@length, cur_state_1)';
tranSysState = cell(1,length(cur_state)+1);
tranSysState{1} = cumsum(tranState(:,1));
for j=1:length(statesz)
    tranSysState{1+j} = tranState(:,1+(1+sum(statesz(1:(j-1)))):(1+sum(statesz(1:j))));
end
arvRates = zeros(size(u,1),sn.nstateful,R);
depRates = zeros(size(u,1),sn.nstateful,R);
dlyRates = zeros(size(u,1),sn.nstateful,R);
startRates = zeros(size(u,1),sn.nstateful,R);
preemptRates = zeros(size(u,1),sn.nstateful,R);

pi = zeros(1,size(u,1));
for s=1:size(u,1)
    pi(s) = sum(tranState(uj==s,1));
end
SSq = SSq(:,ui)'; % we restrict to unique states in the simulation

for ind=1:sn.nnodes
    if sn.isstateful(ind)
        isf = sn.nodeToStateful(ind);
        if sn.isstation(ind)
            ist = sn.nodeToStation(ind);
            %K = sn.phasessz(ist,:);
            %Ks = sn.phaseshift(ist,:);
        end
        for s=1:size(u,1)
            for r=1:R
                arvRates(s,isf,r) = arvRatesSamples(ui(s),isf,r); % for each unique state, one (any) sample of the rate is enough here
                depRates(s,isf,r) = depRatesSamples(ui(s),isf,r); % for each unique state, one (any) sample of the rate is enough here
                % the tag rates are a deterministic function of the state too,
                % so one visit gives them exactly, as for the two above
                startRates(s,isf,r) = startRatesSamples(ui(s),isf,r);
                preemptRates(s,isf,r) = preemptRatesSamples(ui(s),isf,r);
            end
        end
        % see _kb/09-ldes-and-cache.md (SSA: the merge rate needs EVERY visit).
        % Unlike a DEP rate, the merge rate is random given the state, so one
        % visit is a single Bernoulli draw whose variance does not shrink with
        % the sample count.
        if cacheVarW(isf) >= 0
            nsmp = numel(uj); % dlyRatesSamples is preallocated, tranState is not
            visitCnt = accumarray(uj, 1, [size(u,1) 1]);
            for r=1:R
                dlySum = accumarray(uj, dlyRatesSamples(1:nsmp,isf,r), [size(u,1) 1]);
                dlyRates(:,isf,r) = dlySum ./ max(visitCnt,1);
            end
        end
    end
end
pi = pi/sum(pi);
%sn.nservers = init_nserver; % restore Inf at delay nodes
end

function print_progress(options,samples_collected,cur_time)
if LineConsole.isActive()
    % the console owns the line: report a decimated progress row instead of
    % the in-place counter, which a paged log cannot rewrite
    every = max(1,round(options.samples/20));
    if mod(samples_collected, every) == 0
        LineConsole.iter(samples_collected/every, ...
            'simulated %d of %g samples (%.0f%%), simulated time %.4g', ...
            samples_collected, options.samples, ...
            100*samples_collected/options.samples, cur_time);
    end
    return
end
if options.verbose && ~batchStartupOptionUsed
    % ONE REWRITTEN FIELD, not a fixed-width one. LineStatus rewinds by
    % the width it actually wrote, so a counter that only grows needs no
    % padding at all and leaves no trailing blanks -- a fixed %-9d field
    % showed its pad as "SSA samples: 100000   ". It also cannot desync
    % the way a hardcoded run of backspaces does once the count outgrows
    % the field. line_printf closes the row, so the completion banner
    % terminates it without help.
    if options.verbose == 2 || (samples_collected > 0 && mod(samples_collected,1e2) == 0)
        LineStatus.set('SSA samples: %d', samples_collected);
    end
end
end
function rows = sub_tagrows(isf, tagrow)
% ROWS=SUB_TAGROWS(ISF,TAGROW) expand a per-class tag count into one
% [statefulIndex class] row per tagged job, so the trace can report each of
% them separately. A count is a whole number on every path but the merged
% destinations of a G-network signal, where it is an expectation; the trace
% enumerates events, so it takes the integer part and the exact fractional
% value stays in the rate counters.
rows = zeros(0,2);
for r = find(tagrow(:)' > 0)
    for c = 1:tagrow(r)
        rows(end+1,:) = [isf, r]; %#ok<AGROW>
    end
end
end

function x = fcr_regionpop(sn, cur_state, members)
% X=FCR_REGIONPOP(SN,CUR_STATE,MEMBERS) per-class population of a finite
% capacity region given the current state cells
x = zeros(1, sn.nclasses);
for i = members
    ind_i = sn.stationToNode(i);
    isf_i = sn.stationToStateful(i);
    [~, nir_i] = State.toMarginalAggr(sn, ind_i, cur_state{isf_i});
    x = x + nir_i(:)';
end
end

function tf = fcr_violates(xn, ccap, gcap, memcap, sz, A, b)
% TF=FCR_VIOLATES(...) true if population vector xn breaks any admission
% constraint of the region
tf = any(xn > ccap) || sum(xn) > gcap || (xn * sz(:) > memcap);
if ~tf && ~isempty(A)
    tf = any(A * xn(:) > b(:));
end
end

function [cur_state, fcrBuf, eventCache] = fcr_release(sn, cur_state, fcrBuf, ...
    fcrMembers, fcrClassCap, fcrGlobalCap, fcrMemCap, fcrSz, fcrA, fcrb, ...
    isSimulation, eventCache, aectx)
% FCR_RELEASE strict-FIFO head-of-line release of parked region tokens:
% admit heads while the admission constraints permit, applying the arrival
% to the destination station state
K = sn.nclasses;
progress = true;
while progress
    progress = false;
    for f = 1:length(fcrBuf)
        if isempty(fcrBuf{f})
            continue
        end
        x = fcr_regionpop(sn, cur_state, fcrMembers{f});
        tok = fcrBuf{f}(1);
        dest = floor((tok-1)/K) + 1;
        r = mod(tok-1, K) + 1;
        xn = x; xn(r) = xn(r) + 1;
        if fcr_violates(xn, fcrClassCap{f}, fcrGlobalCap(f), fcrMemCap(f), fcrSz{f}, fcrA{f}, fcrb{f})
            continue % head-of-line: this region's FIFO stays blocked
        end
        isf_d = sn.nodeToStateful(dest);
        [ns, ~, ~, eventCache] = State.afterEvent(sn, dest, cur_state{isf_d}, EventType.ARV, r, isSimulation, eventCache, aectx);
        if isempty(ns)
            continue % destination cannot accept (e.g. station capacity)
        end
        cur_state{isf_d} = ns(1,:);
        fcrBuf{f}(1) = [];
        progress = true;
    end
end
end
