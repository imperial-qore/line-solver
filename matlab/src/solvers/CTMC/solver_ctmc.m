function [Q,stateSpace,stateSpaceAggr,Dfilt,arvRates,depRates,sn,DfiltAux]=solver_ctmc(sn,options)
% [Q,SS,SSQ,DFILT,ARVRATES,DEPRATES,QN,DFILTAUX]=SOLVER_CTMC(QN,OPTIONS)
%
% DFILTAUX carries the two DERIVED event filtrations, DFILTAUX.start{i,r} and
% DFILTAUX.preempt{i,r}: the (state x state) rate at which a class-r job
% begins holding a server at station i, and the rate at which one is pushed
% back into the buffer there. They are deliberately NOT appended to DFILT: a
% START rides on the SAME arc as the ARV or DEP that causes it, so adding it
% to the eventFilt cell would double-count D1 in @SolverCTMC/sample.m
% (D0 = infGen - sum(eventFilt)) and break the pairing with sn.sync that
% getGenerator.m asserts. See _kb/06-solver-catalog.md.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% impatience and server-count support checks
% The rules live in SOLVER_CTMC_STATE_SUPPORTS, which SolverCTMC.supportsModelMethod
% and SolverSSA ask as well: one body, so the report and this run cannot differ.
[stateOk, stateWhy] = solver_ctmc_state_supports(sn, 'SolverCTMC');
if ~stateOk
    line_error(mfilename, stateWhy);
end
% see _kb/06-solver-catalog.md (CTMC section, support gates) for rationale, incl. the REPLY-signal exception
if isfield(sn,'issignal') && ~isempty(sn.issignal) && any(sn.issignal)
    annihilated = sn.issignal(:)';
    if isfield(sn,'signaltype') && ~isempty(sn.signaltype)
        for rr = 1:sn.nclasses
            if annihilated(rr) && numel(sn.signaltype) >= rr && ~isempty(sn.signaltype{rr}) ...
                    && ~any(isnan(sn.signaltype{rr})) && sn.signaltype{rr} == SignalType.REPLY
                annihilated(rr) = false;
            end
        end
    end
    for ii = 1:sn.nstations
        if sn.sched(ii) ~= SchedStrategy.EXT
            sn.classcap(ii, annihilated) = 0;
        end
    end
end
% see _kb/06-solver-catalog.md (CTMC section, heterogeneous servers) for rationale
for ind = 1:sn.nnodes
    if sn.isstation(ind) && isfield(sn,'nodeparam') && numel(sn.nodeparam) >= ind ...
            && ~isempty(sn.nodeparam{ind}) && isstruct(sn.nodeparam{ind}) ...
            && isfield(sn.nodeparam{ind},'nservertypes') && sn.nodeparam{ind}.nservertypes > 0
        ist = sn.nodeToStation(ind);
        % PAS/OI stations model heterogeneous compatible servers through the OI
        % rank rate (svcRateFun), not a single-class load-dependent scaling.
        if sn.sched(ist) == SchedStrategy.PAS || sn.sched(ist) == SchedStrategy.OI
            continue
        end
        np = sn.nodeparam{ind};
        served = [];
        for r = 1:sn.nclasses
            if ~isempty(sn.proc{ist}{r}) && ~any(any(isnan(sn.proc{ist}{r}{1}))) && sn.rates(ist,r) > 0
                served(end+1) = r; %#ok<AGROW>
            end
        end
        if numel(served) > 1
            line_error(mfilename,'SolverCTMC supports heterogeneous servers only for single-class stations. Use SolverJMT or SolverLDES for multi-class heterogeneous servers.');
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

%% generate state space
%nnodes = sn.nnodes;
nstateful = sn.nstateful;
nclasses = sn.nclasses;
sync = sn.sync;
% True when at least one service or arrival process is a matrix exponential, so
% that the generator legitimately carries negative off-diagonal entries.
hasMEproc = isfield(sn,'isph') && ~isempty(sn.isph) && ~all(sn.isph(:));
A = length(sync);
csmask = sn.csmask;

line_debug('CTMC solver starting: nstateful=%d, nclasses=%d, sync_events=%d', nstateful, nclasses, A);

if ~isfield(options.config, 'hide_immediate')
    options.config.hide_immediate = true;
end

if ~isfield(options.config, 'state_space_gen')
    options.config.state_space_gen = 'default';
end

%% generate state spaces, detailed and aggregate
switch options.config.state_space_gen
    case 'reachable' % does not handle open models yet (no cutoff)
        line_debug('Using reachable state space generation, calling ctmc_ssg_reachability');
        [stateSpace, stateSpaceAggr, stateSpaceHashed,~,sn] = ctmc_ssg_reachability(sn,options);
    case {'default','full'}
        line_debug('Using full state space generation, calling ctmc_ssg');
        [stateSpace, stateSpaceAggr, stateSpaceHashed,~,sn] = ctmc_ssg(sn,options);
end

line_debug('State space generated: %d states', size(stateSpaceHashed,1));

%% Finite Capacity Region handling
% see _kb/06-solver-catalog.md (CTMC section) for rationale
fcrWaitq = isfield(sn,'nregions') && sn.nregions > 0;

%% Global (Whittle) dependence: one evaluation of phi(n) per state
% phi reads the FULL population matrix, so within a state it is a constant that
% factors out of every rate at that state; see _kb/11-conventions-and-gotchas.md
hasGD = isfield(sn,'gdscaling') && ~isempty(sn.gdscaling);
gdFactor = [];
if hasGD
    if fcrWaitq
        line_error(mfilename,'A global dependence (setGlobalDependence) cannot be combined with finite capacity regions: the region generator builds its own transitions and would ignore the scaling.');
    end
    gdFactor = solver_ctmc_gdfactor(sn, stateSpaceAggr, options);
end

%%
if fcrWaitq
    % see _kb/06-solver-catalog.md (CTMC section) for rationale
    [stateSpace,stateSpaceAggr,stateSpaceHashed,Dfilt,sn,basBlockQ,DfiltAux] = solver_ctmc_fcr_waitq(sn,options);
    Q = speye(size(stateSpaceHashed,1)); % the diagonal elements will be removed later
    % see _kb/06-solver-catalog.md (True BAS blocking) for rationale
    Qimm = 0*Q;
else
Q = speye(size(stateSpaceHashed,1)); % the diagonal elements will be removed later
Dfilt = cell(1,A);
for a=1:A
    Dfilt{a} = 0*Q;
end
% Derived START/PREEMPT filtrations, one sparse matrix per (station, class).
% They live outside Dfilt on purpose: their arcs are the same arcs.
DfiltAux = solver_ctmc_auxfilt_init(sn, Q);
% see _kb/06-solver-catalog.md (True BAS blocking) for rationale
basBlockQ = 0*Q;
% see _kb/06-solver-catalog.md (Vanishing states) for rationale
Qimm = 0*Q;
local = sn.nnodes+1; % passive action

% SPN code
% Adj_t = zeros(size(SSh,1),size(SSh,1));
% Adj_m = zeros(size(SSh,1),size(SSh,1));
% if ~isempty(Adj) && ~isempty(ST)
%    edges = adj_to_mat(Adj);
% end

%% for all synchronizations
for a=1:A
    stateCell = cell(nstateful,1);
    %sn.sync{a}.active{1}.print
    for s=1:size(stateSpaceHashed,1)
        %[a,s]
        state = stateSpaceHashed(s,:);
        % SPN code
        %         ustate = stateSpace(s,:);
        %         state_pn = [];
        %         for st=1:length(ustate)
        %             if ~isempty(sn.varsparam{st}) && isfield(sn.varsparam{st}, 'nodeToPlace')
        %                 state_pn(sn.varsparam{st}.nodeToPlace) = ustate(st);
        %             end
        %         end

        % update state cell array and SSq
        for ind = 1:sn.nnodes
            if sn.isstateful(ind)
                isf = sn.nodeToStateful(ind);
                stateCell{isf} = sn.space{isf}(state(isf),:);
                %                if sn.isstation(ind)
                %                    ist = sn.nodeToStation(ind);
                %                    [~,nir] = State.toMarginal(sn,ind,stateCell{isf});
                %                end
            end
        end
        node_a = sync{a}.active{1}.node;
        state_a = state(sn.nodeToStateful(node_a));
        class_a = sync{a}.active{1}.class;
        event_a = sync{a}.active{1}.event;
        [new_state_a, rate_a, ~, start_a, preempt_a] = State.afterEventHashed( sn, node_a, state_a, event_a, class_a);
        % SPN code:
        %[new_state_a, rate_a,~,trans_a, modes_a] = State.afterEventHashed( qn, node_a, state_a, event_a, class_a);
        if hasGD && sn.isstation(node_a) && (event_a == EventType.DEP || event_a == EventType.PHASE)
            % phi(n) is constant across the rows of rate_a at this state, so it
            % factors out; PHASE is scaled too or phase-type service would
            % advance unscaled (see refreshSync.m, which emits it as active)
            rate_a = rate_a * gdFactor(s, sn.nodeToStation(node_a) + (class_a-1)*sn.nstations);
        end

        %% debugging block
        %        if true%options.verbose == 2
        %            line_printf('---\n');
        %            sync{a}.active{1}.print,
        %        end
        %%
        if new_state_a == -1 % hash not found
            continue
        end
        for ia=1:length(new_state_a)
            % A matrix-exponential process embeds in the generator exactly as a
            % phase-type does, except that the off-diagonal entries of D0 and the
            % completion vector -A*e may be negative. Those transitions are part
            % of the balance equations: dropping them leaves the diagonal to
            % absorb their mass and silently answers a different model (an
            % M/CME/1 lost 3% of its mean queue length). The stationary vector is
            % then a signed measure whose aggregates over each phase block are
            % still the exact probabilities. See sn.isph and
            % _kb/04-networkstruct.md.
            if rate_a(ia)~=0 && (rate_a(ia)>0 || hasMEproc)
                % SPN code:
                %if rate_a(ia)>0 || modes_a(ia) > 0
                node_p = sync{a}.passive{1}.node;
                if node_p ~= local
                    % Skip if the active transition hash was not found
                    if new_state_a(ia) == -1
                        continue
                    end
                    state_p = state(sn.nodeToStateful(node_p));
                    class_p = sync{a}.passive{1}.class;
                    event_p = sync{a}.passive{1}.event;

                    % SPN code:
                    %                     enabled = 0;
                    %                     if ia <= length(trans_a)
                    %                         % check if other input places of the transition contains as many token as the multiplicity of the input arcs
                    %                         tr = trans_a(ia);
                    %                         mode = modes_a(ia);
                    %                         bmatrix = sn.varsparam{tr}.back(:,mode);
                    %                         inmatrix = sn.varsparam{tr}.inh(:,mode);
                    %                         enabled = all(state_pn >= bmatrix' & ~any(inmatrix'>0 & inmatrix' <= state_pn));
                    %                     end

                    %prob_sync_p = sync{a}.passive{1}.prob(state_a, state_p)
                    %if prob_sync_p > 0
                    %% debugging block
                    %if options.verbose == 2
                    %    line_printf('---\n');
                    %    sync{a}.active{1}.print,
                    %    sync{a}.passive{1}.print
                    %end
                    %%
                    if node_p == node_a %self-loop
                        [new_state_p, ~, outprob_p, start_p, preempt_p] = State.afterEventHashed( sn, node_p, new_state_a(ia), event_p, class_p);
                    else % departure
                        [new_state_p, ~, outprob_p, start_p, preempt_p] = State.afterEventHashed( sn, node_p, state_p, event_p, class_p);
                    end
                    %  SPN code:
                    %                    if node_p == node_a %self-loop
                    %                         [new_state_p, ~, outprob_p, trans_p, modes_p] = State.afterEventHashed( qn, node_p, new_state_a(ia), event_p, class_p);
                    %                     else % departure
                    %                         [new_state_p, ~, outprob_p, trans_p, modes_p] = State.afterEventHashed( qn, node_p, state_p, event_p, class_p);
                    %                     end
                    for ip=1:size(new_state_p,1)
                        if node_p ~= local
                            if new_state_p ~= -1
                                if sn.isstatedep(node_a,3)
                                    newStateCell = stateCell;
                                    newStateCell{sn.nodeToStateful(node_a)} = sn.space{sn.nodeToStateful(node_a)}(new_state_a(ia),:);
                                    newStateCell{sn.nodeToStateful(node_p)} = sn.space{sn.nodeToStateful(node_p)}(new_state_p(ip),:);
                                    prob_sync_p = sync{a}.passive{1}.prob(stateCell, newStateCell) * outprob_p(ip); %state-dependent
                                else
                                    prob_sync_p = sync{a}.passive{1}.prob * outprob_p(ip);
                                end
                            else
                                prob_sync_p = 0;
                            end
                        end
                        if ~isempty(new_state_a(ia))
                            if node_p == local % local action
                                new_state = state;
                                new_state(sn.nodeToStateful(node_a)) = new_state_a(ia);
                                prob_sync_p = outprob_p(ip);
                            elseif ~isempty(new_state_p)
                                new_state = state;
                                new_state(sn.nodeToStateful(node_a)) = new_state_a(ia);
                                new_state(sn.nodeToStateful(node_p)) = new_state_p(ip);
                            end
                            % SPN code:
                            %       if enabled
                            %                                 ns = find(ismember(SSh(:,[sn.nodeToStateful(node_a),sn.nodeToStateful(node_p)]),[new_state_a(ia),new_state_p(ip)],'rows'));
                            %                                 for ins=1:length(ns)
                            %                                    if ns(ins) > 0 && ~isempty(trans_p)
                            %                                         tr = trans_p(ip);
                            %                                         mode = modes_p(ip);
                            %                                         bmatrix = sn.varsparam{tr}.back(:,mode);
                            %                                         fmatrix = sn.varsparam{tr}.forw(:,mode);
                            %                                         cmatrix = fmatrix - bmatrix;
                            %                                         if isequal(state_pn + cmatrix',SS(ns(ins),3:end))
                            %                                             [ex_a,seq_a] = ST.search(state_pn');
                            %                                             [ex_p,seq_p] = ST.search(SS(ns(ins),3:end)');
                            %                                             if ex_a && ex_p && edges(seq_a, seq_p)
                            %                                                 Adj_m(s, ns(ins)) = modes_p(ip);
                            %                                                 Adj_t(s, ns(ins)) = trans_p(ip);
                            % %                                                 s,ns(ins)
                            %                                                 if ~isnan(rate_a(ia))
                            %                                                     if node_p < local && ~csmask(class_a, class_p) && rate_a(ia) * prob_sync_p >0 && (sn.nodetype(node_p)~=NodeType.Source)
                            %                                                         error('Error: state-dependent routing at node %d (%s) violates the class switching mask (node %d -> node %d, class %d -> class %d).', node_a, sn.nodenames{node_a}, node_a, node_p, class_a, class_p);
                            %                                                     end
                            %                                                     if size(Dfilt{a}) >= [s,ns(ins)] % check needed as D{a} is a sparse matrix
                            %                                                         Dfilt{a}(s,ns(ins)) = Dfilt{a}(s,ns(ins)) + rate_a(ia) * prob_sync_p;
                            %                                                     else
                            %                                                         Dfilt{a}(s,ns(ins)) = rate_a(ia) * prob_sync_p;
                            %                                                     end
                            %                                                 end
                            %                                             end
                            %                                         end
                            %                                    end
                            %                                 end
                            %            else

                            ns = matchrow(stateSpaceHashed, new_state);
                            if ns>0
                                if ~isnan(rate_a)
                                    if node_p < local && ~csmask(class_a, class_p) && rate_a(ia) * prob_sync_p >0 && (sn.nodetype(node_p)~=NodeType.Source)
                                        line_error(mfilename,sprintf('Error: state-dependent routing at node %d (%s) violates the class switching mask (node %s -> node %s, class %s -> class %s).', node_a, sn.nodenames{node_a}, sn.nodenames{node_a}, sn.nodenames{node_p}, sn.classnames{class_a}, sn.classnames{class_p}));
                                    end
                                    if size(Dfilt{a}) >= [s,ns] % check needed as D{a} is a sparse matrix
                                        Dfilt{a}(s,ns) = Dfilt{a}(s,ns) + rate_a(ia) * prob_sync_p;
                                    else
                                        Dfilt{a}(s,ns) = rate_a(ia) * prob_sync_p;
                                    end
                                    % Both halves of the synchronization are tagged: a
                                    % DEP promotes at the sender while the paired ARV
                                    % starts or preempts at the receiver, and each is
                                    % annotated on its own node's successor row.
                                    w_aux = rate_a(ia) * prob_sync_p;
                                    DfiltAux = solver_ctmc_auxfilt_add(DfiltAux, sn, node_a, s, ns, w_aux, start_a, preempt_a, ia);
                                    DfiltAux = solver_ctmc_auxfilt_add(DfiltAux, sn, node_p, s, ns, w_aux, start_p, preempt_p, ip);
                                end
                            end
                            % SPN code:
                            %                            end
                        end
                    end
                    % see _kb/06-solver-catalog.md (True BAS blocking) for rationale
                    R2 = sn.nclasses;
                    if event_a == EventType.DEP && rate_a(ia) > 0 ...
                            && ~isempty(sn.isbasblocking) && numel(sn.isbasblocking) >= node_a ...
                            && sn.isbasblocking(node_a) == 1 ...
                            && all(new_state_p(:) == -1)
                        isfA = sn.nodeToStateful(node_a);
                        curVecA = sn.space{isfA}(state(isfA),:);
                        if curVecA(end) == 0
                            blockedVec = curVecA; blockedVec(end) = 1;
                            blockedIdx = matchrow(sn.space{isfA}, blockedVec);
                            if blockedIdx > 0
                                new_state_b = state;
                                new_state_b(isfA) = blockedIdx;
                                nsb = matchrow(stateSpaceHashed, new_state_b);
                                if nsb > 0
                                    basBlockQ(s,nsb) = basBlockQ(s,nsb) + rate_a(ia);
                                end
                            end
                        end
                    end
                else % node_p == local
                    if ~isempty(new_state_a(ia))
                        new_state = state;
                        new_state(sn.nodeToStateful(node_a)) = new_state_a(ia);
                        prob_sync_p = 1;
                        ns = matchrow(stateSpaceHashed, new_state);
                        if ns>0
                            if ~isnan(rate_a)
                                if size(Dfilt{a}) >= [s,ns] % needed for sparse matrix
                                    Dfilt{a}(s,ns) = Dfilt{a}(s,ns) + rate_a(ia) * prob_sync_p;
                                else
                                    Dfilt{a}(s,ns) = rate_a(ia) * prob_sync_p;
                                end
                                % local action: only the active node can tag
                                DfiltAux = solver_ctmc_auxfilt_add(DfiltAux, sn, node_a, s, ns, rate_a(ia) * prob_sync_p, start_a, preempt_a, ia);
                            end
                        end
                    end
                end
            end
        end
    end
end
end % if fcrWaitq

% see _kb/06-solver-catalog.md (Vanishing states) for rationale
isfjaug = isfield(sn,'fjsync') && ~isempty(sn.fjsync);
immAction = false(1,A);
for a=1:A
    nt_a = sn.nodetype(sync{a}.active{1}.node);
    immAction(a) = (nt_a == NodeType.Router) || (nt_a == NodeType.Fork) || (isfjaug && nt_a == NodeType.Join);
end
for a=1:A
    Q = Q + Dfilt{a};
    if immAction(a)
        Qimm = Qimm + Dfilt{a};
    end
end
% Fold in true-BAS become-blocked transitions (not counted as departures).
Q = Q + basBlockQ;

%% for all global synchronizations (SPN support)
if isfield(sn, 'gsync') && ~isempty(sn.gsync)
    gsyncEvents = sn.gsync;
    G = length(gsyncEvents);

    % Track FIRE completion rates for arvRates/depRates
    Dfilt_gsync_comp = cell(1, G);
    for g = 1:G
        Dfilt_gsync_comp{g} = sparse(size(stateSpaceHashed,1), size(stateSpaceHashed,1));
    end
    immGsync = false(1,G);

    for g = 1:G
        gind = gsyncEvents{g}.active{1}.node;
        isf_transition = sn.nodeToStateful(gind);
        nmodes_g = sn.nodeparam{gind}.nmodes;
        % ENABLE phase moves and firings of a TimingStrategy.IMMEDIATE mode are
        % the two gsync sources emitted at the GlobalConstants.Immediate scale.
        if gsyncEvents{g}.active{1}.event == EventType.ENABLE
            immGsync(g) = true;
        elseif gsyncEvents{g}.active{1}.event == EventType.FIRE
            mode_g = gsyncEvents{g}.active{1}.mode;
            immGsync(g) = isfield(sn.nodeparam{gind},'timing') && ~isempty(sn.nodeparam{gind}.timing) ...
                && mode_g <= numel(sn.nodeparam{gind}.timing) ...
                && sn.nodeparam{gind}.timing(mode_g) == TimingStrategy.IMMEDIATE;
        end

        for s = 1:size(stateSpaceHashed, 1)
            state = stateSpaceHashed(s, :);

            % Build glspace cell array from hashed state
            glspace = cell(nstateful, 1);
            for isf = 1:nstateful
                glspace{isf} = sn.space{isf}(state(isf), :);
            end

            % Process event (both ENABLE and FIRE)
            [outglspace, outrate, outprob, outcomp] = State.afterGlobalEvent(sn, gind, glspace, gsyncEvents{g}, false);

            if isempty(outrate)
                continue;
            end

            for io = 1:length(outrate)
                if outrate(io) == 0
                    continue;
                end

                % Build new hashed state
                new_state = state;

                % Hash transition's new state
                if size(outglspace{isf_transition}, 1) < io
                    continue;
                end
                trans_state = outglspace{isf_transition}(io, :);
                hash_t = matchrow(sn.space{isf_transition}, trans_state);
                if hash_t <= 0, continue; end
                new_state(isf_transition) = hash_t;

                is_comp = false;
                % see _kb/06-solver-catalog.md (CTMC section, SPN FIRE completion flag) for rationale
                if gsyncEvents{g}.active{1}.event == EventType.FIRE
                    if ~isempty(outcomp) && io <= numel(outcomp)
                        is_comp = outcomp(io);
                    end
                    for isf = 1:nstateful
                        if isf ~= isf_transition && ~isequal(glspace{isf}, outglspace{isf})
                            hash_p = matchrow(sn.space{isf}, outglspace{isf});
                            if hash_p <= 0, continue; end
                            new_state(isf) = hash_p;
                        end
                    end
                end
                % For ENABLE events: only Transition state changes, Places unchanged

                ns = matchrow(stateSpaceHashed, new_state);
                if ns > 0
                    prob_val = 1;
                    if ~isempty(outprob) && io <= length(outprob)
                        prob_val = outprob(io);
                    end
                    rate_val = outrate(io) * prob_val;
                    Q(s, ns) = Q(s, ns) + rate_val;
                    if immGsync(g)
                        Qimm(s, ns) = Qimm(s, ns) + rate_val;
                    end
                    if is_comp
                        Dfilt_gsync_comp{g}(s, ns) = Dfilt_gsync_comp{g}(s, ns) + rate_val;
                    end
                end
            end
        end
    end
end

%% for all fork firing synchronizations (native fork-join support)
FJ = 0;
if isfield(sn,'fjsync') && ~isempty(sn.fjsync)
    FJ = length(sn.fjsync);
    Dfilt_fjsync = cell(1,FJ);
    for k=1:FJ
        Dfilt_fjsync{k} = sparse(size(stateSpaceHashed,1), size(stateSpaceHashed,1));
    end
    for k=1:FJ
        for s=1:size(stateSpaceHashed,1)
            state = stateSpaceHashed(s,:);
            glspace = cell(nstateful,1);
            for isf=1:nstateful
                glspace{isf} = sn.space{isf}(state(isf),:);
            end
            [fjStates, fjrate, fjprob] = State.afterFJEvent(sn, sn.fjsync{k}, glspace, false);
            for io=1:length(fjStates)
                if fjprob(io) <= 0
                    continue
                end
                new_state = state;
                skip = false;
                for isf=1:nstateful
                    if ~isequal(glspace{isf}, fjStates{io}{isf})
                        newrow = fjStates{io}{isf};
                        if length(newrow) < size(sn.space{isf},2)
                            newrow = [zeros(1,size(sn.space{isf},2)-length(newrow)), newrow];
                        end
                        hash_p = matchrow(sn.space{isf}, newrow);
                        if hash_p <= 0
                            skip = true;
                            break
                        end
                        new_state(isf) = hash_p;
                    end
                end
                if skip
                    continue
                end
                ns = matchrow(stateSpaceHashed, new_state);
                if ns > 0
                    rate_val = fjrate(io) * fjprob(io);
                    Q(s,ns) = Q(s,ns) + rate_val;
                    Qimm(s,ns) = Qimm(s,ns) + rate_val;
                    Dfilt_fjsync{k}(s,ns) = Dfilt_fjsync{k}(s,ns) + rate_val;
                end
            end
        end
    end
end

%% vanishing-row purge
% see _kb/06-solver-catalog.md (Vanishing states) for rationale
immPurged = [];
if options.config.hide_immediate
    immPurged = ctmc_find_vanishing_states(sn, stateSpaceHashed, nclasses, nstateful, FJ);
    % see _kb/06-solver-catalog.md (Vanishing states) for rationale
    immRows = immPurged;
    if ~isempty(immRows)
        immGap = immRows(full(sum(Qimm(immRows,:),2)) <= 0);
        if ~isempty(immGap)
            line_warning_always(mfilename, 'CTMC: %d vanishing state(s) have no immediate outgoing arc; the vanishing predicate and the immediate-arc tagging disagree, so those rows keep their timed arcs.', numel(immGap));
            immRows = setdiff(immRows, immGap);
        end
    end
    if ~isempty(immRows)
        Q(immRows,:) = Qimm(immRows,:);
        for a=1:A
            if ~immAction(a)
                Dfilt{a}(immRows,:) = 0;
            end
        end
        if exist('Dfilt_gsync_comp','var')
            for g=1:numel(Dfilt_gsync_comp)
                if ~immGsync(g)
                    Dfilt_gsync_comp{g}(immRows,:) = 0;
                end
            end
        end
    end
end

Q = Q - diag(diag(Q));
%SolverCTMC.printInfGen(Q,stateSpace)
%%
arvRates = zeros(size(stateSpaceHashed,1),nstateful,nclasses);
depRates = zeros(size(stateSpaceHashed,1),nstateful,nclasses);
for a=1:A
    % active
    node_a = sync{a}.active{1}.node;
    class_a = sync{a}.active{1}.class;
    event_a = sync{a}.active{1}.event;
    % passive
    node_p = sync{a}.passive{1}.node;
    class_p = sync{a}.passive{1}.class;
    if event_a == EventType.DEP
        node_a_sf = sn.nodeToStateful(node_a);
        node_p_sf = sn.nodeToStateful(node_p);
        for s=1:size(stateSpaceHashed,1)
            depRates(s,node_a_sf,class_a) = depRates(s,node_a_sf,class_a) + sum(Dfilt{a}(s,:));
            arvRates(s,node_p_sf,class_p) = arvRates(s,node_p_sf,class_p) + sum(Dfilt{a}(s,:));
        end
    end
end

%% Compute arrival/departure rates for gsync FIRE completion events
if isfield(sn, 'gsync') && ~isempty(sn.gsync)
    gsyncEvents = sn.gsync;
    G = length(gsyncEvents);
    for g = 1:G
        if gsyncEvents{g}.active{1}.event == EventType.FIRE
            for j = 1:length(gsyncEvents{g}.passive)
                pev = gsyncEvents{g}.passive{j};
                [pev_node, pev_class] = ind2sub([sn.nnodes, nclasses], pev.node);
                if pev_node > sn.nnodes || ~sn.isstateful(pev_node)
                    continue;
                end
                pev_isf = sn.nodeToStateful(pev_node);
                if pev.event == EventType.PRE
                    for s = 1:size(stateSpaceHashed, 1)
                        depRates(s, pev_isf, pev_class) = depRates(s, pev_isf, pev_class) + sum(Dfilt_gsync_comp{g}(s,:));
                    end
                elseif pev.event == EventType.POST
                    for s = 1:size(stateSpaceHashed, 1)
                        arvRates(s, pev_isf, pev_class) = arvRates(s, pev_isf, pev_class) + sum(Dfilt_gsync_comp{g}(s,:));
                    end
                end
            end
        end
    end
end

%% Compute arrival/departure rates for fork firing synchronizations
if FJ > 0
    for k=1:FJ
        fjentry = sn.fjsync{k};
        isf_fork = sn.nodeToStateful(fjentry.fork);
        for s=1:size(stateSpaceHashed,1)
            rowsum = sum(Dfilt_fjsync{k}(s,:));
            if rowsum > 0
                depRates(s, isf_fork, fjentry.class) = depRates(s, isf_fork, fjentry.class) + rowsum;
                for b=1:length(fjentry.branchheads)
                    isf_bh = sn.nodeToStateful(fjentry.branchheads(b));
                    arvRates(s, isf_bh, fjentry.auxclasses(b)) = arvRates(s, isf_bh, fjentry.auxclasses(b)) + rowsum;
                end
            end
        end
    end
end

zero_row = find(sum(Q,2)==0);
zero_col = find(sum(Q,1)==0);

%%
% in case the last column of Q represent a state for a transient class, it
% is possible that no transitions go back to it, although it is valid for
% the system to be initialized in that state. So we need to fill-in the
% zeros at the end.
Q(:,end+1:end+(size(Q,1)-size(Q,2)))=0;
Q(zero_row,zero_row) = -eye(length(zero_row)); % can this be replaced by []?
Q(zero_col,zero_col) = -eye(length(zero_col));
for a=1:A
    Dfilt{a}(:,end+1:end+(size(Dfilt{a},1)-size(Dfilt{a},2)))=0;
end
for i=1:size(DfiltAux.start,1)
    for r=1:size(DfiltAux.start,2)
        DfiltAux.start{i,r}(:,end+1:end+(size(DfiltAux.start{i,r},1)-size(DfiltAux.start{i,r},2)))=0;
        DfiltAux.preempt{i,r}(:,end+1:end+(size(DfiltAux.preempt{i,r},1)-size(DfiltAux.preempt{i,r},2)))=0;
    end
end

if options.verbose == VerboseLevel.DEBUG || GlobalConstants.Verbose == VerboseLevel.DEBUG
    SolverCTMC.printInfGen(Q,stateSpace);
end
Q = ctmc_makeinfgen(Q);

%% drop states unreachable from the initial state
% see _kb/06-solver-catalog.md (CTMC section, unreachable-state pruning) for rationale
if ~isempty(sn.state) && all(~cellfun(@isempty, sn.state))
    % The per-station initial rows carry only as many buffer slots as the
    % initial population needs, while the enumerated local space is sized for
    % the full capacity. Left-pad each row to its space width (empty buffer
    % slots pad the left, so the server-phase and local-variable tail stays
    % aligned) before matching; without this the lookup fails and the pruning
    % below is silently skipped, leaving any enumerated-but-unreachable
    % absorbing state to break ctmc_solve.
    initRow = [];
    for isf = 1:sn.nstateful
        row_isf = sn.state{isf};
        row_isf = row_isf(1,:);
        w_isf = size(sn.space{isf},2);
        if numel(row_isf) < w_isf
            row_isf = [zeros(1, w_isf-numel(row_isf)), row_isf];
        end
        initRow = [initRow, row_isf]; %#ok<AGROW>
    end
    initState = matchrow(stateSpace, initRow);
    if initState > 0
        % any nonzero off-diagonal is an arc: an ME embeds with negative ones
        adj = abs(Q - diag(diag(Q))) > GlobalConstants.ArcTol;
        reach = false(size(Q,1),1);
        reach(initState) = true;
        frontier = initState;
        while ~isempty(frontier)
            nxt = find(any(adj(frontier,:),1))';
            nxt = nxt(~reach(nxt));
            reach(nxt) = true;
            frontier = nxt;
        end
        % Checked BEFORE pruning: afterwards the retained chain is irreducible by
        % construction and the ambiguity is no longer observable. A state that
        % cannot reach the initial class must lead into some other closed class,
        % so the backward closure of `reach` is exactly the test, at the cost of
        % one transposed sweep rather than a full SCC decomposition.
        if ~all(reach)
            feeds = reach;
            frontierB = find(reach)';
            while ~isempty(frontierB)
                prv = find(any(adj(:,frontierB),2))';
                prv = prv(~feeds(prv));
                feeds(prv) = true;
                frontierB = prv;
            end
            if ~all(feeds)
                line_warning_always(mfilename, ...
                    ['CTMC: the generator has more than one closed communicating class (%d of %d ' ...
                    'states cannot reach the one containing the initial state), so its stationary ' ...
                    'distribution is not unique. Results are those of the class the initial state ' ...
                    'selects, and depend on it: an order-preserving discipline under routing that ' ...
                    'admits no overtaking freezes the relative service order at time zero, so ' ...
                    'declaring the classes in a different order gives a different, equally valid ' ...
                    'answer. Product-form results (MVA, NC) normalize over the whole space instead ' ...
                    'and will not agree.'], sum(~feeds), numel(feeds));
            end
        end
        if ~all(reach)
            line_debug('CTMC: %d of %d states unreachable from the initial state, dropped', ...
                sum(~reach), numel(reach));
            keep = find(reach);
            immPurged = []; % state indices shift, the predicate must be re-evaluated
            Q = Q(keep,keep);
            Q = ctmc_makeinfgen(Q);
            stateSpace = stateSpace(keep,:);
            stateSpaceAggr = stateSpaceAggr(keep,:);
            stateSpaceHashed = stateSpaceHashed(keep,:);
            arvRates = arvRates(keep,:,:);
            depRates = depRates(keep,:,:);
            for a=1:A
                Dfilt{a} = Dfilt{a}(keep,keep);
            end
            for i=1:size(DfiltAux.start,1)
                for r=1:size(DfiltAux.start,2)
                    DfiltAux.start{i,r} = DfiltAux.start{i,r}(keep,keep);
                    DfiltAux.preempt{i,r} = DfiltAux.preempt{i,r}(keep,keep);
                end
            end
            for k=1:FJ
                Dfilt_fjsync{k} = Dfilt_fjsync{k}(keep,keep);
            end
            if exist('Dfilt_gsync_comp','var')
                for g=1:numel(Dfilt_gsync_comp)
                    Dfilt_gsync_comp{g} = Dfilt_gsync_comp{g}(keep,keep);
                end
            end
        end
    end
end

%% now remove immediate transitions
% we first determine states in stateful nodes where there is an immediate
% job in the node

if options.config.hide_immediate % if want to remove immediate transitions
    % see _kb/06-solver-catalog.md (Vanishing states, Design Y positive list) for rationale
    if isempty(immPurged)
        imm = ctmc_find_vanishing_states(sn, stateSpaceHashed, nclasses, nstateful, FJ);
    else
        imm = immPurged; % already evaluated on this state space by the purge above
    end
    nonimm = setdiff(1:size(Q,1),imm);
    stateSpace(imm,:) = [];
    stateSpaceAggr(imm,:) = [];
    %    full(Q)
    [Q,~,Q12,~,Q22] = ctmc_stochcomp(Q, nonimm);
    %    full(Q)
    if FJ > 0 || ~isempty(imm)
        % see _kb/06-solver-catalog.md (Vanishing states, rate complement) for rationale
        arvRates = zeros(length(nonimm),nstateful,nclasses);
        depRates = zeros(length(nonimm),nstateful,nclasses);
        for a=1:A
            % active
            node_a = sync{a}.active{1}.node;
            class_a = sync{a}.active{1}.class;
            event_a = sync{a}.active{1}.event;
            % passive
            node_p = sync{a}.passive{1}.node;
            class_p = sync{a}.passive{1}.class;
            if event_a == EventType.DEP
                node_a_sf = sn.nodeToStateful(node_a);
                node_p_sf = sn.nodeToStateful(node_p);
                % see _kb/06-solver-catalog.md (Vanishing states, rate complement) for rationale
                if immAction(a) && isfjaug && sn.nodetype(node_a) == NodeType.Join
                    r_a = solver_ctmc_ratecomplement(Dfilt{a}, nonimm, imm, Q12, Q22);
                else
                    r_a = full(sum(Dfilt{a}(nonimm,:),2));
                end
                depRates(:,node_a_sf,class_a) = depRates(:,node_a_sf,class_a) + r_a;
                arvRates(:,node_p_sf,class_p) = arvRates(:,node_p_sf,class_p) + r_a;
            end
        end
        % see _kb/06-solver-catalog.md (Vanishing states, rate complement) for rationale
        if isfield(sn, 'gsync') && ~isempty(sn.gsync) && exist('Dfilt_gsync_comp','var')
            gsyncEvents_rc = sn.gsync;
            for g = 1:length(gsyncEvents_rc)
                if gsyncEvents_rc{g}.active{1}.event ~= EventType.FIRE
                    continue;
                end
                r_g = solver_ctmc_ratecomplement(Dfilt_gsync_comp{g}, nonimm, imm, Q12, Q22);
                for j = 1:length(gsyncEvents_rc{g}.passive)
                    pev = gsyncEvents_rc{g}.passive{j};
                    [pev_node, pev_class] = ind2sub([sn.nnodes, nclasses], pev.node);
                    if pev_node > sn.nnodes || ~sn.isstateful(pev_node)
                        continue;
                    end
                    pev_isf = sn.nodeToStateful(pev_node);
                    if pev.event == EventType.PRE
                        depRates(:, pev_isf, pev_class) = depRates(:, pev_isf, pev_class) + r_g;
                    elseif pev.event == EventType.POST
                        arvRates(:, pev_isf, pev_class) = arvRates(:, pev_isf, pev_class) + r_g;
                    end
                end
            end
        end
        % fork firings: departure of the parent class at the Fork, one
        % sibling arrival per branch head in the tag's auxiliary classes
        for k=1:FJ
            fjentry = sn.fjsync{k};
            isf_fork = sn.nodeToStateful(fjentry.fork);
            r_k = solver_ctmc_ratecomplement(Dfilt_fjsync{k}, nonimm, imm, Q12, Q22);
            depRates(:,isf_fork,fjentry.class) = depRates(:,isf_fork,fjentry.class) + r_k;
            for b=1:length(fjentry.branchheads)
                isf_bh = sn.nodeToStateful(fjentry.branchheads(b));
                arvRates(:,isf_bh,fjentry.auxclasses(b)) = arvRates(:,isf_bh,fjentry.auxclasses(b)) + r_k;
            end
        end
    end
    for a=1:A
        % stochastic complement for action a
        Q21a = Dfilt{a}(imm,nonimm);
        Ta = (-Q22) \ Q21a;
        Ta = Q12*Ta;
        Dfilt{a} = Dfilt{a}(nonimm,nonimm)+Ta;
    end
    % The derived filtrations are complemented exactly like Dfilt: a service
    % start that lands on a vanishing state would otherwise be dropped, and
    % the START rate would silently undercount by that path.
    for i=1:size(DfiltAux.start,1)
        for r=1:size(DfiltAux.start,2)
            Tsa = Q12*((-Q22) \ DfiltAux.start{i,r}(imm,nonimm));
            DfiltAux.start{i,r} = DfiltAux.start{i,r}(nonimm,nonimm)+Tsa;
            Tpa = Q12*((-Q22) \ DfiltAux.preempt{i,r}(imm,nonimm));
            DfiltAux.preempt{i,r} = DfiltAux.preempt{i,r}(nonimm,nonimm)+Tpa;
        end
    end
    % recompute arvRates and depRates
    %     arvRates = zeros(size(stateSpace,1),nstateful,nclasses);
    %     depRates = zeros(size(stateSpace,1),nstateful,nclasses);
    %     for a=1:A
    %         % active
    %         node_a = sync{a}.active{1}.node;
    %         class_a = sync{a}.active{1}.class;
    %         event_a = sync{a}.active{1}.event;
    %         % passive
    %         node_p = sync{a}.passive{1}.node;
    %         class_p = sync{a}.passive{1}.class;
    %         if event_a == EventType.DEP
    %             node_a_sf = sn.nodeToStateful(node_a);
    %             node_p_sf = sn.nodeToStateful(node_p);
    %             for s=1:size(stateSpace,1)
    %                 depRates(s,node_a_sf,class_a) = depRates(s,node_a_sf,class_a) + sum(Dfilt{a}(s,:));
    %                 arvRates(s,node_p_sf,class_p) = arvRates(s,node_p_sf,class_p) + sum(Dfilt{a}(s,:));
    %             end
    %         end
    %     end
end
end

%% Local functions
function imm = ctmc_find_vanishing_states(sn, stateSpaceHashed, nclasses, nstateful, FJ)
% Indices of the vanishing (zero-sojourn) global states: Router/Fork
% pass-through occupancy, firable Join sibling sets, SPN markings from
% which an ENABLE event moves the Transition row, and markings enabling a
% TimingStrategy.IMMEDIATE mode. Extracted so the same predicate drives
% both the vanishing-row purge and the stochastic complementation below.
    isImmediatePassThrough = @(nt) (nt == NodeType.Router || nt == NodeType.Fork);

    imm = [];
    for ind = 1:sn.nnodes
        if sn.isstateful(ind) && ~sn.isstation(ind) && isImmediatePassThrough(sn.nodetype(ind))
            isf = sn.nodeToStateful(ind);
            imm_st = find(sum(sn.space{isf}(:,1:nclasses),2)>0);
            imm = [imm; find(arrayfun(@(a) any(a==imm_st),stateSpaceHashed(:,isf)))];
        end
    end
    % see _kb/06-solver-catalog.md (Vanishing states) for rationale
    if FJ > 0
        for ind = 1:sn.nnodes
            if sn.nodetype(ind) == NodeType.Join
                isf = sn.nodeToStateful(ind);
                firable_rows = [];
                origcl = sn.nodeparam{ind}.fj.origclasses;
                for row=1:size(sn.space{isf},1)
                    for r=origcl(:)'
                        [ospace_j] = State.afterEventJoin(sn, ind, sn.space{isf}(row,:), EventType.DEP, r, false, [], NaN);
                        if ~isempty(ospace_j)
                            firable_rows(end+1) = row; %#ok<AGROW>
                            break
                        end
                    end
                end
                if ~isempty(firable_rows)
                    imm = [imm; find(arrayfun(@(a) any(a==firable_rows),stateSpaceHashed(:,isf)))]; %#ok<AGROW>
                end
            end
        end
    end
    % Transition immediate states: states where any ENABLE event would change the state
    if isfield(sn, 'gsync') && ~isempty(sn.gsync)
        gsyncEvents_sc = sn.gsync;
        for s = 1:size(stateSpaceHashed, 1)
            if any(s == imm)
                continue; % already marked
            end
            state_sc = stateSpaceHashed(s, :);
            glspace_sc = cell(nstateful, 1);
            for isf = 1:nstateful
                glspace_sc{isf} = sn.space{isf}(state_sc(isf), :);
            end
            for g = 1:length(gsyncEvents_sc)
                if gsyncEvents_sc{g}.active{1}.event == EventType.ENABLE
                    gind_sc = gsyncEvents_sc{g}.active{1}.node;
                    isf_t_sc = sn.nodeToStateful(gind_sc);
                    orig_row_sc = glspace_sc{isf_t_sc};
                    [outgl_sc, outrate_sc, ~] = State.afterGlobalEvent(sn, gind_sc, glspace_sc, gsyncEvents_sc{g}, false);
                    % see _kb/06-solver-catalog.md (CTMC section, SPN FIRE completion flag) for rationale
                    rowchanged = false;
                    if ~isempty(outrate_sc)
                        og_sc = outgl_sc{isf_t_sc};
                        for io=1:size(og_sc,1)
                            if outrate_sc(io) > 0 && ~isequal(og_sc(io,:), orig_row_sc)
                                rowchanged = true;
                                break;
                            end
                        end
                    end
                    if rowchanged
                        imm = [imm; s]; %#ok<AGROW>
                        break;
                    end
                end
            end
        end
    end
    % see _kb/06-solver-catalog.md (Vanishing states) for rationale
    if isfield(sn, 'gsync') && ~isempty(sn.gsync)
        gsyncEvents_im = sn.gsync;
        for g = 1:length(gsyncEvents_im)
            if gsyncEvents_im{g}.active{1}.event ~= EventType.FIRE
                continue;
            end
            gind_im = gsyncEvents_im{g}.active{1}.node;
            mode_im = gsyncEvents_im{g}.active{1}.mode;
            if ~isfield(sn.nodeparam{gind_im}, 'timing') || isempty(sn.nodeparam{gind_im}.timing)
                continue;
            end
            if mode_im > numel(sn.nodeparam{gind_im}.timing) ...
                    || sn.nodeparam{gind_im}.timing(mode_im) ~= TimingStrategy.IMMEDIATE
                continue;
            end
            for s = 1:size(stateSpaceHashed, 1)
                if any(s == imm)
                    continue; % already marked
                end
                glspace_im = cell(nstateful, 1);
                for isf = 1:nstateful
                    glspace_im{isf} = sn.space{isf}(stateSpaceHashed(s, isf), :);
                end
                [~, outrate_im, ~] = State.afterGlobalEvent(sn, gind_im, glspace_im, gsyncEvents_im{g}, false);
                if ~isempty(outrate_im) && any(outrate_im > 0)
                    imm = [imm; s]; %#ok<AGROW>
                end
            end
        end
    end

    imm = unique(imm);
end
