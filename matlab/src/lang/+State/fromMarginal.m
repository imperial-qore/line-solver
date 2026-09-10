function space = fromMarginal(sn, ind, n, options)
% FROMMARGINAL Generate state space with specific marginal queue lengths
%
% @brief Creates state space where a specific node has given marginal queue lengths
% @param sn Network structure or Network object
% @param ind Node index for which to set marginal queue lengths
% @param n Vector of jobs per class at the specified node
% @param options Optional structure with configuration parameters
% @return space Generated state space satisfying marginal constraints
%
% This function generates all possible network states where the specified
% node (ind) has exactly n(r) jobs of class r, for all classes. It is
% essential for state space analysis and marginal probability computations.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<4 %~exist('options','var')
    options.force = false;
end
if isa(sn,'Network')
    sn=sn.getStruct();
end


% generate states such that the marginal queue-lengths are as in vector n
%  n(r): number of jobs at the station in class r
R = sn.nclasses;
S = sn.nservers;
state = [];
space = [];

% ind: node index
ist = sn.nodeToStation(ind);
%isf = sn.nodeToStateful(ind);

% Synchronous call (REPLY signal): the node holds one server per job that has
% left for its callee and is waiting for the reply. Those servers are not
% derivable from the marginal N, so enumerate the held-server counts here and
% build the rest of the state with the REMAINING servers -- with b servers held,
% only S-b jobs can be in service, a configuration the plain enumeration never
% produces. Recurse on a struct with the block cleared, then append its columns,
% which are the last ones in the local-variable layout (State.replyBlockInfo).
if isfield(sn,'replyblock') && ~isempty(sn.replyblock) && size(sn.replyblock,1) >= ind ...
        && any(sn.replyblock(ind,:) > 0)
    ist_r = sn.nodeToStation(ind);
    rclasses = find(sn.replyblock(ind,:) > 0);
    snb = sn;
    snb.replyblock(ind,:) = 0;
    snb.nvars(ind, (2*R+2):(3*R+1)) = 0;
    bspace = zeros(1,0);
    for r = rclasses
        bspace = State.cartesian(bspace, (0:S(ist_r))');
    end
    subspaces = cell(size(bspace,1),1);
    bkept = cell(size(bspace,1),1);
    maxw = 0;
    for bi = 1:size(bspace,1)
        b = bspace(bi,:);
        if sum(b) > S(ist_r)
            continue
        end
        snb.nservers(ist_r) = S(ist_r) - sum(b);
        if nargin < 4
            subspace = State.fromMarginal(snb, ind, n);
        else
            subspace = State.fromMarginal(snb, ind, n, options);
        end
        if isempty(subspace)
            continue
        end
        subspaces{bi} = subspace;
        bkept{bi} = b;
        maxw = max(maxw, size(subspace,2));
    end
    % Held servers push jobs into the buffer, so the sub-spaces have different
    % buffer widths. The buffer is RIGHT-aligned (empty slots pad the left), so
    % widen the narrow rows on the left before stacking them.
    space = [];
    for bi = 1:size(bspace,1)
        if isempty(subspaces{bi})
            continue
        end
        subspace = subspaces{bi};
        if size(subspace,2) < maxw
            subspace = [zeros(size(subspace,1), maxw-size(subspace,2)), subspace]; %#ok<AGROW>
        end
        space = [space; subspace, repmat(bkept{bi}, size(subspace,1), 1)]; %#ok<AGROW>
    end
    space = unique(space,'rows');
    space = space(end:-1:1,:);
    return
end

if isfield(sn,'isfjaugmented') && sn.isfjaugmented && sn.nodetype(ind) == NodeType.Join
    % FJ-augmented struct: the join state is the per-class count vector
    % of buffered jobs/siblings, deterministic given the marginals
    space = n(:)';
    return
end

if sn.isstation(ind) && any(sn.procid(ist,:)==ProcessType.MAP | sn.procid(sn.nodeToStation(ind),:)==ProcessType.MMPP2)
    if sn.sched(ist) ~= SchedStrategy.FCFS && sn.nodetype(ind) ~= NodeType.Source
        line_error(mfilename,'Non-FCFS MAP stations are not supported.')
    end
end

if sn.isstateful(ind) && ~sn.isstation(ind)
    if sn.nodetype(ind) == NodeType.Transition
        % Transition state format is per-mode, not per-class
        isf = sn.nodeToStateful(ind);
        if ~isempty(sn.space) && isf <= length(sn.space) && ~isempty(sn.space{isf})
            space = sn.space{isf};
        else
            % Generate initial state only (all servers idle)
            nmodes_t = sn.nodeparam{ind}.nmodes;
            nmodeservers_t = sn.nodeparam{ind}.nmodeservers;
            nmodeservers_t(isinf(nmodeservers_t)) = GlobalConstants.MaxInt();
            firingphases_t = sn.nodeparam{ind}.firingphases;
            firingphases_t(isnan(firingphases_t)) = 1;
            space = [nmodeservers_t, zeros(1, sum(firingphases_t)), zeros(size(nmodeservers_t))];
        end
        return
    end
    for r=1:R
        init_r = State.spaceClosedSingle(1,n(r));
        state = State.cartesian(state,init_r);
    end
    space = State.cartesian(space,state);
    return
end

phases = zeros(1,R);
for r=1:R
    if isempty(sn.proc{ist}{r})
        phases(r) = 0;
    elseif isfield(sn,'markidx') && ~isempty(sn.markidx) ...
            && ist <= size(sn.markidx,1) && sn.markidx(ist,r) > 1
        % Marked (MMAP) non-carrier class: modulating chain lives in the carrier's phase block -- see _kb/04-networkstruct.md
        phases(r) = 1;
    else
        phases(r) = length(sn.proc{ist}{r}{1});
    end
end
if (sn.sched(ist) ~= SchedStrategy.EXT) && any(n>sn.classcap(ist,:))
    return
end

% generate local-state space
switch sn.nodetype(ind)
    case {NodeType.Queue, NodeType.Delay, NodeType.Source, NodeType.Place}
        isRetrialStation = isfield(sn,'retrialProc') && ~isempty(sn.retrialProc) ...
            && ist > 0 && any(~cellfun(@isempty, sn.retrialProc(ist,:)));
        if isRetrialStation
            % Retrial station: enumerate every (in-service,orbit) split, incl. idle-server states -- see _kb/04-networkstruct.md
            r = find(n>0, 1);
            if isempty(r)
                space = zeros(1, sum(phases));
            else
                maxorbit = n(r);
                space = [];
                for csrv = 0:min(n(r), S(ist))
                    orbit = n(r) - csrv;
                    buf = [zeros(1, maxorbit-orbit), r*ones(1, orbit)];
                    srv = [];
                    for cls = 1:R
                        if cls == r
                            srv = State.cartesian(srv, State.spaceClosedSingle(phases(cls), csrv));
                        else
                            srv = State.cartesian(srv, State.spaceClosedSingle(phases(cls), 0));
                        end
                    end
                    space = [space; repmat(buf, size(srv,1), 1), srv];
                end
            end
        else
        switch sn.sched(ist)
            case SchedStrategy.EXT
                for r=1:R
                    if ~isempty(sn.proc) && ~isempty(sn.proc{ist}{r}) && any(any(isnan(sn.proc{ist}{r}{1}))) % disabled
                        init_r = 0*ones(1,phases(r));
                    elseif isfield(sn,'markidx') && ~isempty(sn.markidx) ...
                            && ist <= size(sn.markidx,1) && sn.markidx(ist,r) > 1
                        % marked non-carrier class: single always-zero column
                        init_r = 0;
                    else
                        init_r = State.spaceClosedSingle(phases(r),1);
                    end
                    state = State.cartesian(state,init_r);
                end
                space = State.cartesian(space,state); %server part
                space = [Inf*ones(size(space,1),1),space]; % attach infinite buffer before servers
            case {SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO, SchedStrategy.LPS}
                % in these policies we only track the jobs in the servers
                for r=1:R
                    init_r = State.spaceClosedSingle(phases(r),n(r));
                    state = State.cartesian(state,init_r);
                end
                space = State.cartesian(space,state);
            case SchedStrategy.POLLING
                % Un-ordered per-class buffers and a single server, as in SIRO,
                % but NOT work-conserving over the station: the server may be
                % idle while jobs wait, because it is walking towards a buffer
                % (a switchover) or is parked. Both the empty-facility and the
                % one-job-in-service configurations are therefore enumerated;
                % the controller columns appended after the routing variables
                % below discard the combinations the discipline cannot occupy.
                if S(ist) ~= 1
                    line_error(mfilename,'Polling stations must have a single server.');
                end
                space = [n, zeros(1, sum(phases))]; % service facility empty
                for p=1:R
                    if n(p) > 0
                        bufp = n; bufp(p) = bufp(p) - 1;
                        srvp = [];
                        for cls=1:R
                            srvp = State.cartesian(srvp, State.spaceClosedSingle(phases(cls), double(cls==p)));
                        end
                        space = [space; repmat(bufp, size(srvp,1), 1), srvp];
                    end
                end
            case {SchedStrategy.SIRO, SchedStrategy.LEPT, SchedStrategy.SEPT, SchedStrategy.SRPT, SchedStrategy.SRPTPRIO, SchedStrategy.SETF, SchedStrategy.FSP}
                % Unordered buffer + in-server jobs -- see _kb/04-networkstruct.md
                % build list of job classes in the node, with repetition
                if sum(n) <= S(ist)
                    for r=1:R
                        init_r = State.spaceClosedSingle(phases(r),n(r));
                        state = State.cartesian(state,init_r);
                    end
                    space = State.cartesian(space,[zeros(size(state,1),R),state]);
                else
                    si = multichoosecon(n,S(ist)); % jobs of class r that are running
                    mi_buf = repmat(n,size(si,1),1) - si; % jobs of class r in buffer
                    for k=1:size(si,1)
                        % determine number of classes r jobs running in phase j
                        kstate=[];
                        for r=1:R
                            init_r = State.spaceClosedSingle(phases(r),si(k,r));
                            kstate = State.cartesian(kstate,init_r);
                        end
                        state = [repmat(mi_buf(k,:),size(kstate,1),1), kstate];
                        space = [space; state];
                    end
                end
            case {SchedStrategy.FCFSPI, SchedStrategy.FCFSPR, SchedStrategy.FCFSPIPRIO, SchedStrategy.FCFSPRPRIO, SchedStrategy.LCFSPI, SchedStrategy.LCFSPR, SchedStrategy.LCFSPIPRIO, SchedStrategy.LCFSPRPRIO, SchedStrategy.EDF}
                sizeEstimator = multinomialln(n) - gammaln(sum(n)) + gammaln(1+sn.cap(ist));
                sizeEstimator = round(sizeEstimator/log(10));
                if sizeEstimator > 3
                    if ~isfield(options,'force') || options.force == false
                        %line_warning(mfilename,sprintf('Marginal state space size is in the order of thousands of states. Computation may be slow.',sizeEstimator));
                    end
                end

                if sum(n) == 0
                    % Preempt-resume buffer is (class,phase) pairs; empty buffer must be even-width -- see _kb/04-networkstruct.md
                    space = zeros(1,2+sum(phases));
                    space = sub_routevars(sn, ind, R, space);
                    space = sub_trailingvars(sn, ind, R, space, n);
                    return
                end
                % Ordered buffer + in-server jobs -- see _kb/04-networkstruct.md

                % build list of job classes in the node, with repetition
                vi = [];
                for r=1:R
                    if n(r)>0
                        vi=[vi, r*ones(1,n(r))];
                    end
                end

                % gen permutation of their positions in the waiting buffer
                mi = multiset_perms(vi);
                % now generate server states
                if isempty(mi)
                    mi_buf = zeros(1,max(0,sum(n)-S(ist)));
                    state = zeros(1,R);
                    state = State.cartesian(state,[mi_buf,state]);
                else
                    mi = mi(:,(end-min(sum(n),sn.cap(ist))+1):end); % n(r) may count more than once elements within the same chain
                    mi = unique(mi,'rows');
                    % mi_buf: class of job in buffer position i (0=empty)
                    mi_buf = [zeros(size(mi,1),min(sum(n),sn.cap(ist))-S(ist)-size(mi(:,1:end-S(ist)),2)), mi(:,1:end-S(ist))];
                    if isempty(mi_buf)
                        mi_buf = zeros(size(mi,1),1);
                    end
                    mi_buf_kstate = [];
                    %if mi_buf(1)>0
                    % generate job phases for all buffer states
                    %for k=1:size(mi_buf,1)
                    %    mi_buf_kstate(end+1:end+size(bkstate,1),1:size(bkstate,2)) = bkstate;
                    %end
                    %end
                    % mi_srv: class of job running in server i
                    mi_srv = mi(:,max(size(mi,2)-S(ist)+1,1):end);
                    % si: number of class r jobs that are running
                    si =[];
                    for k=1:size(mi_srv,1)                
                        %si(k,1:R) = hist(mi_srv(k,:),1:R); % deprecated
                        si(k, 1:R) = histcounts(mi_srv(k, :), 1:(R+1));
                    end
                    %si = unique(si,'rows');
                    for k=1:size(si,1)
                        % determine number of class r jobs running in phase
                        % j in server state mi_srv(kjs,:) and build
                        % state
                        kstate=[];
                        for r=1:R
                            kstate = State.cartesian(kstate,State.spaceClosedSingle(phases(r),si(k,r)));
                        end
                        % generate job phases for all buffer states
                        bkstate = [];
                        for j=mi_buf(k,:) % for each job in the buffer
                            if j>0
                                bkstate = State.cartesian(bkstate,[1:phases(j)]');
                            else
                                bkstate = 0;
                            end
                        end
                        bufstate_tmp = State.cartesian(mi_buf(k,:), bkstate);
                        % here interleave positions of class and phases in
                        % buf
                        bufstate = zeros(size(bufstate_tmp));
                        bufstate(:,1:2:end)=bufstate_tmp(:,1:size(mi_buf,2));
                        bufstate(:,2:2:end)=bufstate_tmp(:,(size(mi_buf,2)+1):end);
                        state = [state; State.cartesian(bufstate, kstate)];
                    end
                end
                space = state;
            case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRIO, SchedStrategy.LCFS, SchedStrategy.LCFSPRIO, SchedStrategy.EDD}
                sizeEstimator = multinomialln(n) - gammaln(sum(n)) + gammaln(1+sn.cap(ist));
                sizeEstimator = round(sizeEstimator/log(10));
                if sizeEstimator > 3
                    if ~isfield(options,'force') || options.force == false
                        %line_warning(mfilename,sprintf('Marginal state space size is in the order of thousands of states. Computation may be slow.',sizeEstimator));
                    end
                end

                if sum(n) == 0
                    space = zeros(1,1+sum(phases));
                    if sn.nodetype(ind) ~= NodeType.Source
                        for r=1:R
                            switch sn.procid(sn.nodeToStation(ind),r)
                                case {ProcessType.MAP, ProcessType.MMPP2}
                                    space = State.cartesian(space, [1:sn.phases(ind,r)]');
                            end
                        end
                    end
                    % Routing vars precede the node block in the nvars layout.
                    space = sub_routevars(sn, ind, R, space);
                    space = sub_trailingvars(sn, ind, R, space, n);
                    return
                end
                % Ordered buffer + in-server jobs -- see _kb/04-networkstruct.md

                % build list of job classes in the node, with repetition
                vi = [];
                for r=1:R
                    if n(r)>0
                        vi=[vi, r*ones(1,n(r))];
                    end
                end

                % gen permutation of their positions in the waiting buffer
                mi = multiset_perms(vi);
                % now generate server states
                if isempty(mi)
                    mi_buf = zeros(1,max(0,sum(n)-S(ist)));
                    state = zeros(1,R);
                    state = State.cartesian(state,[mi_buf,state]);
                else
                    mi = mi(:,(end-min(sum(n),sn.cap(ist))+1):end); % n(r) may count more than once elements within the same chain
                    mi = unique(mi,'rows');
                    % mi_buf: class of job in buffer position i (0=empty)
                    mi_buf = [zeros(size(mi,1),min(sum(n),sn.cap(ist))-S(ist)-size(mi(:,1:end-S(ist)),2)), mi(:,1:end-S(ist))];
                    if isempty(mi_buf)
                        mi_buf = zeros(size(mi,1),1);
                    end
                    % mi_srv: class of job running in server i
                    mi_srv = mi(:,max(size(mi,2)-S(ist)+1,1):end);
                    % si: number of class r jobs that are running
                    si =[];
                    for k=1:size(mi_srv,1)
                        %si(k,1:R) = hist(mi_srv(k,:),1:R); % deprecated
                        si(k, 1:R) = histcounts(mi_srv(k, :), 1:(R+1));
                    end
                    %si = unique(si,'rows');
                    for k=1:size(si,1)
                        % determine number of class r jobs running in phase
                        % j in server state mi_srv(k,:) and build
                        % state
                        kstate=[];
                        map_cols = [];

                        for r=1:R
                            init_r = State.spaceClosedSingle(phases(r),si(k,r));
                            if sn.procid(sn.nodeToStation(ind),r) == ProcessType.MAP || sn.procid(sn.nodeToStation(ind),r) == ProcessType.MMPP2
                                if si(k,r) == 0
                                    init_r = State.cartesian(init_r, [1:phases(r)]');
                                else
                                    if S(ist) == 1
                                        % Single server case (original logic)
                                        init_r = State.cartesian(init_r, 0);
                                        for i=1:size(init_r,1)
                                            if init_r(i,end) == 0
                                                init_r(i,end) = find(init_r(i,:));
                                            end
                                        end
                                    else
                                        % Multiserver case: FCFS-aware phase selection  
                                        % Use original approach but bias toward occupied phases
                                        phase_list = [];
                                        
                                        for i=1:size(init_r,1)
                                            phase_dist = init_r(i, 1:phases(r));
                                            
                                            % Find phases that have jobs (occupied phases)
                                            occupied_phases = find(phase_dist > 0);
                                            
                                            if ~isempty(occupied_phases)
                                                % For FCFS, include occupied phases plus adjacent phases
                                                % to account for phase transitions during service
                                                extended_phases = [];
                                                for op = occupied_phases
                                                    extended_phases = [extended_phases, op];
                                                    % Add adjacent phases for smooth transitions
                                                    if op > 1
                                                        extended_phases = [extended_phases, op-1];
                                                    end
                                                    if op < phases(r)
                                                        extended_phases = [extended_phases, op+1];
                                                    end
                                                end
                                                extended_phases = unique(extended_phases);
                                                phase_list = [phase_list; extended_phases'];
                                            else
                                                % No jobs - include all phases (original behavior)
                                                phase_list = [phase_list; [1:phases(r)]'];
                                            end
                                        end
                                        
                                        % Remove duplicates and create cartesian product
                                        phase_list = unique(phase_list);
                                        init_r = State.cartesian(init_r, phase_list);
                                    end
                                end
                            end
                            kstate = State.cartesian(kstate,init_r);
                            if sn.procid(sn.nodeToStation(ind),r) == ProcessType.MAP || sn.procid(sn.nodeToStation(ind),r) == ProcessType.MMPP2
                                map_cols(end+1) = size(kstate,2);
                            end
                        end
                        kstate = kstate(:,[setdiff(1:size(kstate,2),map_cols),map_cols]);
                        state = [state; repmat(mi_buf(k,:),size(kstate,1),1), kstate];
                    end
                end
                space = state;
            case SchedStrategy.PAS
                % PAS/OI: local state is the full ordered class-index list, left-aligned and zero-padded -- see _kb/04-networkstruct.md
                W = sn.cap(ist);
                if isinf(W)
                    line_error(mfilename,'PAS stations require finite capacity for state-space generation.');
                end
                if sum(n) == 0
                    space = zeros(1, W);
                elseif sum(n) > W
                    space = zeros(0, W); % infeasible: exceeds total capacity
                else
                    vi = [];
                    for r=1:R
                        if n(r)>0
                            vi=[vi, r*ones(1,n(r))];
                        end
                    end
                    mi = multiset_perms(vi);
                    space = [mi, zeros(size(mi,1), W - size(mi,2))];
                end
            case {SchedStrategy.SJF, SchedStrategy.LJF}
                % in these policies the state space includes continuous
                % random variables for the service times
                line_error(mfilename,'The scheduling policy does not admit a discrete state space.\n');
        end
        end % if isRetrialStation
        space = sub_routevars(sn, ind, R, space);
        % The polling controller trails the routing variables, matching the
        % nvars column order (modulation, routing, node block).
        space = State.pollingSpace(sn, ind, space);
        % True BAS blocked marker (nvars col 2*R+1): gate on sn.isbasblocking, not the station's own drop rule (BUG-83) -- see _kb/04-networkstruct.md
        if ~isempty(sn.isbasblocking) && numel(sn.isbasblocking) >= ind ...
                && sn.isbasblocking(ind) == 1 && sum(n) > 0
            space = State.cartesian(space, [0;1]);
        end
        % Server breakdown status (nvars col 2*R+1): enumerated for EVERY marginal incl. empty -- see _kb/04-networkstruct.md
        if isfield(sn,'hasbreakdown') && ~isempty(sn.hasbreakdown) && numel(sn.hasbreakdown) >= ind ...
                && sn.hasbreakdown(ind) == 1
            space = State.cartesian(space, [0;1]);
        end
    case NodeType.Cache
        switch sn.sched(ist)
            case SchedStrategy.INF
                % in this policies we only track the jobs in the servers
                for r=1:R
                    init_r = State.spaceClosedSingle(phases(r),n(r));
                    state = State.cartesian(state,init_r);
                end
                space = State.cartesian(space,state);
        end
        for r=1:R
            switch sn.routing(ind,r)
                case RoutingStrategy.RROBIN
                    space = State.cartesian(space, sn.nodeparam{ind}{r}.outlinks(:));
            end
        end
end
space = unique(space,'rows'); % do not comment, required to sort empty state as first
space = space(end:-1:1,:); % so that states with jobs in phase 1 comes earlier
end

function space = sub_trailingvars(sn, ind, R, space, n)
% Append the shared trailing local-variable column (nvars col 2*R+1) on the
% empty-station early-return paths, which bypass the general appends at the end
% of fromMarginal. Exactly one feature owns that column per station, which
% refreshLocalVars enforces.
if isfield(sn,'hasbreakdown') && ~isempty(sn.hasbreakdown) ...
        && numel(sn.hasbreakdown) >= ind && sn.hasbreakdown(ind) == 1
    % Server status enumerated for the empty station too (unlike BAS/polling), else queue length inflates by ~1 job -- see _kb/04-networkstruct.md
    space = State.cartesian(space, [0;1]);
elseif size(sn.nvars,2) >= 2*R+1 && sn.nvars(ind, 2*R+1) == 1
    % True BAS: keep the empty-station state width consistent (blocked=0); a
    % job can only be held at the server when the station is non-empty.
    space = State.cartesian(space, 0);
end
end

function space = sub_routevars(sn, ind, R, space)
% Append the round-robin routing-variable columns for node IND to SPACE.
% Every branch of fromMarginal must produce rows carrying these columns,
% including the empty-station early returns: an empty state emitted without
% the pointer column misaligns in fromMarginalBounds and silently drops the
% empty configuration from the state space, which inflates the station QLen
% by about one job (the chain can then never empty the station).
for r=1:sn.nclasses
    switch sn.routing(ind,r)
        case RoutingStrategy.RROBIN
            % RR slot holds destination node index -- enumerate over outlinks.
            space = State.cartesian(space, sn.nodeparam{ind}{r}.outlinks(:));
        case RoutingStrategy.WRROBIN
            % WRR slot holds POSITION in weighted_outlinks.
            np = sn.nodeparam{ind}{r};
            if isfield(np, 'weighted_outlinks') && ~isempty(np.weighted_outlinks)
                positions = (1:length(np.weighted_outlinks))';
            else
                positions = (1:length(np.outlinks))';
            end
            space = State.cartesian(space, positions);
    end
end
end
