function [outspace, outrate, outprob, eventCache, outstart, outpreempt] = afterEvent(sn, ind, inspace, event, class, isSimulation, eventCache, ctx, noPromote)
% [OUTSPACE, OUTRATE, OUTPROB, EVENTCACHE, OUTSTART, OUTPREEMPT] =  AFTEREVENT(QN, IND, INSPACE, EVENT, CLASS, ISSIMULATION, EVENTCACHE, CTX, NOPROMOTE)
%
% OUTSTART and OUTPREEMPT are (rows of OUTSPACE) x (classes) integer matrices
% holding the START and PREEMPT tags of each successor arc: how many class-r
% jobs begin holding a server on that arc, and how many are pushed back into
% the buffer by it. They are annotations on the arcs above, not events: no
% rate, probability or state depends on them, and a caller that asks for only
% the first four outputs pays nothing for them.
%
% CTX (optional): loop-invariant context precomputed by State.afterEventInit
% on the SAME sn passed here (after any caller-side rewrite of
% nservers/cap/classcap, see solver_ssa preamble). Hot callers pass it to
% skip the per-call setup below; semantics are identical.
%
% NOPROMOTE (optional, default false): forwarded to State.afterEventStation.
% When true, a DEP at an FCFS-family station does not promote a waiting job
% into the vacated server. Set only for the departure half of an
% immediate-feedback self-loop (sn.immfeed) so the fed-back job holds the
% server rather than re-queueing. It is part of the cache key below so that
% immediate-feedback and ordinary departures never collide in the cache.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 9 || isempty(noPromote)
    noPromote = false;
end

% Initialized here so that every exit carries them, exactly as
% State.afterEventStation initializes its own outputs: a node type that
% reaches none of the handlers below has no successor and hence no tags.
outspace = [];
outrate = [];
outprob = [];
outstart = [];
outpreempt = [];

% Event cache for faster simulation
if isSimulation && nargin >= 7 && isobject(eventCache)
    vector = [ind, event, class, double(noPromote), inspace];
    key = mat2str(vector);
    if isKey(eventCache, key)
        cachedResult = eventCache{key};
        outprob = cachedResult{1};
        outspace = cachedResult{2};
        outrate = cachedResult{3};
        % Nodes with no service facility (Router, Fork, Cache, Transition,
        % Join) cache three fields: they can start nothing, so their tags are
        % the zero rows tagPad fills in below.
        if numel(cachedResult) >= 5
            outstart = cachedResult{4};
            outpreempt = cachedResult{5};
        else
            outstart = [];
            outpreempt = [];
        end
        outstart = State.tagPad(outstart, size(outspace,1), sn.nclasses);
        outpreempt = State.tagPad(outpreempt, size(outspace,1), sn.nclasses);
        if size(outspace,1) > 1
            tot_rate = sum(outrate);
            cum_rate = cumsum(outrate) / tot_rate;
            firing_ctr = 1 + max([0,find( rand > cum_rate' )]); % select action
            outspace = outspace(firing_ctr,:);
            outrate = sum(outrate);
            outprob = outprob(firing_ctr,:);
            outstart = outstart(firing_ctr,:);
            outpreempt = outpreempt(firing_ctr,:);
        end
        return
    end
else
    key = NaN;
    eventCache = [];
end
% Join synchronization stations on FJ-augmented structs carry a plain
% per-class count state handled by a dedicated branch (they are stations,
% but must bypass the buffer/phase slicing and afterEventStation below)
if isfield(sn,'isfjaugmented') && sn.isfjaugmented && sn.nodetype(ind) == NodeType.Join
    [outspace, outrate, outprob, eventCache] = State.afterEventJoin(sn, ind, inspace, event, class, isSimulation, eventCache, key);
    % a Join holds no server: nothing starts or is preempted there
    outstart = zeros(size(outspace,1), sn.nclasses);
    outpreempt = zeros(size(outspace,1), sn.nclasses);
    return
end

% Else continue to main body below
usectx = nargin >= 8 && ~isempty(ctx);
if usectx
    % Loop-invariant setup precomputed once by State.afterEventInit
    M = ctx.M;
    R = ctx.R;
    S = ctx.S;
    phasessz = ctx.phasessz;
    phaseshift = ctx.phaseshift;
    pie = ctx.pie;
    isf = sn.nodeToStateful(ind);
    if sn.isstation(ind)
        ismkvmod = ctx.ismkvmod(ind);
        ismkvmodclass = ctx.ismkvmodclass{ind};
    end
    lldscaling = ctx.lldscaling;
    lldlimit = ctx.lldlimit;
    cdscaling = ctx.cdscaling;
else
    M = sn.nstations;
    R = sn.nclasses;
    S = sn.nservers;
    phasessz = sn.phasessz;
    phaseshift = sn.phaseshift;
    pie = sn.pie;

    % ind: node index
    isf = sn.nodeToStateful(ind);

    if sn.isstation(ind)
        ismkvmod = any(sn.procid(sn.nodeToStation(ind),:)==ProcessType.MAP | sn.procid(sn.nodeToStation(ind),:)==ProcessType.MMPP2);
        ismkvmodclass = zeros(R,1);
        for r=1:R
            ismkvmodclass(r) = any(sn.procid(sn.nodeToStation(ind),r)==ProcessType.MAP | sn.procid(sn.nodeToStation(ind),r)==ProcessType.MMPP2);
        end
    end

    lldscaling = sn.lldscaling;
    if isempty(lldscaling)
        lldlimit = max(sum(sn.nclosedjobs),1);
        lldscaling = ones(M,lldlimit);
    else
        lldlimit = size(lldscaling,2);
    end

    % sn.cdscaling carries a class-dependence handle only for the stations that
    % declare one; the others are left empty (see getLimitedClassDependence).
    % afterEventStation indexes every station unconditionally, and an empty entry
    % would be silently INDEXED rather than called -- cdscaling{ist}(nir) on []
    % raises "Array indices must be positive integers" whenever some class has
    % zero jobs. Fill the gaps with the neutral scaling, as afterEventInit does.
    cdscaling = sn.cdscaling;
    if isempty(cdscaling)
        cdscaling = cell(M,1);
    end
    neutral = @(ni) 1;
    for i = 1:M
        if i > numel(cdscaling) || isempty(cdscaling{i})
            cdscaling{i} = neutral;
        end
    end
    % Fold joint-dependence handles (sn.jdscaling, non-product-form eta_i) into
    % the effective per-station handle, exactly as afterEventInit does.
    jdscaling = sn.jdscaling;
    if ~isempty(jdscaling)
        for i = 1:M
            if i <= numel(jdscaling) && ~isempty(jdscaling{i})
                cdh = cdscaling{i};
                jdh = jdscaling{i};
                cdscaling{i} = @(ni) cdh(ni) .* jdh(ni);
            end
        end
    end
end

hasOnlyExp = false; % true if all service processes are exponential
if sn.isstation(ind)
    ist = sn.nodeToStation(ind);
    K = phasessz(ist,:);
    Ks = phaseshift(ist,:);
    if max(K)==1
        hasOnlyExp = true;
    end
    if usectx
        mu = ctx.mu;
        phi = ctx.phi;
        proc = ctx.proc;
        capacity = ctx.capacity;
        classcap = ctx.classcap;
    else
        mu = sn.mu;
        phi = sn.phi;
        proc = sn.proc;
        capacity = sn.cap;
        classcap = sn.classcap;
    end
    if K(class) == 0 % if this class is not accepted at the resource
        eventCache{key} = {outprob, outspace, outrate, outstart, outpreempt};
        return
    end
    V = sum(sn.nvars(ind,:));
    % Place nodes: state format is [buffer(R), server(sum(K))] after ARV
    if sn.nodetype(ind) == NodeType.Place
        space_var = zeros(size(inspace,1), 0);  % proper dimensions for concatenation
        state_len = size(inspace, 2);
        expected_len = R + sum(K);
        if state_len == expected_len
            % State already has [buffer, server] format
            space_buf = inspace(:, 1:R);
            space_srv = inspace(:, (R+1):end);
        elseif state_len == R
            % Initial state: just buffer counts
            space_buf = inspace;
            space_srv = zeros(size(inspace,1), sum(K));
        else
            % Fallback for unexpected formats
            space_buf = inspace;
            space_srv = zeros(size(inspace,1), sum(K));
        end
    else
        if sn.sched(ist) == SchedStrategy.EXT
            % Source state layout from fromMarginal: [Inf_buffer, phases...]
            % plus, when the source uses round-robin routing, trailing local
            % variable columns holding the outlink pointers. MAP arrival
            % variables are not materialized (the phase lives in the service
            % slots), so only the routing variables are sliced off.
            V = sum(sn.nvars(ind,(R+1):(2*R)));
            if V > 0
                space_var = inspace(:,(end-V+1):end);
                space_srv = inspace(:,(end-sum(K)-V+1):(end-V));
                space_buf = inspace(:,1:(end-sum(K)-V));
            else
                space_buf = inspace(:, 1:(size(inspace,2)-sum(K)));
                space_srv = inspace(:, (size(inspace,2)-sum(K)+1):end);
                space_var = zeros(size(inspace,1), 0);
                V = 0;
            end
        elseif sn.sched(ist) == SchedStrategy.PAS
            % PAS/OI layout: [ordered-list(cap) | routing vars(V)]; no server split.
            space_var = inspace(:,(end-V+1):end);
            space_srv = zeros(size(inspace,1), 0);
            space_buf = inspace(:,1:(end-V));
        else
            space_var = inspace(:,(end-V+1):end); % local state variables
            space_srv = inspace(:,(end-sum(K)-V+1):(end-V)); % server state
            space_buf = inspace(:,1:(end-sum(K)-V)); % buffer state
        end
    end
elseif sn.isstateful(ind)
    V = sum(sn.nvars(ind,:));
    % in this case service is always immediate so sum(K)=1
    space_var = inspace(:,(end-V+1):end); % local state variables
    if sn.nodetype(ind) == NodeType.Transition
        K = sn.nodeparam{ind}.firingphases;
        nmodes = sn.nodeparam{ind}.nmodes;
        % Handle NaN firingphases (non-phase-type distributions like Pareto)
        % Infer phase count from D0 matrix size
        if any(isnan(K))
            K = zeros(1, nmodes);
            for m = 1:nmodes
                if iscell(sn.nodeparam{ind}.firingproc) && ~isempty(sn.nodeparam{ind}.firingproc{m})
                    K(m) = size(sn.nodeparam{ind}.firingproc{m}{1}, 1);
                else
                    K(m) = 1;
                end
            end
        end
        Ks = [0,cumsum(K,2)];
        space_buf = inspace(:,1:nmodes); % idle servers count put in buf
        space_srv = inspace(:,(nmodes+1):(nmodes+sum(K))); % enabled servers' phases
        % Handle both state formats: with and without fired component
        expected_len_with_fired = 2*nmodes + sum(K);
        expected_len_without_fired = nmodes + sum(K);
        if size(inspace, 2) >= expected_len_with_fired
            space_fired = inspace(:,(nmodes+sum(K)+1):(2*nmodes+sum(K))); % servers that just fired
        elseif size(inspace, 2) == expected_len_without_fired
            % Legacy format without fired component - initialize to zeros
            space_fired = zeros(size(inspace,1), nmodes);
        else
            line_error(mfilename, 'Unexpected state vector length for Transition node');
        end
    else
        space_buf = []; % buffer state
        space_srv = inspace(:,(end-R-V+1):(end-V)); % server state
        space_fired = []; % only for Transition nodes
    end
else % stateless node
    space_var = [];
    space_srv = [];
    space_buf = [];
end

if sn.isstation(ind)
    [outspace, outrate, outprob, eventCache, outstart, outpreempt] = State.afterEventStation(sn, ind, inspace, event, class, isSimulation, eventCache, ...
        M, R, S, phasessz, phaseshift, pie, isf, ismkvmod, ismkvmodclass, lldscaling, lldlimit, cdscaling, ...
        hasOnlyExp, ist, K, Ks, mu, phi, proc, capacity, classcap, V, space_buf, space_srv, space_var, key, noPromote);
elseif sn.isstateful(ind)
    switch sn.nodetype(ind)
        case NodeType.Router
            [outspace, outrate, outprob, eventCache] = State.afterEventRouter(sn, ind, event, class, isSimulation, eventCache, space_buf, space_srv, space_var, key);
        case NodeType.Fork
            [outspace, outrate, outprob, eventCache] = State.afterEventFork(sn, ind, event, class, isSimulation, eventCache, space_buf, space_srv, space_var, key);
        case NodeType.Cache
            [outspace, outrate, outprob, eventCache] = State.afterEventCache(sn, ind, event, class, isSimulation, eventCache, R, space_buf, space_srv, space_var, key);
        case NodeType.Transition
            [outspace, outrate, outprob, eventCache] = State.afterEventTransition(sn, ind, inspace, K, Ks, event, class, isSimulation, eventCache, R, space_buf, space_srv, space_fired, space_var, key);
    end % switch nodeType
end
% Pad the tail: a tag site only writes the rows it touches, and the nodes above
% that hold no server write nothing at all. From here on the two matrices have
% exactly one row per successor, so a caller can index them with the same row
% indices as outspace.
outstart = State.tagPad(outstart, size(outspace,1), sn.nclasses);
outpreempt = State.tagPad(outpreempt, size(outspace,1), sn.nclasses);
end
