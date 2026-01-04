function [outspace, outrate, outprob, eventCache] = afterEvent(sn, ind, inspace, event, class, isSimulation, eventCache)
% [OUTSPACE, OUTRATE, OUTPROB] =  AFTEREVENT(QN, IND, INSPACE, EVENT, CLASS, ISSIMULATION, EVENTCACHE)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Event cache for faster simulation
if isSimulation && nargin >= 7 && isobject(eventCache)
    vector = [ind, event, class, inspace];
    key = mat2str(vector);
    if eventCache.isKey(key)
        cachedResult = eventCache(key);
        outprob = cachedResult{1};
        outspace = cachedResult{2};
        outrate = cachedResult{3};
        if size(outspace,1) > 1
            tot_rate = sum(outrate);
            cum_rate = cumsum(outrate) / tot_rate;
            firing_ctr = 1 + max([0,find( rand > cum_rate' )]); % select action
            outspace = outspace(firing_ctr,:);
            outrate = sum(outrate);
            outprob = outprob(firing_ctr,:);
        end
        return
    end
else
    key = NaN;
    eventCache = [];
end
% Else continue to main body below
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

cdscaling = sn.cdscaling;
if isempty(cdscaling)
    cdscaling = cell(M,1);
    cdscaling{1} = @(ni) 1;
    for i=2:M % faster this way
        cdscaling{i} = cdscaling{1};
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
    mu = sn.mu;
    phi = sn.phi;
    proc = sn.proc;
    capacity = sn.cap;
    classcap = sn.classcap;
    if K(class) == 0 % if this class is not accepted at the resource
        eventCache(key) = {outprob, outspace,outrate};
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
        space_var = inspace(:,(end-V+1):end); % local state variables
        space_srv = inspace(:,(end-sum(K)-V+1):(end-V)); % server state
        space_buf = inspace(:,1:(end-sum(K)-V)); % buffer state
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
    [outspace, outrate, outprob, eventCache] = State.afterEventStation(sn, ind, inspace, event, class, isSimulation, eventCache, ...
        M, R, S, phasessz, phaseshift, pie, isf, ismkvmod, ismkvmodclass, lldscaling, lldlimit, cdscaling, ...
        hasOnlyExp, ist, K, Ks, mu, phi, proc, capacity, classcap, V, space_buf, space_srv, space_var, key);
elseif sn.isstateful(ind)
    switch sn.nodetype(ind)
        case NodeType.Router
            [outspace, outrate, outprob, eventCache] = State.afterEventRouter(sn, ind, event, class, isSimulation, eventCache, space_buf, space_srv, space_var, key);
        case NodeType.Cache
            [outspace, outrate, outprob, eventCache] = State.afterEventCache(sn, ind, event, class, isSimulation, eventCache, R, space_buf, space_srv, space_var, key);
        case NodeType.Transition
            [outspace, outrate, outprob, eventCache] = State.afterEventTransition(sn, ind, inspace, K, Ks, event, class, isSimulation, eventCache, R, space_buf, space_srv, space_fired, space_var, key);
    end % switch nodeType
end
end
