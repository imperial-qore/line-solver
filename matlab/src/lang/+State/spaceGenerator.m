function [SS,SSh,sn,Adj,ST] = spaceGenerator(sn, cutoff, options)
% SPACEGENERATOR Generate complete state space for queueing network analysis
%
% @brief Creates the complete state space including all possible network states
% @param sn Network structure representing the queueing network
% @param cutoff Population cutoff limits for open classes (scalar or matrix)
% @param options Optional configuration structure for state generation
% @return SS Complete state space matrix
% @return SSh Hashed state space for efficient lookups
% @return sn Updated network structure with state space information
% @return Adj Adjacency matrix for state transitions (SPN support)
% @return ST State transition information (SPN support)
%
% The state space generator creates all possible states for the queueing
% network, including those not reachable from the initial state. This is
% essential for steady-state analysis and performance metric computation.
% For open classes, a cutoff parameter limits the maximum population.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
N = sn.njobs';
Np = N;

% Draft SPN support
Adj = [];
ST = [];

%%
if ~exist('cutoff','var') && any(isinf(Np)) % if the model has open classes
    line_error(mfilename,'Unspecified cutoff for open classes in state space generator.');
end

if isscalar(cutoff)
    cutoff = cutoff * ones(sn.nstations, sn.nclasses);
end

[~, sn, capacityc] = State.spaceGeneratorNodes(sn, cutoff, options);

%%
isOpenClass = isinf(Np);
isClosedClass = ~isOpenClass;
for r=1:sn.nclasses %cut-off open classes to finite capacity
    if isOpenClass(r)
        Np(r) = max(capacityc(:,r)); % if replaced by sum stateMarg_i can exceed capacity
    end
end

nstatefulp = sn.nstateful - sum(sn.nodetype == NodeType.Source); % M without sources

% see _kb/04-networkstruct.md (initDefault.m/spaceGenerator.m) for rationale
if isfield(options,'ctmc_max_states') && ~isempty(options.ctmc_max_states)
    maxStates = options.ctmc_max_states;
else
    maxStates = 3e6;
end

% THE POSITION ENUMERATION IS CAPACITY-BOUND, not just the per-node spaces.
% capacityc already zeroes a (node, class) pair the class never visits, and
% spaceGeneratorNodes uses it to prune sn.space -- but the lattice below used to
% distribute every class over every slot with no bound and let the
% `any(stateMarg_i > capacityc)` test further down reject the impossible rows one
% at a time, at one State.fromMarginal call each. On the class-switching chain
% Source->Q1(A)->Q2(B)->Q3(C) that is (cutoff+1)^(3*3) candidates for
% (cutoff+1)^3 reachable states. slotcap is capacityc in the SLOT order the
% position vector uses, which is the non-source stateful nodes in stateful order.
slotcap = zeros(nstatefulp, sn.nclasses);
nprecsrc = 0;
slotok = true;
for ind = 1:sn.nnodes
    if sn.nodetype(ind) == NodeType.Source
        nprecsrc = nprecsrc + 1;
        continue
    end
    if ~sn.isstateful(ind)
        continue
    end
    slot = sn.nodeToStateful(ind) - nprecsrc;
    if slot < 1 || slot > nstatefulp
        slotok = false;
        break
    end
    slotcap(slot,:) = capacityc(ind,:);
end
if ~slotok
    slotcap = [];   % a slot went unmapped: do not bound what we cannot see
end

n = pprod(Np);
chainStationPos=[];
while n>=0
    % Cooperative wall-clock budget checkpoint (options.timeout): state-space
    % enumeration can dominate the solve time, and a partial space would give
    % silently wrong results, so exceeding the budget must abort the solve.
    if nargin>=3 && lineTimeoutExceeded(options)
        line_error(mfilename,'State space generation exceeded the wall-clock time budget (options.timeout=%gs).', options.timeout);
    end
    % this is rather inefficient since n ignores chains and thus it can
    % generate state spaces such as
    %   J =
    %
    %      0     0     0     0
    %      0     0     0     1
    %      0     0     1     0
    %      0     1     0     0
    %      1     0     0     0
    %      0     0     0     1
    %      0     0     1     0
    %      0     1     0     0
    %      1     0     0     0
    %      0     0     0     2
    %      0     0     1     1
    %      0     0     2     0
    %      0     1     0     1
    %      1     0     0     1
    %      0     1     1     0
    %      1     0     1     0
    %      0     2     0     0
    %      1     1     0     0
    %      2     0     0     0
    % that are then in the need for a call to unique
    if all(isOpenClass) | (Np(isClosedClass) == n(isClosedClass)) %#ok<OR2>
        chainStationPos = [chainStationPos; State.spaceClosedMultiCS(nstatefulp,n,sn.chains,slotcap)];
        if size(chainStationPos,1) > maxStates
            line_error(mfilename,'State space too large: population lattice exceeds ctmc_max_states=%g. Increase options.ctmc_max_states or use a different solver.', maxStates);
        end
    end
    n = pprod(n,Np);
end

chainStationPos = unique(chainStationPos,'rows');

netstates = cell(size(chainStationPos,1), sn.nstateful);
for j=1:size(chainStationPos,1)
    % Cooperative wall-clock budget checkpoint: the per-marginal local state
    % enumeration below (State.fromMarginal per node) dominates for large
    % buffers, so the budget must also be enforced at this granularity.
    if nargin>=3 && lineTimeoutExceeded(options)
        line_error(mfilename,'State space generation exceeded the wall-clock time budget (options.timeout=%gs).', options.timeout);
    end
    for ind=1:sn.nnodes
        if sn.nodetype(ind) == NodeType.Source
            isf = sn.nodeToStateful(ind);
            state_i = State.fromMarginal(sn,ind,[]);
            netstates{j,isf} = State.getHash(sn,ind,state_i);
        elseif sn.isstation(ind)
            isf = sn.nodeToStateful(ind);
            stateMarg_i = chainStationPos(j,(isf-sum(sn.nodetype(1:ind-1) == NodeType.Source)):nstatefulp:end);
            if any(stateMarg_i > capacityc(ind,:))
                netstates{j,isf} = State.getHash(sn,ind,[]);
            else
                state_i = State.fromMarginal(sn,ind,stateMarg_i);
                netstates{j,isf} = State.getHash(sn,ind,state_i);
            end
        elseif sn.isstateful(ind)
            isf = sn.nodeToStateful(ind);
            stateMarg_i = chainStationPos(j,(isf-sum(sn.nodetype(1:ind-1) == NodeType.Source)):nstatefulp:end);
            state_i = sn.space{isf};
            if any(stateMarg_i > capacityc(ind,:))
                netstates{j,isf} = State.getHash(sn,ind,[]);
            elseif sn.nodetype(ind) == NodeType.Cache
                % For cache nodes, we need to handle states with jobs in multiple classes
                % (InitClass, HitClass, MissClass) since cache nodes support class switching
                state_i = state_i(findrows(state_i(:,1:length(stateMarg_i)),stateMarg_i),:);
                np = sn.nodeparam{ind};
                if isfield(np,'retrievalSystemCapacity') && np.retrievalSystemCapacity > 0 ...
                        && isfield(np,'retrievalClasses') && ~isempty(np.retrievalClasses)
                    % see _kb/04-networkstruct.md (initDefault.m/spaceGenerator.m) for rationale
                    rc  = np.retrievalClasses;
                    tcc = np.totalCacheCapacity;
                    lvs = length(stateMarg_i);     % per-class server presence width
                    rsqi = np.retrievalSystemQueueIndices;
                    itemsInQueue = [];             % 1-based item ids currently at a queue
                    validMarginal = true;
                    if ~isempty(rsqi)
                        aks = keys(rsqi);
                        for ak = 1:numel(aks)
                            arrivalClass = aks(ak);            % int32, 0-based arrival class
                            classCol = double(arrivalClass) + 1;
                            qNodes = rsqi{arrivalClass};
                            for qq = 1:numel(qNodes)
                                qIdx = qNodes(qq);
                                qOff = sn.nodeToStateful(qIdx) - sum(sn.nodetype(1:qIdx-1) == NodeType.Source);
                                for item = 1:size(rc,1)
                                    if classCol > size(rc,2), continue; end
                                    rClass = rc(item, classCol);
                                    if rClass < 1, continue; end
                                    col = qOff + (rClass-1)*nstatefulp;
                                    if col > size(chainStationPos,2) || chainStationPos(j,col) == 0, continue; end
                                    % item at a queue: no retrieval job may be at the cache server,
                                    % and the item can be at only one queue
                                    if stateMarg_i(rClass) > 0 || any(itemsInQueue == item)
                                        validMarginal = false; break;
                                    end
                                    itemsInQueue(end+1) = item; %#ok<AGROW>
                                end
                                if ~validMarginal, break; end
                            end
                            if ~validMarginal, break; end
                        end
                    end
                    if ~validMarginal
                        state_i = zeros(0, size(state_i,2));
                    else
                        keep = false(size(state_i,1),1);
                        for row = 1:size(state_i,1)
                            st = state_i(row,:);
                            validRow = true;
                            itemsInRSstate = [];
                            % only block A (one column per item) records in-flight
                            % fetches; block B holds the merged secondary requests
                            for col = (lvs+tcc+1):(lvs+tcc+size(rc,1))
                                if st(col) == 0, continue; end
                                item = col - (lvs+tcc);
                                for jac = 1:size(rc,2)
                                    rClass = rc(item, jac);
                                    if rClass < 1, continue; end
                                    % item in the bitmap but not at a queue requires a
                                    % retrieval job ready to depart/arrive at the cache server
                                    if ~any(itemsInQueue == item) && st(rClass) == 0
                                        validRow = false; break;
                                    end
                                end
                                if ~validRow, break; end
                                itemsInRSstate(end+1) = item; %#ok<AGROW>
                            end
                            if validRow
                                % every item at a queue must be recorded in the bitmap
                                for it = itemsInQueue
                                    if ~any(itemsInRSstate == it)
                                        validRow = false; break;
                                    end
                                end
                            end
                            keep(row) = validRow;
                        end
                        state_i = state_i(keep,:);
                    end
                end
                netstates{j,isf} = State.getHash(sn,ind,state_i);
            elseif sn.nodetype(ind) == NodeType.Transition
                % Transition states are per-mode, not per-class; include all states
                netstates{j,isf} = State.getHash(sn,ind,state_i);
            else
                state_i = state_i(findrows(state_i(:,1:length(stateMarg_i)),stateMarg_i),:);
                netstates{j,isf} = State.getHash(sn,ind,state_i);
            end
        end
    end
end

ctr = 0;
%SS = sparse([]);
SS = [];
SSh = [];
tochk = 0;
for j=1:size(chainStationPos,1)
    % for each network state
    v = {netstates{j,:}};
    % cycle over lattice
    vN = cellfun(@length,v)-1;
    n = pprod(vN);
    while n >=0
        % Cooperative wall-clock budget checkpoint: this cartesian composition
        % over the per-node state lattice is the combinatorial hot loop, so the
        % budget is checked here too (amortized every 8192 iterations).
        tochk = tochk + 1;
        if tochk >= 1024
            tochk = 0;
            if nargin>=3 && lineTimeoutExceeded(options)
                line_error(mfilename,'State space generation exceeded the wall-clock time budget (options.timeout=%gs).', options.timeout);
            end
        end
        u={}; h={};
        skip = false;
        for isf=1:length(n)
            h{isf} = v{isf}(1+n(isf));
            if h{isf} < 0
                skip = true;
                break
            end
            u{isf} = sn.space{isf}(v{isf}(1+n(isf)),:);
        end
        if skip == false
            ctr = ctr + 1; % do not move
            if ctr > maxStates
                line_error(mfilename,'State space too large: composed states exceed ctmc_max_states=%g. Increase options.ctmc_max_states or use a different solver.', maxStates);
            end
            SS(ctr,:)=cell2mat(u);
            SSh(ctr,:)=cell2mat(h);
        end
        n = pprod(n,vN);
    end
end
[SS,IA] = unique(SS,'rows');
SSh = SSh(IA,:);
end
