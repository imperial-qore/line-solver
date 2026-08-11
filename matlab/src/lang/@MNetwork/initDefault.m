function initDefault(self, nodes)
% INITDEFAULT(NODES)

% open classes empty
% closed classes initialized at ref station
% running jobs are allocated in class id order until all
% servers are busy

%refreshStruct(self);  % we force update of the model before we initialize

sn = self.getStruct(false);

R = sn.nclasses;
N = sn.njobs';
if nargin < 2
    nodes = 1:self.getNumberOfNodes;
end

% see _kb/04-networkstruct.md (initDefault.m/spaceGenerator.m) for rationale
nplace = zeros(sn.nstations, R);
totplace = zeros(sn.nstations, 1);
for r=find(isfinite(N))'
    refist = sn.refstat(r);
    if sn.nodetype(sn.stationToNode(refist)) == NodeType.Place
        nplace(refist,r) = N(r);
        totplace(refist) = totplace(refist) + N(r);
        continue
    end
    remaining = N(r);
    for jst=[refist, setdiff(1:sn.nstations, refist)]
        if remaining == 0
            break
        end
        if sn.sched(jst) == SchedStrategy.EXT || sn.nodetype(sn.stationToNode(jst)) == NodeType.Place
            continue
        end
        avail = min(sn.classcap(jst,r) - nplace(jst,r), sn.cap(jst) - totplace(jst));
        take = min(remaining, max(0, avail));
        nplace(jst,r) = nplace(jst,r) + take;
        totplace(jst) = totplace(jst) + take;
        remaining = remaining - take;
    end
    if remaining > 0
        line_error(mfilename, sprintf('Cannot place the population of class %d: total station capacity is insufficient.', r));
    end
end

for ind=nodes
    if sn.isstation(ind)
        ist = sn.nodeToStation(ind);
        n0 = nplace(ist,:); % number of jobs in the initial state
        s0 = zeros(1,length(N)); % number of active servers in the initial state
        s = sn.nservers(ist); % allocate
        for r=find(isfinite(N))' % for all closed classes
            s0(r) = min(n0(r),s);
            s = s - s0(r);
        end
        switch sn.nodetype(ind)
            case NodeType.Cache
                state_i = State.fromMarginalAndStarted(sn,ind,n0(:)',s0(:)');
                % Cache state width = totalCacheCapacity + (per-item retrieval bitmap).
                % Initialize the cache region with items 1..totalCacheCapacity and the
                % retrieval-system region (one column per item) with zeros (nothing
                % being retrieved). The bitmap is omitted when no retrieval system.
                if isfield(sn.nodeparam{ind}, 'totalCacheCapacity')
                    tcc = sn.nodeparam{ind}.totalCacheCapacity;
                else
                    tcc = sn.nvars(ind,2*R+1);
                end
                rbw = 0;
                if isfield(sn.nodeparam{ind},'retrievalSystemCapacity') ...
                        && sn.nodeparam{ind}.retrievalSystemCapacity > 0
                    rbw = sn.nodeparam{ind}.nitems;
                end
                state_i = [state_i, 1:tcc, zeros(1,rbw)]; %#ok<AGROW>
            case NodeType.Place
                if sum(self.nodes{ind}.state)>0
                    % if the user pre-loaded manually some jobs, keep them
                    state_i = self.nodes{ind}.state;
                else
                    state_i = zeros(1,self.getNumberOfClasses);
                    for r=1:sn.nclasses
                        if sn.refstat(r) == ist
                            state_i(r) = sn.njobs(r);
                        else
                            state_i(r) = 0;
                        end
                    end
                end
            otherwise
                if sn.isstation(ind) && sn.sched(ist) == SchedStrategy.PAS
                    % see _kb/04-networkstruct.md (initDefault.m/spaceGenerator.m) for rationale
                    userState = self.nodes{ind}.getState();
                    hasUser = ~isempty(userState) && any(userState(:) > 0);
                    sg = [];
                    if numel(sn.nodeparam) >= ind && isstruct(sn.nodeparam{ind}) ...
                            && isfield(sn.nodeparam{ind}, 'swapGraph')
                        sg = sn.nodeparam{ind}.swapGraph;
                    end
                    hasSwap = ~isempty(sg) && any(sg(:) ~= 0);
                    isClosed = any(isfinite(N));
                    if hasUser
                        nUser = zeros(1, R);
                        present = userState(userState > 0);
                        for r = 1:R, nUser(r) = sum(present == r); end
                        sU = zeros(1, R); ss = sn.nservers(ist);
                        for r = 1:R, sU(r) = min(nUser(r), ss); ss = ss - sU(r); end
                        space_i = State.fromMarginalAndStarted(sn, ind, nUser, sU);
                        W = size(space_i, 2);
                        urow = zeros(1, W);
                        cols = min(numel(userState), W);
                        urow(1:cols) = userState(1:cols);
                        % keep the user placement as the initial state (row 1)
                        state_i = unique([urow; space_i], 'rows', 'stable');
                    elseif hasSwap && isClosed && sum(n0) > 1
                        line_error(mfilename, sprintf(['A closed pass-and-swap station with a non-empty swapping ' ...
                            'graph requires an explicit initial job placement (station %d). Call setState on the ' ...
                            'station with the ordered class list (oldest first) before solving.'], ind));
                    else
                        state_i = State.fromMarginalAndStarted(sn, ind, n0(:)', s0(:)');
                    end
                else
                    state_i = State.fromMarginalAndStarted(sn,ind,n0(:)',s0(:)');
                    if sn.isstation(ind)
                        for r=1:sn.nclasses
                            switch sn.procid(sn.nodeToStation(ind),r)
                                case {ProcessType.MAP, ProcessType.MMPP2}
                                    % Markov-modulated service: append the phase-restart
                                    % slot tracked in sn.nvars (see refreshLocalVars)
                                    %state_i = State.cartesian(state_i, [1:sn.phases(i,r)]');
                                    state_i = State.cartesian(state_i, 1);
                            end
                        end
                    end
                end
        end
        for r=1:sn.nclasses
            switch sn.routing(ind,r)
                case {RoutingStrategy.RROBIN, RoutingStrategy.WRROBIN}
                    % start from first connected queue
                    state_i = [state_i, find(sn.connmatrix(ind,:),1)];
            end
        end
        if sn.sched(ist) == SchedStrategy.POLLING
            % The polling controller trails the routing variables, matching the
            % nvars column order used by State.fromMarginal.
            srvclass0 = find(s0 > 0, 1);
            if isempty(srvclass0)
                srvclass0 = 0;
            end
            nbuf0 = n0;
            if srvclass0 > 0
                nbuf0(srvclass0) = nbuf0(srvclass0) - 1;
            end
            state_i = [state_i, State.pollingInit(sn, ind, nbuf0, srvclass0)];
        end
        if isempty(state_i)
            line_error(mfilename,sprintf('Default initialization failed on station %d.',ind));
        end
    elseif sn.isstateful(ind) % not a station
        switch sn.nodetype(ind)
            case NodeType.Cache
                % [class counts | cache contents (items 1..tcc) | per-item retrieval
                % bitmap (zeros, nothing being retrieved)]. The bitmap (one column per
                % item) is omitted when no retrieval system is configured.
                tcc = self.nodes{ind}.totalCacheCapacity;
                if self.nodes{ind}.retrievalSystemCapacity > 0
                    rbw = self.nodes{ind}.items.nitems;
                else
                    rbw = 0;
                end
                state_i = [zeros(1,self.getNumberOfClasses), 1:tcc, zeros(1,rbw)];
            case NodeType.Router
                state_i = zeros(1, self.getNumberOfClasses);
                for r=1:sn.nclasses
                    switch sn.routing(ind,r)
                        case RoutingStrategy.RROBIN
                            % RR slot holds destination value; start at first connected queue
                            state_i = [state_i, find(sn.connmatrix(ind,:),1)]; %#ok<AGROW>
                        case RoutingStrategy.WRROBIN
                            % WRR slot holds POSITION in weighted_outlinks; start at 1
                            state_i = [state_i, 1]; %#ok<AGROW>
                    end
                end
            case NodeType.Fork
                % Stateful Fork (FJ tag-augmented copies only): per-class
                % count of parent jobs held before the fork firing
                state_i = zeros(1, self.getNumberOfClasses);
            case NodeType.Transition
                % Differently from a server, the first nmodes states in a
                % transitions count the servers that are not enabled and
                % the last nodes count servers that just fired.
                % This is required as the local state does not encode the
                % buffer hence the enabling condition is not available.
                state_i = sn.nodeparam{ind}.nmodeservers;
                % For infinite servers, use large finite value for state (SSA needs finite states)
                state_i(isinf(state_i)) = GlobalConstants.MaxInt();
                % For non-Markovian distributions, firingphases is NaN - treat as 1 phase
                firingphases = sn.nodeparam{ind}.firingphases;
                firingphases(isnan(firingphases)) = 1;
                state_i = [state_i, zeros(1, sum(firingphases)), zeros(size(state_i))]; %#ok<AGROW>
            otherwise
                state_i = [];
        end
        %line_error(mfilename,'Default initialization not available on stateful node %d.',i);
    end

    if sn.isstateful(ind) % not a station
        if size(state_i,1)==1
            self.nodes{ind}.setStateSpace(state_i);
            self.nodes{ind}.setStatePrior(1);
            self.nodes{ind}.setState(state_i);
        elseif size(state_i,1)>1
            prior_state_i = zeros(1,size(state_i,1)); prior_state_i(1) = 1;
            self.nodes{ind}.setStateSpace(state_i);
            self.nodes{ind}.setStatePrior(prior_state_i);
            self.nodes{ind}.setState(state_i(1,:));
        else
            self.nodes{ind}.setStateSpace([]);
            self.nodes{ind}.setStatePrior([]);
            self.nodes{ind}.setState([]);
        end
    end
end

if self.isStateValid % problem with example_initState_2
    self.hasState = true;
else
    line_error(mfilename,sprintf('Default initialization failed.'));
end
end