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

for ind=nodes
    if sn.isstation(ind)
        ist = sn.nodeToStation(ind);
        n0 = zeros(1,length(N)); % number of jobs in the initial state
        s0 = zeros(1,length(N)); % number of active servers in the initial state
        s = sn.nservers(ist); % allocate
        for r=find(isfinite(N))' % for all closed classes
            if sn.nodeToStation(ind) == sn.refstat(r)
                n0(r) = N(r);
            end
            s0(r) = min(n0(r),s);
            s = s - s0(r);
        end
        switch sn.nodetype(ind)
            case NodeType.Cache
                state_i = State.fromMarginalAndStarted(sn,ind,n0(:)',s0(:)');
                state_i = [state_i, 1:sn.nvars(ind,2*R+1)]; %#ok<AGROW>
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
                state_i = State.fromMarginalAndStarted(sn,ind,n0(:)',s0(:)');
                if sn.isstation(ind)
                    for r=1:sn.nclasses
                        switch sn.procid(sn.nodeToStation(ind),r)
                            case ProcessType.MAP
                                %state_i = State.cartesian(state_i, [1:sn.phases(i,r)]');
                                state_i = State.cartesian(state_i, 1);
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
        if isempty(state_i)
            line_error(mfilename,sprintf('Default initialization failed on station %d.',ind));
        end
    elseif sn.isstateful(ind) % not a station
        switch sn.nodetype(ind)
            case NodeType.Cache
                state_i = [zeros(1,self.getNumberOfClasses), 1:sum(self.nodes{ind}.itemLevelCap)];
            case NodeType.Router
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