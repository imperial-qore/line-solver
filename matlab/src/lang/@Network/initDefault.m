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
        n0 = zeros(1,length(N));
        s0 = zeros(1,length(N));
        s = sn.nservers(sn.nodeToStation(ind)); % allocate
        for r=find(isfinite(N))' % for all closed classes
            if sn.nodeToStation(ind) == sn.refstat(r)
                n0(r) = N(r);
            end
            s0(r) = min(n0(r),s);
            s = s - s0(r);
        end
        state_i = State.fromMarginalAndStarted(sn,ind,n0(:)',s0(:)');
        switch sn.nodetype(ind)
            case NodeType.Cache
                state_i = [state_i, 1:sn.nvars(ind,2*R+1)];
            case NodeType.Place
                state_i = 0; % for now PNs are single class
            otherwise
                if sn.isstation(ind)
                    for r=1:sn.nclasses
                        switch sn.procid(sn.nodeToStation(ind),r)
                            case ProcessType.ID_MAP
                                %state_i = State.decorate(state_i, [1:sn.phases(i,r)]');
                                state_i = State.decorate(state_i, 1);
                        end
                    end
                end
        end
        for r=1:sn.nclasses
            switch sn.routing(ind,r)
                case {RoutingStrategy.ID_RROBIN, RoutingStrategy.ID_WRROBIN}
                    % start from first connected queue
                    state_i = [state_i, find(sn.rt(ind,:),1)];
            end
        end
        if isempty(state_i)
            line_error(mfilename,sprintf('Default initialization failed on station %d.',ind));
        else
            %state_i = state_i(1,:); % to change: this effectively disables priors
            self.nodes{ind}.setState(state_i);
            prior_state_i = zeros(1,size(state_i,1)); prior_state_i(1) = 1;
            self.nodes{ind}.setStatePrior(prior_state_i);
        end
    elseif sn.isstateful(ind) % not a station
        switch class(self.nodes{ind})
            case 'Cache'
                state_i = zeros(1,self.getNumberOfClasses);
                state_i = [state_i, 1:sum(self.nodes{ind}.itemLevelCap)];
                self.nodes{ind}.setState(state_i);
            case 'Router'
                self.nodes{ind}.setState([0,1]);
            otherwise
                self.nodes{ind}.setState([]);
        end
        %line_error(mfilename,'Default initialization not available on stateful node %d.',i);
    end
end

if self.isStateValid % problem with example_initState_2
    self.hasState = true;
else
    line_error(mfilename,sprintf('Default initialization failed.'));
end
end