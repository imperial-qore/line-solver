function initFromMarginalAndStarted(self, n, s, options) % n(i,r) : number of jobs of class r in node i
% INITFROMMARGINALANDSTARTED(N, S, OPTIONS) % N(I,R) : NUMBER OF JOBS OF CLASS R IN NODE I

if nargin<3 %~exist('options','var')
    options = Solver.defaultOptions;
end

if ~self.hasStruct
    self.refreshStruct;
end
sn = getStruct(self);

if self.getNumberOfStations() < self.getNumberOfNodes()
    if size(n,1)==self.getNumberOfStations()
        nnodes = zeros(sn.nnodes,size(n,2));
        snodes = zeros(sn.nnodes,size(s,2));
        for ist=1:sn.nstations
            ind = sn.stationToNode(ist);
            nnodes(ind,:) = n(ist,:);
            snodes(ind,:) = s(ist,:);
        end
        n = nnodes;
        s = snodes;
    elseif size(n,1)==self.getNumberOfNodes()
        % no-op
    else
        line_error(mfilename,'The supplied matrix of marginal states does not have the correct number of rows. One either one per station or one per node.');
    end
end

[isvalidn] = State.isValid(sn, n, s);
if ~isvalidn
    line_error(mfilename,'Initial state is not valid.');
end
for ind=1:self.getNumberOfNodes
    if self.nodes{ind}.isStateful
        ist = sn.nodeToStation(ind);
        state_i = State.fromMarginalAndStarted(sn,ind,n(ist,:),s(ist,:));
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
        if isempty(self.nodes{ind}.getState)
            line_error(mfilename,sprintf('Invalid state assignment for station %d\n',ind));
        end
    end
end
self.hasState = true;
end
