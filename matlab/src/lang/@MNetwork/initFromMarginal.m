function initFromMarginal(self, n, options) % n(i,r) : number of jobs of class r in node or station i (autodetected)
% INITFROMMARGINAL(N, OPTIONS) % N(I,R) : NUMBER OF JOBS OF CLASS R IN NODE OR STATION I

%global GlobalConstants.CoarseTol

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
        for ist=1:sn.nstations
            ind = sn.stationToNode(ist);
            nnodes(ind,:) = n(ist,:);
        end
        n = nnodes;
    elseif size(n,1)==self.getNumberOfNodes()
        % no-op
    else
        line_error(mfilename,'The supplied matrix of marginal states does not have the correct number of rows. One either one per station or one per node.');
    end
end

[isvalidn] = State.isValid(sn, n, [], options);
if ~isvalidn
    %         line_error(mfilename,'The specified state does not have the correct number of jobs.');
    line_warning(mfilename,'Initial state not contained in the state space. Trying to recover.\n');
    n = round(n);
    [isvalidn] = State.isValid(sn, n, [], options);
    if ~isvalidn
        line_error(mfilename,'Cannot recover - stopping.');
    end
end
for ind=1:sn.nnodes
    if sn.isstateful(ind)
        switch sn.nodetype(ind)
            case NodeType.Place
                state_i = n(ind,:);
            case NodeType.Cache
                if max(abs(n(ind,:) - round(n(ind,:)))) < GlobalConstants.CoarseTol
                    state_i = cacheWithContent(self, sn, ind, State.fromMarginal(sn,ind,round(n(ind,:))));
                else % the state argument is purposedly given fractional, e.g. for Fluid solver initialization
                    state_i = n(ind,:);
                end
            otherwise
                if max(abs(n(ind,:) - round(n(ind,:)))) < GlobalConstants.CoarseTol
                    state_i = State.fromMarginal(sn,ind,round(n(ind,:)));
                else % the state argument is purposedly given fractional, e.g. for Fluid solver initialization
                    state_i = n(ind,:);
                end
        end
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
            line_error(mfilename,sprintf('Invalid state assignment for station %d.',ind));
        end
    end
end
self.hasState = true;
end

function state_i = cacheWithContent(self, sn, ind, state_i)
% A marginal fixes only the jobs a Cache holds. Its item content and retrieval
% blocks trail them in the node state, are no function of the marginal, and
% State.fromMarginal does not produce them; dropping them left a state narrower
% than the state space, so every CTMC start from it failed. Carry the content
% the cache holds now, else the initDefault layout; retrieval blocks restart empty.
R = sn.nclasses;
if isfield(sn.nodeparam{ind}, 'totalCacheCapacity')
    tcc = sn.nodeparam{ind}.totalCacheCapacity;
else
    tcc = sn.nvars(ind,2*R+1);
end
rbw = 0; rpw = 0;
if isfield(sn.nodeparam{ind},'retrievalSystemCapacity') ...
        && sn.nodeparam{ind}.retrievalSystemCapacity > 0
    rbw = sn.nodeparam{ind}.nitems;
    rpw = numel(State.cacheRetrievalClassMap(sn, ind));
end
w = size(state_i,2);
cur = self.nodes{ind}.getState();
if size(cur,2) == w + tcc + rbw + rpw
    items = cur(1, w+1:w+tcc);
else
    items = 1:tcc;
end
state_i = [state_i, repmat([items, zeros(1,rbw+rpw)], size(state_i,1), 1)];
end
