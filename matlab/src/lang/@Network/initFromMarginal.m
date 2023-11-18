function initFromMarginal(self, n, options) % n(i,r) : number of jobs of class r in node i
% INITFROMMARGINAL(N, OPTIONS) % N(I,R) : NUMBER OF JOBS OF CLASS R IN NODE I

%global GlobalConstants.CoarseTol

if nargin<3 %~exist('options','var')
    options = Solver.defaultOptions;
end

if ~self.hasStruct
    self.refreshStruct;
end
sn = getStruct(self);

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
        ist = sn.nodeToStation(ind);
        switch sn.nodetype(ind)
            case NodeType.Place
                self.nodes{ind}.setState(sum(n(ist,:))); % must be single class token
            otherwise
                if max(abs(n(ist,:) - round(n(ist,:)))) < GlobalConstants.CoarseTol
                    self.nodes{ind}.setState(State.fromMarginal(sn,ind,round(n(ist,:))));
                else % the state argument is purposedly given fractional, e.g. for Fluid solver initialization
                    self.nodes{ind}.setState(n(ist,:));
                end
        end
        if isempty(self.nodes{ind}.getState)
            line_error(mfilename,sprintf('Invalid state assignment for station %d.',ind));
        end
    end
end
self.hasState = true;
end