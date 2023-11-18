function initFromMarginalAndStarted(self, n, s, options) % n(i,r) : number of jobs of class r in node i
% INITFROMMARGINALANDSTARTED(N, S, OPTIONS) % N(I,R) : NUMBER OF JOBS OF CLASS R IN NODE I

sn = getStruct(self);
[isvalidn] = State.isValid(sn, n, s);
if ~isvalidn
    line_error(mfilename,'Initial state is not valid.');
end
for ind=1:self.getNumberOfNodes
    if self.nodes{ind}.isStateful
        ist = sn.nodeToStation(ind);
        self.nodes{ind}.setState(State.fromMarginalAndStarted(sn,ind,n(ist,:),s(ist,:)));
        if isempty(self.nodes{ind}.getState)
            line_error(mfilename,sprintf('Invalid state assignment for station %d\n',ind));
        end
    end
end
self.hasState = true;
end
