function initFromMarginalAndRunning(self, n, s, options) % n(i,r) : number of jobs of class r in node i
% INITFROMMARGINALANDRUNNING(N, S, OPTIONS) % N(I,R) : NUMBER OF JOBS OF CLASS R IN NODE I

sn = getStruct(self);
[isvalidn] = State.isValid(sn, n, s);
if ~isvalidn
    line_error(mfilename,'Initial state is not valid.');
end
for i=1:self.getNumberOfNodes
    if self.nodes{i}.isStateful
        self.nodes{i}.setState(State.fromMarginalAndRunning(sn,i,n(i,:),s(i,:)));
        if isempty(self.nodes{i}.getState)
            line_error(mfilename,sprintf('Invalid state assignment for station %d\n',i));
        end
    end
end
self.hasState = true;
end