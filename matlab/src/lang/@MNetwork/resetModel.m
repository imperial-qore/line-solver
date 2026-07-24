function resetModel(self, resetState)
% RESETMODEL(RESETSTATE, RESETHANDLES)
%
% If RESETSTATE is true, the model requires re-initialization
% of its state

resetHandles(self);

if self.hasStruct
    if isempty(self.sn)
        rtorig = [];
    else
        rtorig = self.sn.rtorig; % save linked routing table
    end
    self.sn = [];
    self.sn.rtorig = rtorig;
    self.hasStruct = false;
end

if nargin == 2 && resetState
    self.hasState = false;
end
for ind = 1:length(self.getNodes)
    self.nodes{ind}.reset();
end
end
