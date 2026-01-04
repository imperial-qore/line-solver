function [state, prior, space] = getState(self)
% [STATE, PRIORSPACE, SPACE] = GETSTATE() % GET STATE

if ~self.hasInitState()
    self.initDefault();
end

nodes = self.nodes;
prior = {};
state = {};
space = {};
for ind=1:length(self.nodes)
    if nodes{ind}.isStateful
        state{end+1,1} = nodes{ind}.getState();
        prior{end+1,1} = nodes{ind}.getStatePrior();
        space{end+1,1} = nodes{ind}.getStateSpace();
    end
end
end