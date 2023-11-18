function [initialState, priorInitialState] = getState(self) % get initial state
% [INITIALSTATE, PRIORINITIALSTATE] = GETSTATE() % GET INITIAL STATE

if ~self.hasInitState
    self.initDefault;
end
nodes = self.nodes;
initialState = {};
priorInitialState = {};
for ind=1:length(self.nodes)
    if nodes{ind}.isStateful
        initialState{end+1,1} = nodes{ind}.getState();
        priorInitialState{end+1,1} = nodes{ind}.getStatePrior();
    end
end
end