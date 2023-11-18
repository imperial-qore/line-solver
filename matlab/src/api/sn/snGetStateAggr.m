function [initialStateAggr] = snGetStateAggr(sn) % get initial state
% [INITIALSTATEAGGR] = GETSTATEAGGR() % GET INITIAL STATE

initialState = sn.state;
initialStateAggr = cell(size(initialState));
for isf=1:length(initialStateAggr)
    ind = sn.statefulToNode(isf);
    [~,initialStateAggr{isf}] = State.toMarginalAggr(sn, ind, initialState{isf});
end
end