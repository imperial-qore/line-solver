function [initialStateAggr] = getStateAggr(self) % get initial state
% [INITIALSTATEAGGR] = GETSTATEAGGR() % GET INITIAL STATE

initialStateAggr = snGetStateAggr(self.getStruct);
end