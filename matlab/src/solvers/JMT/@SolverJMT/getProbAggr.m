function Pr = getProbAggr(self, node, state_a)
% PR = GETPROBSTATEAGGR(NODE, STATE_A)
if GlobalConstants.DummyMode
    Pr = NaN;
    return
end
sn = self.getStruct;
if nargin<3 %~exist('state_a','var')
    state_a = sn.state{sn.nodeToStation(node.index)};
end
stationStateAggr = self.sampleAggr(node);
rows = findrows(stationStateAggr.state, state_a);
t = stationStateAggr.t;
dt = [diff(t);0];
Pr = sum(dt(rows))/sum(dt);
end