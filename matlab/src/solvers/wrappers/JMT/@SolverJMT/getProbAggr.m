function Pr = getProbAggr(self, node, state_a)
% PR = GETPROBAGGR(NODE, STATE_A)
%
% Probability of a SPECIFIC per-class job distribution at a station.
% Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) for given state.
%
% Compare with getProbMarg: returns total queue-length distribution,
% i.e., P(n total jobs) summed over all class combinations.
%
% Input:
%   node    - Node object
%   state_a - Per-class job counts, e.g., [2,1] = 2 class-1, 1 class-2
%
% Output:
%   Pr - Scalar probability in [0,1] (estimated via simulation)
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