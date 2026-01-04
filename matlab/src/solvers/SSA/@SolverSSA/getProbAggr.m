function ProbAggr = getProbAggr(self, node, state)
% PROBAGGR = GETPROBAGGR(NODE, STATE)
%
% Probability of a SPECIFIC per-class job distribution at a station.
% Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) for given state.
%
% Compare with getProbMarg: returns total queue-length distribution,
% i.e., P(n total jobs) summed over all class combinations.
%
% Input:
%   node  - Node object
%   state - Per-class job counts, e.g., [2,1] = 2 class-1, 1 class-2
%
% Output:
%   ProbAggr - Scalar probability in [0,1] (estimated via simulation)

if GlobalConstants.DummyMode
    ProbAggr = NaN;
    return
end

% we do not use probSysState as that is for joint states
TranSysStateAggr = self.sampleSysAggr;
sn = self.getStruct;
isf = sn.nodeToStateful(node.index);
TSS = cell2mat({TranSysStateAggr.t,TranSysStateAggr.state{isf}});
TSS(:,1)=[TSS(1,1);diff(TSS(:,1))];
if nargin<3 %~exist('state','var')
    state = sn.state{isf};
end
rows = findrows(TSS(:,2:end), state);
if ~isempty(rows)
    ProbAggr = sum(TSS(rows,1))/sum(TSS(:,1));
else
    line_warning(mfilename,'The state was not seen during the simulation.\n');
    ProbAggr = 0;
end
end