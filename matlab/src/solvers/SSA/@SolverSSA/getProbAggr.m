function ProbAggr = getProbAggr(self, node, state)
% PROBSTATEAGGR = GETPROBSTATEAGGR(NODE, STATE)

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