function Prob = getProb(self, node, state)
% PROBSTATE = GETPROBSTATE(NODE, STATE)

if GlobalConstants.DummyMode
    Prob = NaN;
    return
end

% we do not use probSysState as that is for joint states
[~, tranSysState] = self.runAnalyzer;
sn = self.getStruct;
isf = sn.nodeToStateful(node.index);
TSS = cell2mat({tranSysState{1},tranSysState{1+isf}});
TSS(:,1)=[TSS(1,1);diff(TSS(:,1))];
if nargin<3 %~exist('state','var')
    if size(sn.state{isf},1)>1
        error('There are multiple station states, choose an initial state as a parameter to getProb.');
    end
    state = sn.state{isf};
end
% add padding of zeros for FCFS stations
state = [zeros(1,size(TSS(:,2:end),2)-size(state,2)),state];
rows = findrows(TSS(:,2:end), state);
if ~isempty(rows)
    Prob = sum(TSS(rows,1))/sum(TSS(:,1));
else
    line_warning(mfilename,'The state was not seen during the simulation.\n');
    Prob = 0;
end
end