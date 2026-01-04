function probSys = getProbSys(self)
% PROBSYS = GETPROBSYS()
if GlobalConstants.DummyMode
    probSys = NaN;
    return
end
TranSysState = self.sampleSys;
TSS = cell2mat([TranSysState.t,TranSysState.state(:)']);
TSS(:,1)=[TSS(1,1);diff(TSS(:,1))];
sn = self.getStruct;
% fill-in for FCFS states
for isf=1:size(TranSysState.state,2)
    if size(sn.state{isf},1)>1
        error('There are multiple station states, choose an initial state as a parameter to getProb.');
    end
    sn.state{isf} = [zeros(1,size(TranSysState.state{isf},2)-size(sn.state{isf},2)),sn.state{isf}];
end
state = cell2mat(sn.state');
rows = findrows(TSS(:,2:end), state);
if ~isempty(rows)
    probSys = sum(TSS(rows,1))/sum(TSS(:,1));
else
    line_warning(mfilename,'The state was not seen during the simulation.\n');
    probSys = 0;
end
end