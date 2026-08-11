function probSysStateAggr = getProbSysAggr(self)
% PROBSYSSTATEAGGR = GETPROBSYSSTATEAGGR()
if GlobalConstants.DummyMode
    probSysStateAggr = NaN;
    return
end
sn = self.getStruct;
TranSysStateAggr = self.sampleSysAggr;
TSS = cell2mat([TranSysStateAggr.t,TranSysStateAggr.state(:)']);
TSS(:,1)=[diff(TSS(:,1));0];
state = sn.state;
nir = zeros(sn.nstateful,sn.nclasses);
for isf=1:sn.nstateful
    ind = sn.statefulToNode(isf);
    [~,nir(isf,:)] = State.toMarginal(sn, ind, state{isf});
end
nir = nir';
rows = findrows(TSS(:,2:end), nir(:)');
if ~isempty(rows)
    probSysStateAggr = sum(TSS(rows,1))/sum(TSS(:,1));
else
    line_warning(mfilename,'The state was not seen during the simulation.\n');
    probSysStateAggr = 0;
end
end