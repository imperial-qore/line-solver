function S = sample(self, node, numSamples)
% S = SAMPLE(NODE, NUMSAMPLES)
options = self.getOptions;
options.force = true;
if isempty(self.result) || ~isfield(self.result,'infGen')
    runAnalyzer(self);
end
[infGen, eventFilt] = getGenerator(self);
stateSpace = getStateSpace(self);
%stateSpaceAggr = getStateSpaceAggr(self);

initState = sn.state;
nst = cumsum([1,cellfun(@length,initState)']);
s0 = cell2mat(initState(:)');

% set initial state
pi0 = zeros(1,size(stateSpace,1));
pi0(matchrow(stateSpace,s0))=1;

% filter all CTMC events as a marked Markovian arrival process
D1 = cellsum(eventFilt);
D0 = infGen-D1;
MMAP = mmap_normalize([{D0},{D1},eventFilt(:)']);

% now sampel the MMAP
[sjt,event,~,~,sts] = mmap_sample(MMAP,numSamples, pi0);

sn = self.getStruct;
S = struct();
S.handle = node;
S.t = cumsum([0,sjt(1:end-1)']');
ind = node.index;
isf = sn.nodeToStateful(ind);
S.state = stateSpace(sts,(nst(isf):nst(isf+1)-1));

S.event = {};
%nodeEvent = false(length(event),1);
%nodeTS = zeros(length(event),1);
for e = 1:length(event)
    for a=1:length(sn.sync{event(e)}.active)
        S.event{end+1} = sn.sync{event(e)}.active{a};
        S.event{end}.t = S.t(e);
%        if  sn.sync{event(e)}.active{a}.node == ind 
%            nodeEvent(e) = true;
%            nodeTS(e) = tranSysState.t(e);
%        end
    end
    for p=1:length(sn.sync{event(e)}.passive)
        S.event{end+1} = sn.sync{event(e)}.passive{p};
        S.event{end}.t = S.t(e);
%        if  sn.sync{event(e)}.passive{p}.node == ind 
%            nodeEvent(e) = true;
%            nodeTS(e) = tranSysState.t(e);
%        end
    end
end

%tranSysState.state = tranSysState.state([1;find(nodeTS>0)],:);
%tranSysState.t = unique(nodeTS);
%tranSysState.event = tranSysState.event(nodeEvent)';
S.isaggregate = false;

end