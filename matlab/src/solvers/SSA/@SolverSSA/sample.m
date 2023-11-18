function sampleNodeState = sample(self, node, numSamples, markActivePassive)
% TRANNODESTATE = SAMPLE(NODE)

if GlobalConstants.DummyMode
    sampleNodeState = NaN;
    return
end

options = self.getOptions;

if nargin<4
    markActivePassive = false;
end

if nargin>=3 %exist('numSamples','var')
    options.samples = numSamples;
else
    numSamples = options.samples;
end
switch options.method
    case {'default','serial'}
        [~, tranSystemState, tranSync] = self.runAnalyzer(options);
        event = tranSync;
        sn = self.getStruct;
        isf = sn.nodeToStateful(node.index);
        sampleNodeState = struct();
        sampleNodeState.handle = node;
        sampleNodeState.t = tranSystemState{1};
        sampleNodeState.state = tranSystemState{1+isf};
        
        sn = self.getStruct;
        sampleNodeState.event = {};
        for e = 1:length(event)
            for a=1:length(sn.sync{event(e)}.active)
                sampleNodeState.event{end+1} = sn.sync{event(e)}.active{a};
                sampleNodeState.event{end}.t = sampleNodeState.t(e);
            end
            for p=1:length(sn.sync{event(e)}.passive)
                sampleNodeState.event{end+1} = sn.sync{event(e)}.passive{p};
                sampleNodeState.event{end}.t = sampleNodeState.t(e);
            end
        end
        sampleNodeState.isaggregate = false;

    otherwise
        line_error(mfilename,'sample is not available in SolverSSA with the chosen method.');
end
%sampleNodeState.t = [0; sampleNodeState.t(2:end)];

if markActivePassive
    apevent = cell(1,length(sampleNodeState.t)-1);
    for ti = 1:length(apevent)
        apevent{ti} = struct('active',[],'passive',[]);
    end
    for e=1:length(sampleNodeState.event)
        ti = find(sampleNodeState.event{e}.t == sampleNodeState.t);
        if ~isempty(ti) && ti<length(sampleNodeState.t)
        switch sampleNodeState.event{e}.event
            case EventType.ID_ARV
                apevent{ti}.passive = sampleNodeState.event{e};
            otherwise
                apevent{ti}.active = sampleNodeState.event{e};
        end
        end
    end
    sampleNodeState.event = apevent';
end

end