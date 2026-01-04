function stationStateAggr = sampleAggr(self, node, numSamples, markActivePassive)
% SAMPLE = SAMPLEAGGR(NODE, NUMSAMPLES)

if GlobalConstants.DummyMode
    stationStateAggr = NaN;
    return
end

if nargin<2 %~exist('node','var')
    line_error(mfilename,'sampleAggr requires to specify a station.');
end

%if exist('numsamples','var')
    %line_warning(mfilename,'SolveSSA does not support the numsamples parameter, use instead the samples option upon instantiating the solver.');
%end

options = self.getOptions;
switch options.method
    case {'default','serial'}
        options.samples = numSamples;
        options.force = true;
        options.method = 'serial'; % nrm does not support tran*
        [~, tranSystemState, event] = self.runAnalyzer(options);
        sn = self.getStruct;
        isf = sn.nodeToStateful(node.index);
        [~,nir]=State.toMarginal(sn,sn.statefulToNode(isf),tranSystemState{1+isf});
        stationStateAggr = struct();
        stationStateAggr.handle = node;
        stationStateAggr.t = tranSystemState{1};
        stationStateAggr.state = nir;
        sn = self.getStruct;
        stationStateAggr.event = {};
        for e = 1:length(event)
            for a=1:length(sn.sync{event(e)}.active)
                stationStateAggr.event{end+1} = sn.sync{event(e)}.active{a};
                stationStateAggr.event{end}.t = stationStateAggr.t(e);
            end
            for p=1:length(sn.sync{event(e)}.passive)
                stationStateAggr.event{end+1} = sn.sync{event(e)}.passive{p};
                stationStateAggr.event{end}.t = stationStateAggr.t(e);
            end
        end       
        stationStateAggr.isaggregate = true;
    otherwise
        line_error(mfilename,'sampleAggr is not available in SolverSSA with the chosen method.');
end
%stationStateAggr.t = [0; stationStateAggr.t(2:end)];
if markActivePassive
    apevent = cell(1,length(stationStateAggr.t)-1);
    for ti = 1:length(apevent)
        apevent{ti} = struct('active',[],'passive',[]);
    end
    for e=1:length(stationStateAggr.event)
        ti = find(stationStateAggr.event{e}.t == stationStateAggr.t);
        if ~isempty(ti) && ti<length(stationStateAggr.t)
        switch stationStateAggr.event{e}.event
            case EventType.ARV
                apevent{ti}.passive = stationStateAggr.event{e};
            otherwise
                apevent{ti}.active = stationStateAggr.event{e};
        end
        end
    end
    stationStateAggr.event = apevent';
end
end