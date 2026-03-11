function tranSysState = sampleSys(self, numEvents, markActivePassive)
% TRANSYSSTATE = SAMPLESYS(NUMSAMPLES)
options = self.getOptions;

if GlobalConstants.DummyMode
    tranSysState = NaN;
    return
end

if nargin>=2 %exist('numEvents','var')
    options.samples = numEvents;
else
    numEvents = options.samples;
end

if nargin<3
    markActivePassive = false;
end

switch options.method
    case {'default','serial'}
        options.method = 'serial'; % nrm does not support tran*
        [~, tranSystemState, tranSync] = self.runAnalyzer(options);
		
        tranSysState = struct();
        tranSysState.handle = self.model.getStatefulNodes';
        tranSysState.t = tranSystemState{1};
        tranSysState.state = {tranSystemState{2:end}};
        tranSysState.event = tranSync;
        event = tranSysState.event;
        
        for i=1:size(tranSysState.state,2)
            if size(tranSysState.state{i},1) > numEvents
                tranSysState.t = tranSystemState(1:numEvents);
                tranSysState.state = tranSysState.state{i}(1:numEvents,:);
            end
        end
        sn = self.getStruct;
        tranSysState.event = {};
        for e = 1:length(event)
            for a=1:length(sn.sync{event(e)}.active)
                tranSysState.event{end+1} = sn.sync{event(e)}.active{a};
                tranSysState.event{end}.t = tranSysState.t(e);
            end
            for p=1:length(sn.sync{event(e)}.passive)
                tranSysState.event{end+1} = sn.sync{event(e)}.passive{p};
                tranSysState.event{end}.t = tranSysState.t(e);
            end
        end
        tranSysState.isaggregate = false;                
    otherwise
        line_error(mfilename,'sampleSys is not available in SolverSSA with the chosen method.');
end

if markActivePassive
    apevent = cell(1,length(tranSysState.t)-1);
    for ti = 1:length(apevent)
        apevent{ti} = struct('active',[],'passive',[]);
    end
    for e=1:length(tranSysState.event)
        ti = find(tranSysState.event{e}.t == tranSysState.t);
        if ~isempty(ti) && ti<length(tranSysState.t)
            switch tranSysState.event{e}.event
                case EventType.ARV
                    apevent{ti}.passive = tranSysState.event{e};
                otherwise
                    apevent{ti}.active = tranSysState.event{e};
            end
        end
    end
    tranSysState.event = apevent';
end

end